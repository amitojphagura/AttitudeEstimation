/*
  ORIENT testbed  -  Mahony attitude filter on ESP32 + MPU6050
  Orientation Reference Inertial Estimation Testbed

  Structurally identical to esp32_madgwick_fixed.ino: same MPU6050 driver, same
  startup gyro bias calibration, same 100 Hz scheduler, same web server, same
  reset behaviour, same units end to end. Only the filter object changes, so any
  difference you measure between the two sketches is the filter and nothing else.

  Library needed beyond the Madgwick sketch's set:
    Mahony            (Library Manager, search "Mahony", by Arduino / Paul Stoffregen)

  Units, and they are the same trap as MadgwickAHRS:
    updateIMU() wants gyro in DEGREES per second. It converts to rad/s itself.
    getRoll()/getPitch()/getYaw() return DEGREES. Three.js wants radians.

  One Mahony-specific quirk: the library's getYaw() returns yaw + 180 degrees,
  which MadgwickAHRS does not do. yawOffset is primed with 180 below so this
  sketch starts near zero like the others.
*/

#include <Arduino.h>
#include <WiFi.h>
#include <AsyncTCP.h>
#include <ESPAsyncWebServer.h>
#include <MPU6050.h>
#include <Wire.h>
#include <Arduino_JSON.h>
#include <MahonyAHRS.h>
#include "LittleFS.h"

// Replace with your network credentials
const char* ssid = "";
const char* password = "";

AsyncWebServer server(80);
AsyncEventSource events("/events");

JSONVar readings;

// Timing
unsigned long lastTimeTemperature = 0;
unsigned long lastTimeAcc = 0;
unsigned long lastDriftCheck = 0;
unsigned long temperatureDelay = 1000;
unsigned long accelerometerDelay = 200;
unsigned long driftCheckInterval = 5000;

// Filter update scheduler (target 100 Hz)
const unsigned long filterInterval = 10000;  // microseconds
unsigned long lastFilterMicros = 0;
unsigned long lastUpdateMicros = 0;

MPU6050 mpu;
Mahony filter;

int16_t ax, ay, az, gx, gy, gz;
float accX, accY, accZ;

// Gyro bias (raw counts), set at startup
float biasGx = 0, biasGy = 0, biasGz = 0;

// Latest filter output in RADIANS (after offset)
float rollRad = 0, pitchRad = 0, yawRad = 0;

// Reset offsets in DEGREES (captured from filter output)
float rollOffset = 0, pitchOffset = 0, yawOffset = 0;

const float DEG_TO_RADF = 0.0174533;

void calibrateGyro() {
  Serial.println("Calibrating gyro (keep board still for 2 seconds)...");
  delay(2000);

  float sumGx = 0, sumGy = 0, sumGz = 0;
  int calibSamples = 200;
  for (int i = 0; i < calibSamples; i++) {
    mpu.getMotion6(&ax, &ay, &az, &gx, &gy, &gz);
    sumGx += gx;
    sumGy += gy;
    sumGz += gz;
    delay(10);
  }
  biasGx = sumGx / calibSamples;
  biasGy = sumGy / calibSamples;
  biasGz = sumGz / calibSamples;

  Serial.print("Gyro bias - X: ");
  Serial.print(biasGx); Serial.print(" Y: ");
  Serial.print(biasGy); Serial.print(" Z: ");
  Serial.println(biasGz);
}

void initMPU() {
  Wire.begin(21, 22);
  mpu.initialize();
  Serial.println("MPU6050 Found!");

  calibrateGyro();

  filter.begin(100);  // initial guess; real dt is fed each cycle below

  // The Mahony library reports yaw with a fixed +180 degree offset. Prime the
  // reset offset with it so the cube starts near zero rather than backwards.
  yawOffset = 180.0;

  lastUpdateMicros = micros();
  Serial.println("Mahony filter initialized");
}

void initLittleFS() {
  if (!LittleFS.begin()) {
    Serial.println("An error has occurred while mounting LittleFS");
  }
  Serial.println("LittleFS mounted successfully");
}

void initWiFi() {
  WiFi.mode(WIFI_STA);
  WiFi.begin(ssid, password);
  Serial.println("");
  Serial.print("Connecting to WiFi...");
  while (WiFi.status() != WL_CONNECTED) {
    Serial.print(".");
    delay(1000);
  }
  Serial.println("");
  Serial.println(WiFi.localIP());
}

// Read sensor, run one Mahony step using the true elapsed time.
void updateFilter() {
  mpu.getMotion6(&ax, &ay, &az, &gx, &gy, &gz);

  // Accelerometer in g
  float accelX = ax / 16384.0;
  float accelY = ay / 16384.0;
  float accelZ = az / 16384.0;

  // Gyro: subtract bias, convert to DEG/S (the library wants deg/s)
  float gyroX = (gx - biasGx) / 131.0;
  float gyroY = (gy - biasGy) / 131.0;
  float gyroZ = (gz - biasGz) / 131.0;

  // Feed the real sample rate so loop jitter does not distort integration
  unsigned long now = micros();
  float dt = (now - lastUpdateMicros) / 1000000.0;
  lastUpdateMicros = now;
  if (dt > 0.0001 && dt < 0.5) {
    filter.begin(1.0 / dt);  // begin() just sets invSampleFreq = dt
  }

  filter.updateIMU(gyroX, gyroY, gyroZ, accelX, accelY, accelZ);

  // Library returns DEGREES; apply reset offset, then convert to radians
  float rollDeg  = filter.getRoll()  - rollOffset;
  float pitchDeg = filter.getPitch() - pitchOffset;
  float yawDeg   = filter.getYaw()   - yawOffset;

  rollRad  = rollDeg  * DEG_TO_RADF;
  pitchRad = pitchDeg * DEG_TO_RADF;
  yawRad   = yawDeg   * DEG_TO_RADF;
}

String getGyroReadings() {
  // Send filter output (radians) in the existing field names so the
  // cube and cards keep working unchanged.
  readings["gyroX"] = String(rollRad, 4);
  readings["gyroY"] = String(pitchRad, 4);
  readings["gyroZ"] = String(yawRad, 4);
  return JSON.stringify(readings);
}

String getAccReadings() {
  accX = ax / 16384.0;
  accY = ay / 16384.0;
  accZ = az / 16384.0;
  readings["accX"] = String(accX, 4);
  readings["accY"] = String(accY, 4);
  readings["accZ"] = String(accZ, 4);
  return JSON.stringify(readings);
}

String getTemperature() {
  int16_t temp = mpu.getTemperature();
  float temperature = (temp / 340.0) + 36.53;
  return String(temperature, 2);
}

void checkDrift() {
  Serial.print("Roll: ");
  Serial.print(rollRad, 3);
  Serial.print(" rad | Pitch: ");
  Serial.print(pitchRad, 3);
  Serial.print(" rad | Yaw: ");
  Serial.print(yawRad, 3);
  Serial.print(" rad   (deg: ");
  Serial.print(rollRad / DEG_TO_RADF, 1);
  Serial.print(", ");
  Serial.print(pitchRad / DEG_TO_RADF, 1);
  Serial.print(", ");
  Serial.print(yawRad / DEG_TO_RADF, 1);
  Serial.println(")");
}

void setup() {
  Serial.begin(115200);
  initMPU();
  initWiFi();
  initLittleFS();

  server.on("/", HTTP_GET, [](AsyncWebServerRequest *request){
    request->send(LittleFS, "/index.html", "text/html");
  });

  server.serveStatic("/", LittleFS, "/");

  // Reset captures current orientation as the new zero
  server.on("/reset", HTTP_GET, [](AsyncWebServerRequest *request){
    rollOffset  = filter.getRoll();
    pitchOffset = filter.getPitch();
    yawOffset   = filter.getYaw();
    request->send(200, "text/plain", "OK");
  });

  server.on("/resetX", HTTP_GET, [](AsyncWebServerRequest *request){
    rollOffset = filter.getRoll();
    request->send(200, "text/plain", "OK");
  });

  server.on("/resetY", HTTP_GET, [](AsyncWebServerRequest *request){
    pitchOffset = filter.getPitch();
    request->send(200, "text/plain", "OK");
  });

  server.on("/resetZ", HTTP_GET, [](AsyncWebServerRequest *request){
    yawOffset = filter.getYaw();
    request->send(200, "text/plain", "OK");
  });

  events.onConnect([](AsyncEventSourceClient *client){
    if(client->lastId()){
      Serial.printf("Client reconnected! Last message ID that it got is: %u\n", client->lastId());
    }
    client->send("hello!", NULL, millis(), 10000);
  });
  server.addHandler(&events);

  server.begin();
  Serial.println("Server started");
}

void loop() {
  // Run the filter at a steady ~100 Hz using a micros scheduler
  if ((micros() - lastFilterMicros) >= filterInterval) {
    lastFilterMicros = micros();
    updateFilter();
    events.send(getGyroReadings().c_str(), "gyro_readings", millis());
  }

  if ((millis() - lastTimeAcc) > accelerometerDelay) {
    events.send(getAccReadings().c_str(), "accelerometer_readings", millis());
    lastTimeAcc = millis();
  }

  if ((millis() - lastTimeTemperature) > temperatureDelay) {
    events.send(getTemperature().c_str(), "temperature_reading", millis());
    lastTimeTemperature = millis();
  }

  if ((millis() - lastDriftCheck) > driftCheckInterval) {
    checkDrift();
    lastDriftCheck = millis();
  }
}
