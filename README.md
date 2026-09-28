# AttitudeEstimation

An ESP32 + MPU6050 testbed for comparing five attitude estimation filters head to head. Each filter runs in an identical sketch (same sensor driver, same startup gyro calibration, same 100 Hz update rate, same web server) so the filter itself is the only variable. Each one streams live orientation to a small web dashboard with a 3D cube, plus accelerometer and temperature readings.

## Circuit

Wire the MPU6050 to the ESP32 over I2C:

| MPU6050 | ESP32 |
|---|---|
| VCC | 3.3V |
| GND | GND |
| SCL | GPIO22 |
| SDA | GPIO21 |

These are the ESP32's default I2C pins, set in the sketch via `Wire.begin(21, 22)`.

## Filters

Each filter lives in its own folder (sketch, filter header if it needs one, and its web dashboard) so you can grab just one or the whole repo.

* **EKF**: 7 state Extended Kalman Filter. Keeps estimating gyro bias while running. Needs the `BasicLinearAlgebra` library.
* **UKF**: Same 7 state model as EKF, but propagated with sigma points instead of Jacobians. Needs the `BasicLinearAlgebra` library.
* **Madgwick**: Classic gradient descent complementary filter. Lightweight, no extra library.
* **Mahony**: Similar complementary filter with a different correction gain structure. Needs the `MahonyAHRS` library.
* **KalmanSimple**: Two independent per axis linear Kalman filters (angle and bias) for roll and pitch. Yaw has no accelerometer reference so it's bias corrected gyro integration only, and it will drift. No extra library.

## How to use

1. Install the Arduino IDE with ESP32 board support.
2. Install these libraries from the Library Manager: `MPU6050` (I2Cdevlib), `ESPAsyncWebServer`, `AsyncTCP`, `Arduino_JSON`, plus whichever filter specific library your chosen folder needs (`MadgwickAHRS`, `MahonyAHRS`, or `BasicLinearAlgebra`).
3. Wire the MPU6050 as shown above.
4. Open the `.ino` file from whichever filter folder you want to try.
5. Set your WiFi `ssid` and `password` near the top of the file.
6. Upload the folder's `data` subfolder to the ESP32's filesystem using a LittleFS/SPIFFS data upload tool for the Arduino IDE.
7. Upload the sketch itself to the ESP32. Open the Serial Monitor at 115200 baud to see the board's IP address once it connects to WiFi.
8. Open that IP address in a browser to see live orientation, accelerometer readings, temperature, and a 3D cube. The reset buttons zero out roll, pitch, and yaw from the current position. 

## Sources
https://randomnerdtutorials.com/esp32-mpu-6050-accelerometer-gyroscope-arduino/

