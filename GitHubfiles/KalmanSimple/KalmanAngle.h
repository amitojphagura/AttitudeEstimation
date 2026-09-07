/*
  KalmanAngle.h  -  the classic per-axis 2-state linear Kalman filter for
  attitude from a 6-axis IMU, plus a thin AHRS wrapper that exposes the same
  getRoll() / getPitch() / getYaw() interface as MadgwickAHRS, Mahony, EKF7 and
  UKF7 so the ORIENT sketches stay interchangeable.

  Part of the ORIENT testbed (Orientation Reference Inertial Estimation Testbed).

  Per axis
  --------
    state        x = [ angle , bias ]^T          (degrees, degrees per second)
    prediction   angle += dt * (gyro - bias)     bias is constant
    measurement  the tilt angle derived from the accelerometer

  This is the linear-Gaussian baseline of the study: no quaternion, no coupling
  between axes, three scalar tunables. It is what "Kalman filter" means in most
  hobby IMU code, and the point of including it is to show what the linear
  small-angle model costs relative to the quaternion EKF and UKF.

  Known limitations, and they are deliberate
  --------------------------------------------
  1. Roll is driven by gyro X and pitch by gyro Y. Body rates only equal Euler
     rates near level, so error grows with tilt.
  2. Roll and pitch are estimated independently, so the filter has no notion of
     the two being coupled through the same rotation.
  3. Yaw has no measurement at all on a 6-axis IMU. It is bias-corrected gyro
     integration and it will drift. That drift is a result to report, not a bug
     to fix.

  Units in:   gyro in DEGREES per second, accelerometer in g.
  Units out:  DEGREES, matching every other ORIENT filter.

  No external library is required. (rfetick's "Kalman" library from the Library
  Manager was considered and rejected: it pins BasicLinearAlgebra 3.x, which
  conflicts with the 5.x that the EKF and UKF sketches need.)
*/

#ifndef ORIENT_KALMAN_ANGLE_H
#define ORIENT_KALMAN_ANGLE_H

#include <math.h>

class KalmanAngle
{
   public:
    // Process noise on the angle state. Larger -> follows the accelerometer
    // more closely and is noisier.
    float Q_angle = 0.001f;

    // Process noise on the bias state. Larger -> tracks a changing gyro bias
    // faster and wanders more.
    float Q_bias = 0.003f;

    // Accelerometer measurement noise. Larger -> leans on the gyro, smoother
    // but slower to correct drift.
    float R_measure = 0.03f;

    void begin(float startAngle = 0.0f)
    {
        angle = startAngle;
        bias = 0.0f;
        P[0][0] = 0.0f;
        P[0][1] = 0.0f;
        P[1][0] = 0.0f;
        P[1][1] = 0.0f;
    }

    // newAngle in degrees (from the accelerometer), newRate in deg/s (gyro),
    // dt in seconds. Returns the filtered angle in degrees.
    float update(float newAngle, float newRate, float dt)
    {
        // --- Predict ---------------------------------------------------
        rate = newRate - bias;
        angle += dt * rate;

        P[0][0] += dt * (dt * P[1][1] - P[0][1] - P[1][0] + Q_angle);
        P[0][1] -= dt * P[1][1];
        P[1][0] -= dt * P[1][1];
        P[1][1] += Q_bias * dt;

        // --- Update ----------------------------------------------------
        float S = P[0][0] + R_measure;
        float K0 = P[0][0] / S;
        float K1 = P[1][0] / S;

        float y = newAngle - angle;
        angle += K0 * y;
        bias += K1 * y;

        float P00 = P[0][0];
        float P01 = P[0][1];
        P[0][0] -= K0 * P00;
        P[0][1] -= K0 * P01;
        P[1][0] -= K1 * P00;
        P[1][1] -= K1 * P01;

        return angle;
    }

    void setAngle(float a) { angle = a; }
    float getAngle() const { return angle; }
    float getRate() const { return rate; }
    float getBias() const { return bias; }

   private:
    float angle = 0.0f;
    float bias = 0.0f;
    float rate = 0.0f;
    float P[2][2] = {{0, 0}, {0, 0}};
};

// ---------------------------------------------------------------------------

class SimpleKalmanAHRS
{
   public:
    KalmanAngle kalmanRoll;
    KalmanAngle kalmanPitch;

    void begin()
    {
        kalmanRoll.begin(0.0f);
        kalmanPitch.begin(0.0f);
        yawDeg = 0.0f;
        seeded = false;
    }

    // Seed roll and pitch from one accelerometer sample so the filter does not
    // have to converge from level. Yaw starts at zero by definition.
    void setFromAccel(float ax, float ay, float az)
    {
        kalmanRoll.begin(accelRoll(ax, ay, az));
        kalmanPitch.begin(accelPitch(ax, ay, az));
        yawDeg = 0.0f;
        seeded = true;
    }

    // gx,gy,gz in deg/s ; ax,ay,az in g ; dt in seconds.
    void update(float gxDps, float gyDps, float gzDps,
                float ax, float ay, float az, float dt)
    {
        if (dt <= 0.0f || dt > 0.5f) return;

        float rollAcc = accelRoll(ax, ay, az);
        float pitchAcc = accelPitch(ax, ay, az);

        if (!seeded)
        {
            kalmanRoll.begin(rollAcc);
            kalmanPitch.begin(pitchAcc);
            seeded = true;
        }

        // The accelerometer roll wraps at +/-180. When it jumps across that
        // seam, restart the filter's angle there instead of chasing a 360
        // degree innovation.
        if ((rollAcc < -90.0f && kalmanRoll.getAngle() > 90.0f) ||
            (rollAcc > 90.0f && kalmanRoll.getAngle() < -90.0f))
        {
            kalmanRoll.setAngle(rollAcc);
            rollDeg = rollAcc;
        }
        else
        {
            rollDeg = kalmanRoll.update(rollAcc, gxDps, dt);
        }

        // Accelerometer pitch is only defined over +/-90, so past vertical the
        // gyro sign has to be flipped to stay consistent with it.
        if (fabsf(rollDeg) > 90.0f) gyDps = -gyDps;
        pitchDeg = kalmanPitch.update(pitchAcc, gyDps, dt);

        // Yaw: no measurement exists on a 6-axis IMU. Integrate and wrap.
        yawDeg += gzDps * dt;
        while (yawDeg > 180.0f) yawDeg -= 360.0f;
        while (yawDeg < -180.0f) yawDeg += 360.0f;
    }

    float getRoll() const { return rollDeg; }
    float getPitch() const { return pitchDeg; }
    float getYaw() const { return yawDeg; }

    float getRollRadians() const { return rollDeg * 0.017453293f; }
    float getPitchRadians() const { return pitchDeg * 0.017453293f; }
    float getYawRadians() const { return yawDeg * 0.017453293f; }

    // Bias estimated by each axis filter, in deg/s. Yaw has none.
    float getBiasX() const { return kalmanRoll.getBias(); }
    float getBiasY() const { return kalmanPitch.getBias(); }

   private:
    float rollDeg = 0.0f, pitchDeg = 0.0f, yawDeg = 0.0f;
    bool seeded = false;

    static float accelRoll(float ax, float ay, float az)
    {
        (void)ax;
        return atan2f(ay, az) * 57.29577951f;
    }

    static float accelPitch(float ax, float ay, float az)
    {
        return atan2f(-ax, sqrtf(ay * ay + az * az)) * 57.29577951f;
    }
};

#endif  // ORIENT_KALMAN_ANGLE_H
