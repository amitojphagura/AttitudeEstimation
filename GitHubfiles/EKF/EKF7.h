/*
  EKF7.h  -  7-state Extended Kalman Filter attitude estimator for a 6-axis IMU.

  Part of the ORIENT testbed (Orientation Reference Inertial Estimation Testbed).
  Companion to the Madgwick and Mahony baselines: same input units, same output
  units, same Euler convention, so the four filters can be compared directly.

  State (7):   x = [ q0 q1 q2 q3  bx by bz ]^T
                 q  = body -> world rotation quaternion (same convention as the
                      Arduino MadgwickAHRS / Mahony libraries)
                 b  = gyroscope bias in rad/s, estimated online

  Process:     q_dot = 0.5 * W(q) * (omega_meas - b)
               b_dot = 0                      (random walk driven by Q)

  Measurement: the normalised accelerometer vector, predicted as the direction
               of gravity expressed in the body frame:
                 h(q) = [ 2(q1 q3 - q0 q2),
                          2(q0 q1 + q2 q3),
                          q0^2 - q1^2 - q2^2 + q3^2 ]

  Units in:    gyro in DEGREES per second, accelerometer in g (any scale; it is
               normalised internally). This matches how the Madgwick and Mahony
               Arduino libraries are fed, so the sketches stay interchangeable.
  Units out:   getRoll() / getPitch() / getYaw() return DEGREES, matching
               MadgwickAHRS. getRollRadians() etc. return radians.

  Dependency:  BasicLinearAlgebra (Tom Stewart), version 5.x, from the Arduino
               Library Manager.
*/

#ifndef ORIENT_EKF7_H
#define ORIENT_EKF7_H

#include <BasicLinearAlgebra.h>
#include <math.h>

class EKF7
{
   public:
    // ---- Tuning knobs -----------------------------------------------------
    // Gyroscope white noise, rad/s per sqrt(Hz). Drives the attitude block of Q.
    // Larger  -> trusts the accelerometer more, converges faster, noisier.
    float sigmaGyro = 0.02f;

    // Gyro bias random walk, rad/s^2 per sqrt(Hz). Drives the bias block of Q.
    // Larger  -> bias tracks changes faster but wanders more.
    float sigmaBias = 0.006f;

    // Accelerometer measurement noise, in normalised-gravity units.
    // Larger  -> attitude leans on the gyro, drifts more, rejects motion better.
    float sigmaAccel = 0.35f;

    // Adaptive measurement gating. During a rotation the accelerometer sees
    // gravity plus linear acceleration, so |a| departs from 1 g. R is inflated
    // in proportion to that departure instead of hard-gating the update, which
    // avoids the discontinuity that makes hard motion-gating oscillate.
    // Set to 0.0f to disable and get a plain fixed-R EKF.
    float accelGateGain = 12.0f;

    // ---- Lifecycle --------------------------------------------------------
    void begin()
    {
        q0 = 1.0f; q1 = 0.0f; q2 = 0.0f; q3 = 0.0f;
        bx = 0.0f; by = 0.0f; bz = 0.0f;

        P.Fill(0.0f);
        for (int i = 0; i < 4; ++i) P(i, i) = 1e-3f;   // attitude uncertainty
        for (int i = 4; i < 7; ++i) P(i, i) = 1e-4f;   // bias uncertainty

        anglesComputed = false;
    }

    // Seed the attitude from a single accelerometer sample (roll/pitch only;
    // a 6-axis IMU has no yaw reference). Call once after calibration so the
    // filter does not have to converge from level.
    void setFromAccel(float ax, float ay, float az)
    {
        float n = sqrtf(ax * ax + ay * ay + az * az);
        if (n < 1e-6f) return;
        ax /= n; ay /= n; az /= n;

        float roll  = atan2f(ay, az);
        float pitch = atan2f(-ax, sqrtf(ay * ay + az * az));

        float cr = cosf(roll * 0.5f),  sr = sinf(roll * 0.5f);
        float cp = cosf(pitch * 0.5f), sp = sinf(pitch * 0.5f);

        q0 = cr * cp;
        q1 = sr * cp;
        q2 = cr * sp;
        q3 = -sr * sp;
        normaliseQuat();
        anglesComputed = false;
    }

    // ---- One filter step --------------------------------------------------
    // gx,gy,gz in deg/s ; ax,ay,az in g ; dt in seconds.
    void update(float gxDps, float gyDps, float gzDps,
                float ax, float ay, float az, float dt)
    {
        if (dt <= 0.0f || dt > 0.5f) return;

        const float wxm = gxDps * DEG2RAD;
        const float wym = gyDps * DEG2RAD;
        const float wzm = gzDps * DEG2RAD;

        predict(wxm, wym, wzm, dt);
        correct(ax, ay, az);

        anglesComputed = false;
    }

    // ---- Output -----------------------------------------------------------
    float getRoll()  { computeAngles(); return roll  * RAD2DEG; }
    float getPitch() { computeAngles(); return pitch * RAD2DEG; }
    float getYaw()   { computeAngles(); return yaw   * RAD2DEG; }

    float getRollRadians()  { computeAngles(); return roll;  }
    float getPitchRadians() { computeAngles(); return pitch; }
    float getYawRadians()   { computeAngles(); return yaw;   }

    // Estimated gyro bias in deg/s, the quantity a plain Madgwick/Mahony
    // filter cannot recover.
    float getBiasX() const { return bx * RAD2DEG; }
    float getBiasY() const { return by * RAD2DEG; }
    float getBiasZ() const { return bz * RAD2DEG; }

    float getQ0() const { return q0; }
    float getQ1() const { return q1; }
    float getQ2() const { return q2; }
    float getQ3() const { return q3; }

   private:
    static constexpr float DEG2RAD = 0.017453293f;
    static constexpr float RAD2DEG = 57.29577951f;

    float q0 = 1.0f, q1 = 0.0f, q2 = 0.0f, q3 = 0.0f;
    float bx = 0.0f, by = 0.0f, bz = 0.0f;

    BLA::Matrix<7, 7, float> P;

    float roll = 0.0f, pitch = 0.0f, yaw = 0.0f;
    bool anglesComputed = false;

    void normaliseQuat()
    {
        float n = sqrtf(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
        if (n < 1e-9f) { q0 = 1.0f; q1 = q2 = q3 = 0.0f; return; }
        n = 1.0f / n;
        q0 *= n; q1 *= n; q2 *= n; q3 *= n;
    }

    // W(q), the 4x3 matrix with q_dot = 0.5 * W(q) * omega.
    void quatRateJacobian(BLA::Matrix<4, 3, float> &W) const
    {
        W(0, 0) = -q1; W(0, 1) = -q2; W(0, 2) = -q3;
        W(1, 0) =  q0; W(1, 1) = -q3; W(1, 2) =  q2;
        W(2, 0) =  q3; W(2, 1) =  q0; W(2, 2) = -q1;
        W(3, 0) = -q2; W(3, 1) =  q1; W(3, 2) =  q0;
    }

    void predict(float wxm, float wym, float wzm, float dt)
    {
        // Bias-corrected angular rate
        const float wx = wxm - bx;
        const float wy = wym - by;
        const float wz = wzm - bz;

        BLA::Matrix<4, 3, float> W;
        quatRateJacobian(W);

        // --- State propagation: q += dt * 0.5 * W(q) * omega ---------------
        const float h = 0.5f * dt;
        float nq0 = q0 + h * (-q1 * wx - q2 * wy - q3 * wz);
        float nq1 = q1 + h * ( q0 * wx - q3 * wy + q2 * wz);
        float nq2 = q2 + h * ( q3 * wx + q0 * wy - q1 * wz);
        float nq3 = q3 + h * (-q2 * wx + q1 * wy + q0 * wz);
        q0 = nq0; q1 = nq1; q2 = nq2; q3 = nq3;
        normaliseQuat();
        // bias is propagated unchanged

        // --- Jacobian Fmat = d f / d x ----------------------------------------
        BLA::Matrix<7, 7, float> Fmat;
        Fmat.Fill(0.0f);
        for (int i = 0; i < 7; ++i) Fmat(i, i) = 1.0f;

        // d q_dot / d q  =  0.5 * Omega(omega)
        Fmat(0, 1) += -h * wx;  Fmat(0, 2) += -h * wy;  Fmat(0, 3) += -h * wz;
        Fmat(1, 0) +=  h * wx;  Fmat(1, 2) +=  h * wz;  Fmat(1, 3) += -h * wy;
        Fmat(2, 0) +=  h * wy;  Fmat(2, 1) += -h * wz;  Fmat(2, 3) +=  h * wx;
        Fmat(3, 0) +=  h * wz;  Fmat(3, 1) +=  h * wy;  Fmat(3, 2) += -h * wx;

        // d q_dot / d b  =  -0.5 * W(q)
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 3; ++j) Fmat(i, 4 + j) = -h * W(i, j);

        // --- Process noise Q ------------------------------------------------
        BLA::Matrix<7, 7, float> Q;
        Q.Fill(0.0f);

        const float qg = 0.25f * dt * dt * sigmaGyro * sigmaGyro;
        BLA::Matrix<4, 4, float> Qq = W * (~W) * qg;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) Q(i, j) = Qq(i, j);

        const float qb = dt * sigmaBias * sigmaBias;
        for (int i = 4; i < 7; ++i) Q(i, i) = qb;

        P = Fmat * P * (~Fmat) + Q;
        symmetrise();
    }

    void correct(float ax, float ay, float az)
    {
        float norm = sqrtf(ax * ax + ay * ay + az * az);
        if (norm < 1e-6f) return;          // free fall or dead sensor
        const float aMag = norm;           // in g, before normalisation
        const float inv = 1.0f / norm;
        ax *= inv; ay *= inv; az *= inv;

        // Predicted gravity direction in the body frame
        BLA::Matrix<3, 1, float> h;
        h(0) = 2.0f * (q1 * q3 - q0 * q2);
        h(1) = 2.0f * (q0 * q1 + q2 * q3);
        h(2) = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;

        // Measurement Jacobian H (3x7); the bias columns are zero
        BLA::Matrix<3, 7, float> H;
        H.Fill(0.0f);
        H(0, 0) = -2.0f * q2;  H(0, 1) =  2.0f * q3;  H(0, 2) = -2.0f * q0;  H(0, 3) =  2.0f * q1;
        H(1, 0) =  2.0f * q1;  H(1, 1) =  2.0f * q0;  H(1, 2) =  2.0f * q3;  H(1, 3) =  2.0f * q2;
        H(2, 0) =  2.0f * q0;  H(2, 1) = -2.0f * q1;  H(2, 2) = -2.0f * q2;  H(2, 3) =  2.0f * q3;

        // Adaptive R: trust the accelerometer less the further |a| sits from 1 g
        float rScale = 1.0f;
        if (accelGateGain > 0.0f)
        {
            float dev = fabsf(aMag - 1.0f);
            rScale = 1.0f + accelGateGain * dev * dev * 100.0f;
        }
        const float r = sigmaAccel * sigmaAccel * rScale;

        BLA::Matrix<3, 1, float> y;
        y(0) = ax - h(0);
        y(1) = ay - h(1);
        y(2) = az - h(2);

        BLA::Matrix<7, 3, float> PHt = P * (~H);
        BLA::Matrix<3, 3, float> S = H * PHt;
        for (int i = 0; i < 3; ++i) S(i, i) += r;

        BLA::Matrix<3, 3, float> Sinv;
        if (!BLA::Invert(S, Sinv)) return;    // singular; skip this correction

        BLA::Matrix<7, 3, float> K = PHt * Sinv;
        BLA::Matrix<7, 1, float> dx = K * y;

        q0 += dx(0); q1 += dx(1); q2 += dx(2); q3 += dx(3);
        bx += dx(4); by += dx(5); bz += dx(6);
        normaliseQuat();

        // Joseph form: P = (I - K H) P (I - K H)^T + K R K^T.
        // The short form P = (I - K H) P loses positive definiteness within a
        // few thousand single-precision updates and the filter then diverges.
        BLA::Matrix<7, 7, float> IKH = K * H;
        for (int i = 0; i < 7; ++i)
            for (int j = 0; j < 7; ++j) IKH(i, j) = (i == j ? 1.0f : 0.0f) - IKH(i, j);

        P = IKH * P * (~IKH) + (K * (~K)) * r;
        symmetrise();
    }

    void symmetrise()
    {
        for (int i = 0; i < 7; ++i)
        {
            if (P(i, i) < 1e-12f) P(i, i) = 1e-12f;
            for (int j = i + 1; j < 7; ++j)
            {
                float m = 0.5f * (P(i, j) + P(j, i));
                P(i, j) = m;
                P(j, i) = m;
            }
        }
    }

    // Same Euler extraction as the Arduino MadgwickAHRS / Mahony libraries, so
    // roll/pitch/yaw mean exactly the same thing across every ORIENT sketch.
    void computeAngles()
    {
        if (anglesComputed) return;
        roll  = atan2f(q0 * q1 + q2 * q3, 0.5f - q1 * q1 - q2 * q2);
        pitch = asinf(constrainf(-2.0f * (q1 * q3 - q0 * q2), -1.0f, 1.0f));
        yaw   = atan2f(q1 * q2 + q0 * q3, 0.5f - q2 * q2 - q3 * q3);
        anglesComputed = true;
    }

    static float constrainf(float v, float lo, float hi)
    {
        return v < lo ? lo : (v > hi ? hi : v);
    }
};

#endif  // ORIENT_EKF7_H
