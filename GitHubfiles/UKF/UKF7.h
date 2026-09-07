/*
  UKF7.h  -  7-state Unscented Kalman Filter attitude estimator for a 6-axis IMU.

  Part of the ORIENT testbed (Orientation Reference Inertial Estimation Testbed).
  Same state, same measurement model, same units and same Euler convention as
  EKF7.h, so the only difference between the two sketches is how the nonlinearity
  is propagated: analytic Jacobians (EKF) versus sigma points (UKF).

  State (7):   x = [ q0 q1 q2 q3  bx by bz ]^T
  Sigma points: 2n+1 = 15, scaled unscented transform.

  Scaling parameters
  ------------------
  The textbook default alpha = 1e-3 makes (n + lambda) = alpha^2 * (n + kappa)
  about 7e-6 for n = 7, so W0 becomes roughly -1e6. That is fine in double
  precision and useless in the ESP32's single-precision FPU. This filter
  therefore defaults to alpha = 1, kappa = 0, beta = 2, which gives lambda = 0,
  W0_mean = 0, Wi = 1/14 and W0_cov = 2. All weights stay well conditioned in
  float32. Change them only if you know why.

  Units in:    gyro in DEGREES per second, accelerometer in g.
  Units out:   getRoll() / getPitch() / getYaw() in DEGREES, matching
               MadgwickAHRS, Mahony and EKF7.

  Dependency:  BasicLinearAlgebra (Tom Stewart), version 5.x, from the Arduino
               Library Manager.
*/

#ifndef ORIENT_UKF7_H
#define ORIENT_UKF7_H

#include <BasicLinearAlgebra.h>
#include <math.h>

class UKF7
{
   public:
    // ---- Tuning knobs (same meaning and same defaults as EKF7) ------------
    float sigmaGyro = 0.02f;      // rad/s /sqrt(Hz)
    float sigmaBias = 0.006f;     // rad/s^2 /sqrt(Hz)
    float sigmaAccel = 0.35f;     // normalised gravity units
    float accelGateGain = 12.0f;  // 0 disables the adaptive R inflation

    // Unscented transform scaling. See the note above before changing these.
    float alpha = 1.0f;
    float beta = 2.0f;
    float kappa = 0.0f;

    // ---- Lifecycle --------------------------------------------------------
    void begin()
    {
        x.Fill(0.0f);
        x(0) = 1.0f;

        P.Fill(0.0f);
        for (int i = 0; i < 4; ++i) P(i, i) = 1e-3f;
        for (int i = 4; i < 7; ++i) P(i, i) = 1e-4f;

        computeWeights();
        anglesComputed = false;
    }

    void setFromAccel(float ax, float ay, float az)
    {
        float n = sqrtf(ax * ax + ay * ay + az * az);
        if (n < 1e-6f) return;
        ax /= n; ay /= n; az /= n;

        float roll = atan2f(ay, az);
        float pitch = atan2f(-ax, sqrtf(ay * ay + az * az));

        float cr = cosf(roll * 0.5f), sr = sinf(roll * 0.5f);
        float cp = cosf(pitch * 0.5f), sp = sinf(pitch * 0.5f);

        x(0) = cr * cp;
        x(1) = sr * cp;
        x(2) = cr * sp;
        x(3) = -sr * sp;
        normaliseQuat(x);
        anglesComputed = false;
    }

    // ---- One filter step --------------------------------------------------
    void update(float gxDps, float gyDps, float gzDps,
                float ax, float ay, float az, float dt)
    {
        if (dt <= 0.0f || dt > 0.5f) return;

        const float wxm = gxDps * DEG2RAD;
        const float wym = gyDps * DEG2RAD;
        const float wzm = gzDps * DEG2RAD;

        if (!generateSigmaPoints()) return;
        propagate(wxm, wym, wzm, dt);
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

    float getBiasX() const { return x(4) * RAD2DEG; }
    float getBiasY() const { return x(5) * RAD2DEG; }
    float getBiasZ() const { return x(6) * RAD2DEG; }

    float getQ0() const { return x(0); }
    float getQ1() const { return x(1); }
    float getQ2() const { return x(2); }
    float getQ3() const { return x(3); }

    // True when the last Cholesky factorisation failed and the covariance had
    // to be reconditioned. Useful to log while tuning.
    bool lastFactorisationRepaired() const { return repaired; }

   private:
    static const int N = 7;
    static const int NSIGMA = 15;  // 2N + 1

    static constexpr float DEG2RAD = 0.017453293f;
    static constexpr float RAD2DEG = 57.29577951f;

    BLA::Matrix<7, 1, float> x;
    BLA::Matrix<7, 7, float> P;
    BLA::Matrix<7, 7, float> S;    // lower-triangular sqrt of (N + lambda) P
    BLA::Matrix<7, 15, float> X;   // sigma points
    BLA::Matrix<3, 15, float> Z;   // measurement sigma points

    float Wm[NSIGMA], Wc[NSIGMA];
    float lambda = 0.0f;
    bool repaired = false;

    float roll = 0.0f, pitch = 0.0f, yaw = 0.0f;
    bool anglesComputed = false;

    void computeWeights()
    {
        lambda = alpha * alpha * (N + kappa) - N;
        const float denom = N + lambda;
        Wm[0] = lambda / denom;
        Wc[0] = Wm[0] + (1.0f - alpha * alpha + beta);
        for (int i = 1; i < NSIGMA; ++i)
        {
            Wm[i] = 0.5f / denom;
            Wc[i] = Wm[i];
        }
    }

    static void normaliseQuat(BLA::Matrix<7, 1, float> &v)
    {
        float n = sqrtf(v(0) * v(0) + v(1) * v(1) + v(2) * v(2) + v(3) * v(3));
        if (n < 1e-9f) { v(0) = 1.0f; v(1) = v(2) = v(3) = 0.0f; return; }
        n = 1.0f / n;
        for (int i = 0; i < 4; ++i) v(i) *= n;
    }

    // ---- Step 1: sigma points from chol((N + lambda) * P) -----------------
    bool generateSigmaPoints()
    {
        computeWeights();  // pick up any change to alpha/beta/kappa at runtime
        repaired = false;

        const float scale = N + lambda;
        BLA::Matrix<7, 7, float> A = P * scale;

        if (!choleskyInto(A, S))
        {
            // Recondition: bump the diagonal and try once more. Single-precision
            // round-off can push P slightly indefinite after many updates.
            repaired = true;
            for (int i = 0; i < N; ++i) A(i, i) += 1e-6f * scale;

            if (!choleskyInto(A, S))
            {
                // Still not usable: reset the covariance to its startup value
                // rather than let the filter diverge.
                P.Fill(0.0f);
                for (int i = 0; i < 4; ++i) P(i, i) = 1e-3f;
                for (int i = 4; i < 7; ++i) P(i, i) = 1e-4f;
                A = P * scale;
                if (!choleskyInto(A, S)) return false;
            }
        }

        for (int r = 0; r < N; ++r)
        {
            X(r, 0) = x(r);
            for (int c = 0; c < N; ++c)
            {
                float s = S(r, c);   // lower-triangular factor, 0 above diagonal
                X(r, 1 + c) = x(r) + s;
                X(r, 1 + N + c) = x(r) - s;
            }
        }
        return true;
    }

    // CholeskyDecompose works in place and its result object holds a reference
    // to the matrix, so take a copy and hand back a plain lower-triangular
    // matrix the caller owns.
    static bool choleskyInto(BLA::Matrix<7, 7, float> A, BLA::Matrix<7, 7, float> &L)
    {
        auto chol = BLA::CholeskyDecompose(A);
        if (!chol.positive_definite) return false;
        for (int r = 0; r < N; ++r)
            for (int c = 0; c < N; ++c) L(r, c) = chol.L(r, c);
        return true;
    }

    // ---- Step 2: push every sigma point through the process model ----------
    void propagate(float wxm, float wym, float wzm, float dt)
    {
        const float h = 0.5f * dt;

        for (int c = 0; c < NSIGMA; ++c)
        {
            float q0 = X(0, c), q1 = X(1, c), q2 = X(2, c), q3 = X(3, c);
            const float wx = wxm - X(4, c);
            const float wy = wym - X(5, c);
            const float wz = wzm - X(6, c);

            float n0 = q0 + h * (-q1 * wx - q2 * wy - q3 * wz);
            float n1 = q1 + h * ( q0 * wx - q3 * wy + q2 * wz);
            float n2 = q2 + h * ( q3 * wx + q0 * wy - q1 * wz);
            float n3 = q3 + h * (-q2 * wx + q1 * wy + q0 * wz);

            float nn = sqrtf(n0 * n0 + n1 * n1 + n2 * n2 + n3 * n3);
            if (nn < 1e-9f) { n0 = 1.0f; n1 = n2 = n3 = 0.0f; nn = 1.0f; }
            nn = 1.0f / nn;

            X(0, c) = n0 * nn;
            X(1, c) = n1 * nn;
            X(2, c) = n2 * nn;
            X(3, c) = n3 * nn;
            // bias states propagate unchanged
        }

        // q and -q are the same rotation. Flip any sigma quaternion sitting in
        // the opposite hemisphere from the reference point, otherwise the
        // weighted average below is meaningless.
        for (int c = 1; c < NSIGMA; ++c)
        {
            float d = X(0, c) * X(0, 0) + X(1, c) * X(1, 0) +
                      X(2, c) * X(2, 0) + X(3, c) * X(3, 0);
            if (d < 0.0f)
                for (int r = 0; r < 4; ++r) X(r, c) = -X(r, c);
        }

        // Predicted mean
        for (int r = 0; r < N; ++r)
        {
            float s = 0.0f;
            for (int c = 0; c < NSIGMA; ++c) s += Wm[c] * X(r, c);
            x(r) = s;
        }
        normaliseQuat(x);

        // Predicted covariance plus process noise
        P.Fill(0.0f);
        for (int c = 0; c < NSIGMA; ++c)
        {
            float d[N];
            for (int r = 0; r < N; ++r) d[r] = X(r, c) - x(r);
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j) P(i, j) += Wc[c] * d[i] * d[j];
        }
        addProcessNoise(dt);
        symmetrise();
    }

    // Q, built the same way as in EKF7 so the two filters are given identical
    // noise assumptions.
    void addProcessNoise(float dt)
    {
        const float q0 = x(0), q1 = x(1), q2 = x(2), q3 = x(3);
        float W[4][3] = {
            {-q1, -q2, -q3},
            { q0, -q3,  q2},
            { q3,  q0, -q1},
            {-q2,  q1,  q0}};

        const float qg = 0.25f * dt * dt * sigmaGyro * sigmaGyro;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
            {
                float s = 0.0f;
                for (int k = 0; k < 3; ++k) s += W[i][k] * W[j][k];
                P(i, j) += qg * s;
            }

        const float qb = dt * sigmaBias * sigmaBias;
        for (int i = 4; i < 7; ++i) P(i, i) += qb;
    }

    // ---- Step 3: unscented measurement update ------------------------------
    void correct(float ax, float ay, float az)
    {
        float norm = sqrtf(ax * ax + ay * ay + az * az);
        if (norm < 1e-6f) return;
        const float aMag = norm;
        const float inv = 1.0f / norm;
        ax *= inv; ay *= inv; az *= inv;

        // Regenerate sigma points around the predicted state so the measurement
        // spread reflects the predicted covariance.
        if (!generateSigmaPoints()) return;

        for (int c = 0; c < NSIGMA; ++c)
        {
            float q0 = X(0, c), q1 = X(1, c), q2 = X(2, c), q3 = X(3, c);
            float n = sqrtf(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
            if (n > 1e-9f) { n = 1.0f / n; q0 *= n; q1 *= n; q2 *= n; q3 *= n; }
            Z(0, c) = 2.0f * (q1 * q3 - q0 * q2);
            Z(1, c) = 2.0f * (q0 * q1 + q2 * q3);
            Z(2, c) = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;
        }

        float zbar[3];
        for (int r = 0; r < 3; ++r)
        {
            float s = 0.0f;
            for (int c = 0; c < NSIGMA; ++c) s += Wm[c] * Z(r, c);
            zbar[r] = s;
        }

        float rScale = 1.0f;
        if (accelGateGain > 0.0f)
        {
            float dev = fabsf(aMag - 1.0f);
            rScale = 1.0f + accelGateGain * dev * dev * 100.0f;
        }
        const float r = sigmaAccel * sigmaAccel * rScale;

        BLA::Matrix<3, 3, float> Pzz;
        Pzz.Fill(0.0f);
        BLA::Matrix<7, 3, float> Pxz;
        Pxz.Fill(0.0f);

        for (int c = 0; c < NSIGMA; ++c)
        {
            float dz[3];
            for (int i = 0; i < 3; ++i) dz[i] = Z(i, c) - zbar[i];
            float dxv[N];
            for (int i = 0; i < N; ++i) dxv[i] = X(i, c) - x(i);

            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) Pzz(i, j) += Wc[c] * dz[i] * dz[j];
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < 3; ++j) Pxz(i, j) += Wc[c] * dxv[i] * dz[j];
        }
        for (int i = 0; i < 3; ++i) Pzz(i, i) += r;

        BLA::Matrix<3, 3, float> Pzzinv;
        if (!BLA::Invert(Pzz, Pzzinv)) return;

        BLA::Matrix<7, 3, float> K = Pxz * Pzzinv;

        BLA::Matrix<3, 1, float> y;
        y(0) = ax - zbar[0];
        y(1) = ay - zbar[1];
        y(2) = az - zbar[2];

        x += K * y;
        normaliseQuat(x);

        P -= K * Pzz * (~K);
        symmetrise();
    }

    void symmetrise()
    {
        for (int i = 0; i < N; ++i)
        {
            if (P(i, i) < 1e-12f) P(i, i) = 1e-12f;
            for (int j = i + 1; j < N; ++j)
            {
                float m = 0.5f * (P(i, j) + P(j, i));
                P(i, j) = m;
                P(j, i) = m;
            }
        }
    }

    void computeAngles()
    {
        if (anglesComputed) return;
        const float q0 = x(0), q1 = x(1), q2 = x(2), q3 = x(3);
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

#endif  // ORIENT_UKF7_H
