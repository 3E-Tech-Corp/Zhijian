#pragma once

#include "common/types.hpp"
#include "common/cuda_utils.hpp"

namespace zhijian {

// Time integration schemes

// Abstract base class for time integrators
class TimeIntegrator {
public:
    virtual ~TimeIntegrator() = default;

    // Advance solution from U to U + dt * dU/dt
    // n: number of elements * solution points * variables
    virtual void advance(Real* U, const Real* dUdt, int n, Real dt,
                         cudaStream_t stream = 0) = 0;

    // Number of stages
    virtual int numStages() const = 0;

    // Get stage coefficients (for multi-stage methods)
    virtual void getStageCoeffs(int stage, Real& alpha, Real& beta) const = 0;
};

// 3rd-order Strong Stability Preserving Runge-Kutta (SSP-RK3)
// Also known as TVD-RK3
//
// u^(1) = u^n + dt * L(u^n)
// u^(2) = 3/4 * u^n + 1/4 * u^(1) + 1/4 * dt * L(u^(1))
// u^(n+1) = 1/3 * u^n + 2/3 * u^(2) + 2/3 * dt * L(u^(2))
class SSPRK3 {
public:
    SSPRK3() = default;

    // Get number of stages
    static constexpr int numStages() { return 3; }

    // Stage 1: U1 = U0 + dt * R(U0)
    static void stage1_CPU(Real* U1, const Real* U0, const Real* R, int n, Real dt);

    // Stage 2: U2 = 0.75*U0 + 0.25*U1 + 0.25*dt*R(U1)
    static void stage2_CPU(Real* U2, const Real* U0, const Real* U1, const Real* R,
                           int n, Real dt);

    // Stage 3: U = (1/3)*U0 + (2/3)*U2 + (2/3)*dt*R(U2)
    static void stage3_CPU(Real* U, const Real* U0, const Real* U2, const Real* R,
                           int n, Real dt);

    // GPU versions
    static void stage1_GPU(Real* U1, const Real* U0, const Real* R, int n, Real dt,
                           cudaStream_t stream = 0);
    static void stage2_GPU(Real* U2, const Real* U0, const Real* U1, const Real* R,
                           int n, Real dt, cudaStream_t stream = 0);
    static void stage3_GPU(Real* U, const Real* U0, const Real* U2, const Real* R,
                           int n, Real dt, cudaStream_t stream = 0);
};

// Forward Euler (1st order) - mainly for debugging
class ForwardEuler {
public:
    static constexpr int numStages() { return 1; }

    static void advance_CPU(Real* U, const Real* R, int n, Real dt);
    static void advance_GPU(Real* U, const Real* R, int n, Real dt,
                            cudaStream_t stream = 0);
};

// Classical 4th-order Runge-Kutta (RK4) - for comparison
class RK4 {
public:
    static constexpr int numStages() { return 4; }

    // Requires 4 temporary arrays for k1, k2, k3, k4
    static void advance_CPU(Real* U, Real* k1, Real* k2, Real* k3, Real* k4,
                            const Real* R, int n, Real dt,
                            const std::function<void(Real*, Real*)>& compute_rhs);
};

}  // namespace zhijian
