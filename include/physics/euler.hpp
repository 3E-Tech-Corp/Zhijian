#pragma once

#include "common/types.hpp"

namespace zhijian {

// Ideal gas equations of state and flux computation
class IdealGas {
public:
    // Constructor with specific heat ratio
    ZHIJIAN_HD
    explicit IdealGas(Real gamma = 1.4) : gamma_(gamma) {}

    // Compute pressure from conservative variables
    ZHIJIAN_HD
    Real pressure(const State& U) const {
        Real rho = U.rho();
        Real u = U.rhou() / rho;
        Real v = U.rhov() / rho;
        Real E = U.rhoE() / rho;
        Real ke = 0.5 * (u * u + v * v);
        return (gamma_ - 1.0) * rho * (E - ke);
    }

    // Compute temperature from conservative variables
    // T = p / (rho * R) where R = (gamma - 1) * cv
    ZHIJIAN_HD
    Real temperature(const State& U, Real R_gas = 287.0) const {
        return pressure(U) / (U.rho() * R_gas);
    }

    // Compute speed of sound
    ZHIJIAN_HD
    Real soundSpeed(const State& U) const {
        Real p = pressure(U);
        return sqrt(gamma_ * p / U.rho());
    }

    // Compute Mach number
    ZHIJIAN_HD
    Real machNumber(const State& U) const {
        Real rho = U.rho();
        Real u = U.rhou() / rho;
        Real v = U.rhov() / rho;
        Real vel = sqrt(u * u + v * v);
        return vel / soundSpeed(U);
    }

    // Compute velocity components
    ZHIJIAN_HD
    void velocity(const State& U, Real& u, Real& v) const {
        u = U.rhou() / U.rho();
        v = U.rhov() / U.rho();
    }

    // Compute specific enthalpy h = (E + p/rho)
    ZHIJIAN_HD
    Real enthalpy(const State& U) const {
        return (U.rhoE() + pressure(U)) / U.rho();
    }

    // Compute total enthalpy H = rho*h = rhoE + p
    ZHIJIAN_HD
    Real totalEnthalpy(const State& U) const {
        return U.rhoE() + pressure(U);
    }

    // Compute entropy s = p / rho^gamma (up to a constant)
    ZHIJIAN_HD
    Real entropy(const State& U) const {
        Real p = pressure(U);
        return p / pow(U.rho(), gamma_);
    }

    // Convert primitive variables (rho, u, v, p) to conservative
    ZHIJIAN_HD
    State primToConserv(Real rho, Real u, Real v, Real p) const {
        State U;
        U.rho() = rho;
        U.rhou() = rho * u;
        U.rhov() = rho * v;
        U.rhoE() = p / (gamma_ - 1.0) + 0.5 * rho * (u * u + v * v);
        return U;
    }

    // Convert conservative to primitive variables
    ZHIJIAN_HD
    void conservToPrim(const State& U, Real& rho, Real& u, Real& v, Real& p) const {
        rho = U.rho();
        u = U.rhou() / rho;
        v = U.rhov() / rho;
        p = pressure(U);
    }

    // Compute inviscid flux in x-direction
    ZHIJIAN_HD
    State fluxX(const State& U) const {
        Real rho = U.rho();
        Real u = U.rhou() / rho;
        Real v = U.rhov() / rho;
        Real p = pressure(U);

        State F;
        F[0] = U.rhou();                    // rho * u
        F[1] = U.rhou() * u + p;            // rho * u^2 + p
        F[2] = U.rhou() * v;                // rho * u * v
        F[3] = (U.rhoE() + p) * u;          // (rhoE + p) * u
        return F;
    }

    // Compute inviscid flux in y-direction
    ZHIJIAN_HD
    State fluxY(const State& U) const {
        Real rho = U.rho();
        Real u = U.rhou() / rho;
        Real v = U.rhov() / rho;
        Real p = pressure(U);

        State G;
        G[0] = U.rhov();                    // rho * v
        G[1] = U.rhov() * u;                // rho * v * u
        G[2] = U.rhov() * v + p;            // rho * v^2 + p
        G[3] = (U.rhoE() + p) * v;          // (rhoE + p) * v
        return G;
    }

    // Compute normal flux F * nx + G * ny
    ZHIJIAN_HD
    State fluxNormal(const State& U, Real nx, Real ny) const {
        Real rho = U.rho();
        Real u = U.rhou() / rho;
        Real v = U.rhov() / rho;
        Real p = pressure(U);
        Real H = totalEnthalpy(U);
        Real vn = u * nx + v * ny;

        State Fn;
        Fn[0] = rho * vn;
        Fn[1] = rho * u * vn + p * nx;
        Fn[2] = rho * v * vn + p * ny;
        Fn[3] = H * vn;
        return Fn;
    }

    // Compute flux Jacobian eigenvalues for stability analysis
    ZHIJIAN_HD
    void eigenvalues(const State& U, Real nx, Real ny,
                     Real& lambda1, Real& lambda2, Real& lambda3, Real& lambda4) const {
        Real u, v;
        velocity(U, u, v);
        Real vn = u * nx + v * ny;
        Real c = soundSpeed(U);

        lambda1 = vn - c;
        lambda2 = vn;
        lambda3 = vn;
        lambda4 = vn + c;
    }

    // Maximum wave speed for stability
    ZHIJIAN_HD
    Real maxWaveSpeed(const State& U) const {
        Real u, v;
        velocity(U, u, v);
        Real vel = sqrt(u * u + v * v);
        return vel + soundSpeed(U);
    }

    Real gamma() const { return gamma_; }

private:
    Real gamma_;
};

// Riemann solvers for inviscid interface flux
class RiemannSolvers {
public:
    // Rusanov (Local Lax-Friedrichs) flux
    ZHIJIAN_HD
    static State rusanov(const State& UL, const State& UR,
                         Real nx, Real ny, Real gamma) {
        IdealGas gas(gamma);

        State FL = gas.fluxNormal(UL, nx, ny);
        State FR = gas.fluxNormal(UR, nx, ny);

        Real sL = gas.maxWaveSpeed(UL);
        Real sR = gas.maxWaveSpeed(UR);
        Real smax = fmax(sL, sR);

        State F;
        for (int i = 0; i < N_VARS; ++i) {
            F[i] = 0.5 * (FL[i] + FR[i]) - 0.5 * smax * (UR[i] - UL[i]);
        }
        return F;
    }

    // Roe flux with entropy fix
    ZHIJIAN_HD
    static State roe(const State& UL, const State& UR,
                     Real nx, Real ny, Real gamma) {
        IdealGas gas(gamma);

        // Roe averages
        Real rhoL = UL.rho();
        Real rhoR = UR.rho();
        Real uL = UL.rhou() / rhoL;
        Real vL = UL.rhov() / rhoL;
        Real uR = UR.rhou() / rhoR;
        Real vR = UR.rhov() / rhoR;
        Real HL = gas.totalEnthalpy(UL) / rhoL;
        Real HR = gas.totalEnthalpy(UR) / rhoR;

        Real sqrtRhoL = sqrt(rhoL);
        Real sqrtRhoR = sqrt(rhoR);
        Real denom = sqrtRhoL + sqrtRhoR;

        Real rhoRoe = sqrtRhoL * sqrtRhoR;
        Real uRoe = (sqrtRhoL * uL + sqrtRhoR * uR) / denom;
        Real vRoe = (sqrtRhoL * vL + sqrtRhoR * vR) / denom;
        Real HRoe = (sqrtRhoL * HL + sqrtRhoR * HR) / denom;

        Real vnRoe = uRoe * nx + vRoe * ny;
        Real qsqRoe = uRoe * uRoe + vRoe * vRoe;
        Real cRoe = sqrt((gamma - 1.0) * (HRoe - 0.5 * qsqRoe));

        // Eigenvalues
        Real lambda1 = fabs(vnRoe - cRoe);
        Real lambda2 = fabs(vnRoe);
        Real lambda3 = fabs(vnRoe);
        Real lambda4 = fabs(vnRoe + cRoe);

        // Entropy fix (Harten)
        Real eps = 0.1 * cRoe;
        if (lambda1 < eps) lambda1 = (lambda1 * lambda1 + eps * eps) / (2.0 * eps);
        if (lambda4 < eps) lambda4 = (lambda4 * lambda4 + eps * eps) / (2.0 * eps);

        // Wave strengths
        Real drho = rhoR - rhoL;
        Real du = uR - uL;
        Real dv = vR - vL;
        Real dp = gas.pressure(UR) - gas.pressure(UL);
        Real dvn = du * nx + dv * ny;

        Real alpha1 = (dp - rhoRoe * cRoe * dvn) / (2.0 * cRoe * cRoe);
        Real alpha2 = drho - dp / (cRoe * cRoe);
        Real alpha3 = rhoRoe * (dv * nx - du * ny);  // tangential
        Real alpha4 = (dp + rhoRoe * cRoe * dvn) / (2.0 * cRoe * cRoe);

        // Flux
        State FL = gas.fluxNormal(UL, nx, ny);
        State FR = gas.fluxNormal(UR, nx, ny);

        // Right eigenvectors times wave strengths times eigenvalues
        State dF;
        Real tx = -ny, ty = nx;  // tangent direction

        // Wave 1 (u - c)
        dF[0] = lambda1 * alpha1 * 1.0;
        dF[1] = lambda1 * alpha1 * (uRoe - cRoe * nx);
        dF[2] = lambda1 * alpha1 * (vRoe - cRoe * ny);
        dF[3] = lambda1 * alpha1 * (HRoe - cRoe * vnRoe);

        // Wave 2 (u - entropy wave)
        dF[0] += lambda2 * alpha2 * 1.0;
        dF[1] += lambda2 * alpha2 * uRoe;
        dF[2] += lambda2 * alpha2 * vRoe;
        dF[3] += lambda2 * alpha2 * 0.5 * qsqRoe;

        // Wave 3 (u - shear wave)
        dF[0] += lambda3 * alpha3 * 0.0;
        dF[1] += lambda3 * alpha3 * tx;
        dF[2] += lambda3 * alpha3 * ty;
        dF[3] += lambda3 * alpha3 * (uRoe * tx + vRoe * ty);

        // Wave 4 (u + c)
        dF[0] += lambda4 * alpha4 * 1.0;
        dF[1] += lambda4 * alpha4 * (uRoe + cRoe * nx);
        dF[2] += lambda4 * alpha4 * (vRoe + cRoe * ny);
        dF[3] += lambda4 * alpha4 * (HRoe + cRoe * vnRoe);

        State F;
        for (int i = 0; i < N_VARS; ++i) {
            F[i] = 0.5 * (FL[i] + FR[i]) - 0.5 * dF[i];
        }
        return F;
    }

    // HLLC flux
    ZHIJIAN_HD
    static State hllc(const State& UL, const State& UR,
                      Real nx, Real ny, Real gamma) {
        IdealGas gas(gamma);

        Real rhoL = UL.rho();
        Real rhoR = UR.rho();
        Real uL = UL.rhou() / rhoL;
        Real vL = UL.rhov() / rhoL;
        Real uR = UR.rhou() / rhoR;
        Real vR = UR.rhov() / rhoR;
        Real pL = gas.pressure(UL);
        Real pR = gas.pressure(UR);
        Real cL = gas.soundSpeed(UL);
        Real cR = gas.soundSpeed(UR);

        Real vnL = uL * nx + vL * ny;
        Real vnR = uR * nx + vR * ny;

        // Wave speed estimates (pressure-based)
        Real pstar = 0.5 * (pL + pR) - 0.5 * (vnR - vnL) * 0.5 * (rhoL + rhoR) * 0.5 * (cL + cR);
        pstar = fmax(pstar, 0.0);

        Real qL = (pstar <= pL) ? 1.0 : sqrt(1.0 + (gamma + 1.0) / (2.0 * gamma) * (pstar / pL - 1.0));
        Real qR = (pstar <= pR) ? 1.0 : sqrt(1.0 + (gamma + 1.0) / (2.0 * gamma) * (pstar / pR - 1.0));

        Real SL = vnL - cL * qL;
        Real SR = vnR + cR * qR;
        Real SM = (pR - pL + rhoL * vnL * (SL - vnL) - rhoR * vnR * (SR - vnR))
                  / (rhoL * (SL - vnL) - rhoR * (SR - vnR));

        // Compute flux
        if (SL >= 0.0) {
            return gas.fluxNormal(UL, nx, ny);
        } else if (SR <= 0.0) {
            return gas.fluxNormal(UR, nx, ny);
        } else if (SM >= 0.0) {
            // Star left state
            Real factor = rhoL * (SL - vnL) / (SL - SM);
            State Ustar;
            Ustar[0] = factor;
            Ustar[1] = factor * (uL + (SM - vnL) * nx);
            Ustar[2] = factor * (vL + (SM - vnL) * ny);
            Ustar[3] = factor * (UL.rhoE() / rhoL + (SM - vnL) * (SM + pL / (rhoL * (SL - vnL))));

            State FL = gas.fluxNormal(UL, nx, ny);
            State F;
            for (int i = 0; i < N_VARS; ++i) {
                F[i] = FL[i] + SL * (Ustar[i] - UL[i]);
            }
            return F;
        } else {
            // Star right state
            Real factor = rhoR * (SR - vnR) / (SR - SM);
            State Ustar;
            Ustar[0] = factor;
            Ustar[1] = factor * (uR + (SM - vnR) * nx);
            Ustar[2] = factor * (vR + (SM - vnR) * ny);
            Ustar[3] = factor * (UR.rhoE() / rhoR + (SM - vnR) * (SM + pR / (rhoR * (SR - vnR))));

            State FR = gas.fluxNormal(UR, nx, ny);
            State F;
            for (int i = 0; i < N_VARS; ++i) {
                F[i] = FR[i] + SR * (Ustar[i] - UR[i]);
            }
            return F;
        }
    }
};

}  // namespace zhijian
