#pragma once

#include "common/types.hpp"
#include "physics/euler.hpp"

namespace zhijian {

// Viscous flux computation for Navier-Stokes equations
class ViscousFlux {
public:
    // Constructor with physical parameters
    __host__ __device__
    ViscousFlux(Real gamma, Real Re, Real Pr, Real mu_ref = 1.0, Real T_ref = 1.0)
        : gamma_(gamma), Re_(Re), Pr_(Pr), mu_ref_(mu_ref), T_ref_(T_ref) {}

    // Compute dynamic viscosity using Sutherland's law
    __host__ __device__
    Real dynamicViscosity(Real T) const {
        // Sutherland's law: mu = mu_ref * (T/T_ref)^1.5 * (T_ref + S) / (T + S)
        // where S = 110.4 K for air
        const Real S = 110.4 / T_ref_;  // Normalized Sutherland constant
        Real T_norm = T / T_ref_;
        return mu_ref_ * pow(T_norm, 1.5) * (1.0 + S) / (T_norm + S);
    }

    // Compute thermal conductivity
    // k = mu * cp / Pr, where cp = gamma * R / (gamma - 1)
    __host__ __device__
    Real thermalConductivity(Real mu) const {
        Real cp = gamma_ / (gamma_ - 1.0);  // Non-dimensional cp
        return mu * cp / Pr_;
    }

    // Compute viscous stress tensor components
    // tau_xx, tau_xy, tau_yy
    __host__ __device__
    void stressTensor(Real mu, Real dudx, Real dudy, Real dvdx, Real dvdy,
                      Real& tau_xx, Real& tau_xy, Real& tau_yy) const {
        // Using Stokes hypothesis (bulk viscosity = 0)
        Real div_v = dudx + dvdy;
        Real lambda = -2.0 / 3.0 * mu;  // Second viscosity coefficient

        tau_xx = 2.0 * mu * dudx + lambda * div_v;
        tau_yy = 2.0 * mu * dvdy + lambda * div_v;
        tau_xy = mu * (dudy + dvdx);
    }

    // Compute viscous flux in x-direction
    // grad_U = [drho/dx, drhou/dx, drhov/dx, drhoE/dx] (not used directly)
    // We need velocity and temperature gradients
    __host__ __device__
    State fluxX(const State& U, Real dudx, Real dudy, Real dvdx, Real dvdy,
                Real dTdx, Real dTdy) const {
        IdealGas gas(gamma_);
        Real rho = U.rho();
        Real u = U.rhou() / rho;
        Real v = U.rhov() / rho;
        Real T = gas.pressure(U) / rho;  // Using p = rho * R * T with R = 1

        Real mu = dynamicViscosity(T) / Re_;
        Real k = thermalConductivity(mu) / Re_;

        Real tau_xx, tau_xy, tau_yy;
        stressTensor(mu, dudx, dudy, dvdx, dvdy, tau_xx, tau_xy, tau_yy);

        State Fv;
        Fv[0] = 0.0;                                    // No mass diffusion
        Fv[1] = tau_xx;                                 // x-momentum viscous flux
        Fv[2] = tau_xy;                                 // y-momentum viscous flux
        Fv[3] = u * tau_xx + v * tau_xy + k * dTdx;    // Energy viscous flux
        return Fv;
    }

    // Compute viscous flux in y-direction
    __host__ __device__
    State fluxY(const State& U, Real dudx, Real dudy, Real dvdx, Real dvdy,
                Real dTdx, Real dTdy) const {
        IdealGas gas(gamma_);
        Real rho = U.rho();
        Real u = U.rhou() / rho;
        Real v = U.rhov() / rho;
        Real T = gas.pressure(U) / rho;

        Real mu = dynamicViscosity(T) / Re_;
        Real k = thermalConductivity(mu) / Re_;

        Real tau_xx, tau_xy, tau_yy;
        stressTensor(mu, dudx, dudy, dvdx, dvdy, tau_xx, tau_xy, tau_yy);

        State Gv;
        Gv[0] = 0.0;
        Gv[1] = tau_xy;
        Gv[2] = tau_yy;
        Gv[3] = u * tau_xy + v * tau_yy + k * dTdy;
        return Gv;
    }

    // Compute viscous normal flux at interface
    // Using BR2 (Bassi-Rebay 2) or interior penalty method
    __host__ __device__
    State fluxNormal(const State& U, Real nx, Real ny,
                     Real dudx, Real dudy, Real dvdx, Real dvdy,
                     Real dTdx, Real dTdy) const {
        State Fv = fluxX(U, dudx, dudy, dvdx, dvdy, dTdx, dTdy);
        State Gv = fluxY(U, dudx, dudy, dvdx, dvdy, dTdx, dTdy);

        State Fn;
        for (int i = 0; i < N_VARS; ++i) {
            Fn[i] = Fv[i] * nx + Gv[i] * ny;
        }
        return Fn;
    }

    // Compute velocity gradients from conservative variable gradients
    __host__ __device__
    static void velocityGradients(const State& U,
                                  const State& dUdx, const State& dUdy,
                                  Real& dudx, Real& dudy,
                                  Real& dvdx, Real& dvdy) {
        Real rho = U.rho();
        Real u = U.rhou() / rho;
        Real v = U.rhov() / rho;

        // d(rho*u)/dx = rho*du/dx + u*drho/dx
        // => du/dx = (d(rho*u)/dx - u*drho/dx) / rho
        dudx = (dUdx[1] - u * dUdx[0]) / rho;
        dudy = (dUdy[1] - u * dUdy[0]) / rho;
        dvdx = (dUdx[2] - v * dUdx[0]) / rho;
        dvdy = (dUdy[2] - v * dUdy[0]) / rho;
    }

    // Compute temperature gradient from conservative variable gradients
    __host__ __device__
    void temperatureGradient(const State& U,
                             const State& dUdx, const State& dUdy,
                             Real& dTdx, Real& dTdy) const {
        Real rho = U.rho();
        Real u = U.rhou() / rho;
        Real v = U.rhov() / rho;

        IdealGas gas(gamma_);
        Real p = gas.pressure(U);

        // Temperature T = p / rho (with R = 1)
        // dT/dx = (dp/dx - T*drho/dx) / rho

        // dp/dx from p = (gamma-1)*(rhoE - 0.5*rho*(u^2+v^2))
        Real dudx, dudy, dvdx, dvdy;
        velocityGradients(U, dUdx, dUdy, dudx, dudy, dvdx, dvdy);

        Real dpdx = (gamma_ - 1.0) * (dUdx[3] - 0.5 * dUdx[0] * (u*u + v*v)
                                       - rho * (u * dudx + v * dvdx));
        Real dpdy = (gamma_ - 1.0) * (dUdy[3] - 0.5 * dUdy[0] * (u*u + v*v)
                                       - rho * (u * dudy + v * dvdy));

        Real T = p / rho;
        dTdx = (dpdx - T * dUdx[0]) / rho;
        dTdy = (dpdy - T * dUdy[0]) / rho;
    }

private:
    Real gamma_;
    Real Re_;
    Real Pr_;
    Real mu_ref_;
    Real T_ref_;
};

// Interface viscous flux computation using BR2 scheme
class BR2Scheme {
public:
    // Compute common viscous flux at interface
    // UL, UR: states from left and right
    // gradUL, gradUR: gradients from left and right
    // eta: penalty parameter
    __host__ __device__
    static State commonFlux(const State& UL, const State& UR,
                            const State& dUdxL, const State& dUdyL,
                            const State& dUdxR, const State& dUdyR,
                            Real nx, Real ny,
                            Real gamma, Real Re, Real Pr,
                            Real eta = 1.0) {
        // Average state and gradients
        State U_avg;
        State dUdx_avg, dUdy_avg;
        for (int i = 0; i < N_VARS; ++i) {
            U_avg[i] = 0.5 * (UL[i] + UR[i]);
            dUdx_avg[i] = 0.5 * (dUdxL[i] + dUdxR[i]);
            dUdy_avg[i] = 0.5 * (dUdyL[i] + dUdyR[i]);
        }

        ViscousFlux visc(gamma, Re, Pr);

        // Compute velocity and temperature gradients
        Real dudx, dudy, dvdx, dvdy;
        ViscousFlux::velocityGradients(U_avg, dUdx_avg, dUdy_avg,
                                        dudx, dudy, dvdx, dvdy);

        Real dTdx, dTdy;
        visc.temperatureGradient(U_avg, dUdx_avg, dUdy_avg, dTdx, dTdy);

        // Compute viscous flux at average state
        return visc.fluxNormal(U_avg, nx, ny, dudx, dudy, dvdx, dvdy, dTdx, dTdy);
    }
};

}  // namespace zhijian
