#pragma once

#include "common/types.hpp"
#include "physics/euler.hpp"
#include <memory>

namespace zhijian {

// Abstract base class for boundary conditions
class BoundaryCondition {
public:
    virtual ~BoundaryCondition() = default;

    // Get BC type
    virtual BCType type() const = 0;

    // Compute ghost state for inviscid flux computation
    // U_int: interior state
    // nx, ny: outward unit normal
    // Returns ghost state for Riemann solver
    virtual State ghostState(const State& U_int, Real nx, Real ny) const = 0;

    // Compute ghost state gradient for viscous flux computation
    // Returns whether gradients should be reflected
    virtual void ghostGradient(const State& dUdx_int, const State& dUdy_int,
                                Real nx, Real ny,
                                State& dUdx_ghost, State& dUdy_ghost) const {
        // Default: copy interior gradients
        dUdx_ghost = dUdx_int;
        dUdy_ghost = dUdy_int;
    }

    // Check if this BC is viscous (requires gradient computation)
    virtual bool isViscous() const { return false; }

protected:
    Real gamma_ = 1.4;
};

// No-slip wall (adiabatic)
// Velocity = 0 at wall, adiabatic condition for temperature
class NoSlipWallBC : public BoundaryCondition {
public:
    NoSlipWallBC(Real gamma = 1.4) { gamma_ = gamma; }

    BCType type() const override { return BCType::Wall; }

    State ghostState(const State& U_int, Real nx, Real ny) const override {
        // Reflect velocity to enforce no-slip
        IdealGas gas(gamma_);
        Real rho = U_int.rho();
        Real u = U_int.rhou() / rho;
        Real v = U_int.rhov() / rho;
        Real p = gas.pressure(U_int);

        // Ghost state has opposite velocity
        return gas.primToConserv(rho, -u, -v, p);
    }

    void ghostGradient(const State& dUdx_int, const State& dUdy_int,
                        Real nx, Real ny,
                        State& dUdx_ghost, State& dUdy_ghost) const override {
        // For adiabatic wall: dT/dn = 0
        // Copy gradients (BR2 will handle properly)
        dUdx_ghost = dUdx_int;
        dUdy_ghost = dUdy_int;
    }

    bool isViscous() const override { return true; }
};

// Isothermal wall
class IsothermalWallBC : public BoundaryCondition {
public:
    IsothermalWallBC(Real T_wall, Real gamma = 1.4)
        : T_wall_(T_wall) { gamma_ = gamma; }

    BCType type() const override { return BCType::Wall; }

    State ghostState(const State& U_int, Real nx, Real ny) const override {
        IdealGas gas(gamma_);
        Real rho = U_int.rho();
        Real u = U_int.rhou() / rho;
        Real v = U_int.rhov() / rho;

        // Use wall temperature for pressure: p = rho * R * T (R = 1 nondim)
        Real p_ghost = rho * T_wall_;

        return gas.primToConserv(rho, -u, -v, p_ghost);
    }

    bool isViscous() const override { return true; }

private:
    Real T_wall_;
};

// Slip wall (inviscid wall / symmetry)
class SlipWallBC : public BoundaryCondition {
public:
    SlipWallBC(Real gamma = 1.4) { gamma_ = gamma; }

    BCType type() const override { return BCType::SlipWall; }

    State ghostState(const State& U_int, Real nx, Real ny) const override {
        // Reflect normal velocity component
        IdealGas gas(gamma_);
        Real rho = U_int.rho();
        Real u = U_int.rhou() / rho;
        Real v = U_int.rhov() / rho;
        Real p = gas.pressure(U_int);

        // Normal velocity
        Real vn = u * nx + v * ny;

        // Reflected velocity: v_ghost = v_int - 2 * (v . n) * n
        Real u_ghost = u - 2.0 * vn * nx;
        Real v_ghost = v - 2.0 * vn * ny;

        return gas.primToConserv(rho, u_ghost, v_ghost, p);
    }
};

// Symmetry plane
class SymmetryBC : public BoundaryCondition {
public:
    SymmetryBC(Real gamma = 1.4) { gamma_ = gamma; }

    BCType type() const override { return BCType::Symmetry; }

    State ghostState(const State& U_int, Real nx, Real ny) const override {
        // Same as slip wall - reflect normal velocity
        IdealGas gas(gamma_);
        Real rho = U_int.rho();
        Real u = U_int.rhou() / rho;
        Real v = U_int.rhov() / rho;
        Real p = gas.pressure(U_int);

        Real vn = u * nx + v * ny;
        Real u_ghost = u - 2.0 * vn * nx;
        Real v_ghost = v - 2.0 * vn * ny;

        return gas.primToConserv(rho, u_ghost, v_ghost, p);
    }

    void ghostGradient(const State& dUdx_int, const State& dUdy_int,
                        Real nx, Real ny,
                        State& dUdx_ghost, State& dUdy_ghost) const override {
        // Reflect gradient normal component
        dUdx_ghost = dUdx_int;
        dUdy_ghost = dUdy_int;

        // For symmetry: d/dn(tangential) = 0, tangential d/dn(normal) = 0
        // This is handled implicitly by the weak formulation
    }
};

// Far-field (characteristic) boundary condition
class FarFieldBC : public BoundaryCondition {
public:
    FarFieldBC(Real rho_inf, Real u_inf, Real v_inf, Real p_inf, Real gamma = 1.4)
        : rho_inf_(rho_inf), u_inf_(u_inf), v_inf_(v_inf), p_inf_(p_inf) {
        gamma_ = gamma;
        IdealGas gas(gamma_);
        U_inf_ = gas.primToConserv(rho_inf, u_inf, v_inf, p_inf);
        c_inf_ = gas.soundSpeed(U_inf_);
    }

    // Constructor from Mach number and angle of attack
    static std::shared_ptr<FarFieldBC> fromMach(Real Mach, Real AoA_deg,
                                                  Real rho_inf, Real p_inf,
                                                  Real gamma = 1.4) {
        Real AoA = AoA_deg * M_PI / 180.0;
        IdealGas gas(gamma);
        Real c = sqrt(gamma * p_inf / rho_inf);
        Real vel = Mach * c;
        Real u = vel * cos(AoA);
        Real v = vel * sin(AoA);
        return std::make_shared<FarFieldBC>(rho_inf, u, v, p_inf, gamma);
    }

    BCType type() const override { return BCType::FarField; }

    State ghostState(const State& U_int, Real nx, Real ny) const override {
        IdealGas gas(gamma_);

        // Interior state
        Real rho_i = U_int.rho();
        Real u_i = U_int.rhou() / rho_i;
        Real v_i = U_int.rhov() / rho_i;
        Real p_i = gas.pressure(U_int);
        Real c_i = gas.soundSpeed(U_int);

        // Normal velocity (positive = outflow)
        Real vn_i = u_i * nx + v_i * ny;
        Real vn_inf = u_inf_ * nx + v_inf_ * ny;

        // Riemann invariants
        Real R_plus = vn_i + 2.0 * c_i / (gamma_ - 1.0);   // From interior
        Real R_minus = vn_inf - 2.0 * c_inf_ / (gamma_ - 1.0);  // From far-field

        // Boundary values
        Real vn_b = 0.5 * (R_plus + R_minus);
        Real c_b = 0.25 * (gamma_ - 1.0) * (R_plus - R_minus);

        Real rho_b, u_b, v_b, p_b;

        if (vn_b >= 0) {
            // Outflow: use interior tangential velocity and entropy
            Real vt_i = u_i * (-ny) + v_i * nx;  // tangential velocity
            Real s_i = p_i / pow(rho_i, gamma_);  // entropy

            rho_b = pow(c_b * c_b / (gamma_ * s_i), 1.0 / (gamma_ - 1.0));
            p_b = rho_b * c_b * c_b / gamma_;
            u_b = vn_b * nx - vt_i * ny;
            v_b = vn_b * ny + vt_i * nx;
        } else {
            // Inflow: use far-field tangential velocity and entropy
            Real vt_inf = u_inf_ * (-ny) + v_inf_ * nx;
            Real s_inf = p_inf_ / pow(rho_inf_, gamma_);

            rho_b = pow(c_b * c_b / (gamma_ * s_inf), 1.0 / (gamma_ - 1.0));
            p_b = rho_b * c_b * c_b / gamma_;
            u_b = vn_b * nx - vt_inf * ny;
            v_b = vn_b * ny + vt_inf * nx;
        }

        return gas.primToConserv(rho_b, u_b, v_b, p_b);
    }

    const State& farFieldState() const { return U_inf_; }

private:
    Real rho_inf_, u_inf_, v_inf_, p_inf_;
    Real c_inf_;
    State U_inf_;
};

// Subsonic inflow (total pressure and temperature specified)
class SubsonicInflowBC : public BoundaryCondition {
public:
    SubsonicInflowBC(Real p_total, Real T_total, Real nx_flow, Real ny_flow,
                      Real gamma = 1.4)
        : p_total_(p_total), T_total_(T_total), nx_flow_(nx_flow), ny_flow_(ny_flow) {
        gamma_ = gamma;
        // Normalize flow direction
        Real mag = sqrt(nx_flow_ * nx_flow_ + ny_flow_ * ny_flow_);
        nx_flow_ /= mag;
        ny_flow_ /= mag;
    }

    BCType type() const override { return BCType::Inflow; }

    State ghostState(const State& U_int, Real nx, Real ny) const override {
        IdealGas gas(gamma_);

        // Use interior pressure to determine Mach number
        Real p_int = gas.pressure(U_int);

        // Isentropic relations: p/p_total = (1 + (gamma-1)/2 * M^2)^(-gamma/(gamma-1))
        // Solve for M given p_int
        Real ratio = p_int / p_total_;
        Real M2 = 2.0 / (gamma_ - 1.0) * (pow(ratio, -(gamma_ - 1.0) / gamma_) - 1.0);
        M2 = fmax(M2, 0.0);  // Ensure non-negative
        Real M = sqrt(M2);

        // Static temperature from total temperature
        Real T = T_total_ / (1.0 + 0.5 * (gamma_ - 1.0) * M2);

        // Density from ideal gas law: p = rho * R * T (R = 1 nondim)
        Real rho = p_int / T;

        // Speed of sound and velocity
        Real c = sqrt(gamma_ * p_int / rho);
        Real vel = M * c;

        // Velocity components in flow direction
        Real u = vel * nx_flow_;
        Real v = vel * ny_flow_;

        return gas.primToConserv(rho, u, v, p_int);
    }

private:
    Real p_total_, T_total_;
    Real nx_flow_, ny_flow_;
};

// Subsonic outflow (static pressure specified)
class SubsonicOutflowBC : public BoundaryCondition {
public:
    SubsonicOutflowBC(Real p_static, Real gamma = 1.4)
        : p_static_(p_static) { gamma_ = gamma; }

    BCType type() const override { return BCType::Outflow; }

    State ghostState(const State& U_int, Real nx, Real ny) const override {
        IdealGas gas(gamma_);

        // Extrapolate density and velocity from interior
        Real rho = U_int.rho();
        Real u = U_int.rhou() / rho;
        Real v = U_int.rhov() / rho;

        // Use specified static pressure
        return gas.primToConserv(rho, u, v, p_static_);
    }

private:
    Real p_static_;
};

// Periodic boundary condition (handled specially in connectivity)
class PeriodicBC : public BoundaryCondition {
public:
    PeriodicBC(Real gamma = 1.4) { gamma_ = gamma; }

    BCType type() const override { return BCType::Periodic; }

    // For periodic BC, ghost state comes from the paired element
    // This is handled in the solver connectivity, not here
    State ghostState(const State& U_int, Real nx, Real ny) const override {
        return U_int;  // Placeholder - actual implementation in solver
    }
};

// Factory function to create boundary conditions from BCType enum
std::shared_ptr<BoundaryCondition> createBC(BCType type, const SimParams& params);

// Factory function to create boundary conditions from config-file specification
std::shared_ptr<BoundaryCondition> createBCFromSpec(
        const SimParams::BCSpec& spec, const SimParams& params);

}  // namespace zhijian
