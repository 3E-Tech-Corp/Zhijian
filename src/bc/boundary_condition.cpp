#include "bc/boundary_condition.hpp"
#include <algorithm>
#include <iostream>

namespace zhijian {

std::shared_ptr<BoundaryCondition> createBC(BCType type, const SimParams& params) {
    switch (type) {
        case BCType::Wall:
            return std::make_shared<NoSlipWallBC>(params.gamma);
        case BCType::SlipWall:
            return std::make_shared<SlipWallBC>(params.gamma);
        case BCType::Symmetry:
            return std::make_shared<SymmetryBC>(params.gamma);
        case BCType::FarField:
            return FarFieldBC::fromMach(params.Mach_inf, params.AoA,
                                         params.rho_inf, params.p_inf, params.gamma);
        case BCType::Inflow:
            return std::make_shared<SubsonicInflowBC>(
                params.p_inf * 1.2, params.T_inf, 1.0, 0.0, params.gamma);
        case BCType::Outflow:
            return std::make_shared<SubsonicOutflowBC>(params.p_inf, params.gamma);
        case BCType::Periodic:
            return std::make_shared<PeriodicBC>(params.gamma);
        default:
            std::cerr << "Warning: Unknown BC type " << static_cast<int>(type)
                      << ", defaulting to SlipWall\n";
            return std::make_shared<SlipWallBC>(params.gamma);
    }
}

std::shared_ptr<BoundaryCondition> createBCFromSpec(
        const SimParams::BCSpec& spec, const SimParams& params) {
    std::string t = spec.type_str;

    if (t == "wall" || t == "noslip" || t == "no_slip" || t == "no-slip") {
        if (spec.T_wall > 0)
            return std::make_shared<IsothermalWallBC>(spec.T_wall, params.gamma);
        return std::make_shared<NoSlipWallBC>(params.gamma);
    }
    if (t == "slipwall" || t == "slip_wall" || t == "slip" || t == "inviscid_wall")
        return std::make_shared<SlipWallBC>(params.gamma);
    if (t == "symmetry" || t == "sym")
        return std::make_shared<SymmetryBC>(params.gamma);
    if (t == "farfield" || t == "far_field" || t == "far-field" || t == "freestream")
        return FarFieldBC::fromMach(params.Mach_inf, params.AoA,
                                     params.rho_inf, params.p_inf, params.gamma);
    if (t == "inflow" || t == "inlet") {
        Real pt = (spec.p_total > 0) ? spec.p_total : params.p_inf * 1.2;
        Real Tt = (spec.T_total > 0) ? spec.T_total : params.T_inf;
        return std::make_shared<SubsonicInflowBC>(
            pt, Tt, spec.dir_x, spec.dir_y, params.gamma);
    }
    if (t == "outflow" || t == "outlet" || t == "exit") {
        Real ps = (spec.p_static > 0) ? spec.p_static : params.p_inf;
        return std::make_shared<SubsonicOutflowBC>(ps, params.gamma);
    }
    if (t == "periodic")
        return std::make_shared<PeriodicBC>(params.gamma);

    std::cerr << "Warning: Unknown BC type string '" << t
              << "', defaulting to SlipWall\n";
    return std::make_shared<SlipWallBC>(params.gamma);
}

}  // namespace zhijian
