#pragma once

#include <cstdint>
#include <cstddef>
#include <vector>
#include <array>
#include <memory>
#include <string>
#include <map>
#include <cmath>

// CUDA host/device macro - only use CUDA keywords when compiling with nvcc
#ifdef __CUDACC__
#define ZHIJIAN_HD __host__ __device__
#else
#define ZHIJIAN_HD
#endif

namespace zhijian {

// Basic numeric types
using Real = double;
using Index = int64_t;
using LocalIndex = int32_t;

// Maximum polynomial order supported
constexpr int MAX_POLY_ORDER = 5;

// Number of conservative variables for 2D compressible NS
constexpr int N_VARS = 4;  // rho, rho*u, rho*v, rho*E

// Variable indices
enum class Var : int {
    Density = 0,
    MomentumX = 1,
    MomentumY = 2,
    Energy = 3
};

// Element types
enum class ElementType : int {
    Triangle = 0,
    Quadrilateral = 1
};

// Boundary condition types
enum class BCType : int {
    Interior = 0,
    Wall = 1,           // No-slip wall
    SlipWall = 2,       // Slip wall (inviscid)
    Symmetry = 3,       // Symmetry plane
    FarField = 4,       // Far-field (characteristic)
    Inflow = 5,         // Subsonic inflow
    Outflow = 6,        // Subsonic outflow
    Periodic = 7        // Periodic boundary
};

// Flux type for FR scheme
enum class FluxType : int {
    DG = 0,             // Discontinuous Galerkin
    SD = 1,             // Spectral Difference
    HU = 2,             // Huynh's g2 scheme
    GA = 3              // Gauss-Legendre
};

// Riemann solver types
enum class RiemannSolver : int {
    Rusanov = 0,        // Local Lax-Friedrichs
    Roe = 1,            // Roe's approximate solver
    HLLC = 2            // HLLC solver
};

// 2D vector
struct Vec2 {
    Real x, y;

    ZHIJIAN_HD Vec2() : x(0), y(0) {}
    ZHIJIAN_HD Vec2(Real x_, Real y_) : x(x_), y(y_) {}

    ZHIJIAN_HD Vec2 operator+(const Vec2& v) const { return Vec2(x + v.x, y + v.y); }
    ZHIJIAN_HD Vec2 operator-(const Vec2& v) const { return Vec2(x - v.x, y - v.y); }
    ZHIJIAN_HD Vec2 operator*(Real s) const { return Vec2(x * s, y * s); }
    ZHIJIAN_HD Vec2 operator/(Real s) const { return Vec2(x / s, y / s); }
    ZHIJIAN_HD Real dot(const Vec2& v) const { return x * v.x + y * v.y; }
    ZHIJIAN_HD Real norm() const { return sqrt(x * x + y * y); }
    ZHIJIAN_HD Real norm2() const { return x * x + y * y; }
};

// State vector for compressible flow
struct State {
    Real data[N_VARS];

    ZHIJIAN_HD State() {
        for (int i = 0; i < N_VARS; ++i) data[i] = 0;
    }

    ZHIJIAN_HD Real& operator[](int i) { return data[i]; }
    ZHIJIAN_HD const Real& operator[](int i) const { return data[i]; }

    ZHIJIAN_HD Real& rho() { return data[0]; }
    ZHIJIAN_HD Real& rhou() { return data[1]; }
    ZHIJIAN_HD Real& rhov() { return data[2]; }
    ZHIJIAN_HD Real& rhoE() { return data[3]; }

    ZHIJIAN_HD const Real& rho() const { return data[0]; }
    ZHIJIAN_HD const Real& rhou() const { return data[1]; }
    ZHIJIAN_HD const Real& rhov() const { return data[2]; }
    ZHIJIAN_HD const Real& rhoE() const { return data[3]; }

    ZHIJIAN_HD State operator+(const State& s) const {
        State r;
        for (int i = 0; i < N_VARS; ++i) r.data[i] = data[i] + s.data[i];
        return r;
    }

    ZHIJIAN_HD State operator*(Real a) const {
        State r;
        for (int i = 0; i < N_VARS; ++i) r.data[i] = data[i] * a;
        return r;
    }
};

// Flux tensor (2D flux for each variable)
struct Flux2D {
    State Fx;  // Flux in x-direction
    State Fy;  // Flux in y-direction
};

// Simulation parameters
struct SimParams {
    // Physical parameters
    Real gamma = 1.4;           // Ratio of specific heats
    Real Pr = 0.72;             // Prandtl number
    Real Re = 1000.0;           // Reynolds number
    Real Mach_inf = 0.5;        // Free-stream Mach number
    Real T_inf = 300.0;         // Free-stream temperature (K)
    Real p_inf = 101325.0;      // Free-stream pressure (Pa)
    Real rho_inf = 1.225;       // Free-stream density (kg/m^3)
    Real AoA = 0.0;             // Angle of attack (degrees)

    // Numerical parameters
    int poly_order = 3;         // Polynomial order (P)
    FluxType flux_type = FluxType::HU;  // FR flux type
    RiemannSolver riemann = RiemannSolver::Roe;  // Riemann solver

    // Time stepping
    Real CFL = 0.5;             // CFL number
    Real dt = 0.0;              // Time step (0 = auto)
    Real t_final = 1.0;         // Final time
    int max_iter = 100000;      // Maximum iterations

    // Viscous
    bool viscous = true;        // Enable viscous terms

    // Output
    int output_freq = 100;      // Output frequency
    int restart_freq = 1000;    // Restart file frequency
    std::string output_dir = "./output";
    std::string case_name = "simulation";

    // Boundary condition specifications from config file
    // Each entry: tag -> {bc_type_string, key-value parameters}
    struct BCSpec {
        std::string type_str;           // "wall", "slipwall", "farfield", etc.
        Real p_static = 0.0;            // for outflow
        Real p_total = 0.0;             // for inflow
        Real T_total = 0.0;             // for inflow
        Real T_wall = 0.0;              // for isothermal wall
        Real dir_x = 1.0, dir_y = 0.0;  // for inflow direction
    };
    std::map<int, BCSpec> bc_specs;     // tag -> BC specification
};

// MPI partition info
struct PartitionInfo {
    int rank = 0;               // MPI rank
    int num_procs = 1;          // Total MPI processes
    int num_local_elems = 0;    // Number of local elements
    int num_ghost_elems = 0;    // Number of ghost elements
    std::vector<int> neighbor_ranks;  // Neighboring MPI ranks
};

}  // namespace zhijian
