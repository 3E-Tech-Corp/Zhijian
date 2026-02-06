#pragma once

#include "common/types.hpp"
#include "common/cuda_utils.hpp"
#include "mesh/mesh.hpp"
#include "basis/basis.hpp"
#include "physics/euler.hpp"
#include "physics/navierstokes.hpp"
#include <memory>
#include <vector>
#include <functional>
#include <map>

namespace zhijian {

// Forward declarations
class BoundaryCondition;
class MPICommunicator;

// Solution data for an element
struct ElementSolution {
    std::vector<State> sol_pts;     // Solution at solution points
    std::vector<State> flux_pts;    // Solution at flux points (per edge)
    std::vector<State> dUdx;        // x-gradient at solution points
    std::vector<State> dUdy;        // y-gradient at solution points
};

// GPU data structure for the solver
struct GPUSolverData {
    // Solution arrays
    DeviceArray<Real> U;            // Conservative variables [n_elem * n_sp * N_VARS]
    DeviceArray<Real> dUdt;         // Time derivative [n_elem * n_sp * N_VARS]
    DeviceArray<Real> dUdx;         // x-gradient [n_elem * n_sp * N_VARS]
    DeviceArray<Real> dUdy;         // y-gradient [n_elem * n_sp * N_VARS]

    // Inviscid flux at solution points (for flux divergence)
    DeviceArray<Real> Fx_sp;        // x-flux at solution points
    DeviceArray<Real> Fy_sp;        // y-flux at solution points

    // Flux point data
    DeviceArray<Real> U_fp;         // Solution at flux points
    DeviceArray<Real> F_fp;         // Numerical flux at flux points
    DeviceArray<Real> U_ghost;      // Ghost states at boundary faces

    // RK stages
    DeviceArray<Real> U0;           // Initial solution for RK
    DeviceArray<Real> U1;           // Intermediate stage

    // Geometry
    DeviceArray<Real> J;            // Jacobian determinant at solution points
    DeviceArray<Real> Jinv;         // Inverse Jacobian components
    DeviceArray<Real> face_normals; // Face normals
    DeviceArray<Real> face_lengths; // Face lengths

    // Connectivity
    DeviceArray<int> elem_type;     // Element types
    DeviceArray<int> face_info;     // Face connectivity info
    DeviceArray<int> bc_type;       // Boundary condition types

    // Face connectivity for Riemann flux
    DeviceArray<int> face_left_elem;   // Left element index for each face
    DeviceArray<int> face_left_local;  // Left local face index (0-3)
    DeviceArray<int> face_right_elem;  // Right element index (-1 for boundary)
    DeviceArray<int> face_right_local; // Right local face index

    // FR operators (precomputed)
    DeviceArray<Real> diff_xi;      // Differentiation matrix
    DeviceArray<Real> diff_eta;
    DeviceArray<Real> interp_fp;    // Interpolation to flux points
    DeviceArray<Real> corr_deriv;   // Correction function derivatives

    // Residual for convergence check
    DeviceArray<Real> residual;
};

// Flux Reconstruction Solver
class FRSolver {
public:
    FRSolver();
    ~FRSolver();

    // Initialize solver with mesh and parameters
    void initialize(const Mesh& mesh, const SimParams& params);

    // Set initial condition
    void setInitialCondition(const std::function<State(Vec2)>& ic_func);

    // Set boundary conditions
    void setBoundaryCondition(int bc_tag, std::shared_ptr<BoundaryCondition> bc);

    // Advance solution by one time step
    void advance(Real dt);

    // Compute time step based on CFL
    Real computeTimeStep() const;

    // Compute residual norm
    Real computeResidual() const;

    // Get current time
    Real currentTime() const { return time_; }

    // Get iteration count
    int iteration() const { return iter_; }

    // Access solution
    const std::vector<ElementSolution>& solution() const { return solution_; }

    // Copy solution from GPU to host
    void copyToHost();

    // Copy solution from host to GPU
    void copyToDevice();

    // Get simulation parameters
    const SimParams& params() const { return params_; }

    // Get mesh
    const Mesh& mesh() const { return *mesh_; }

    // Compute derived quantities (pressure, Mach, etc.)
    void computeDerivedQuantities(std::vector<std::vector<Real>>& pressure,
                                   std::vector<std::vector<Real>>& mach) const;

private:
    // Compute spatial residual (inviscid + viscous)
    void computeResidual_GPU();

    // Compute inviscid flux contribution
    void computeInviscidFlux_GPU();

    // Compute viscous flux contribution
    void computeViscousFlux_GPU();

    // Compute gradients using BR2
    void computeGradients_GPU();

    // Apply boundary conditions
    void applyBoundaryConditions_GPU();

    // Exchange halo data via MPI
    void exchangeHalo();

    // Initialize FR operators
    void initOperators();

    // Initialize geometry on GPU
    void initGeometry();

    // Allocate GPU memory
    void allocateGPUMemory();

    // Data members
    const Mesh* mesh_;
    SimParams params_;
    std::vector<ElementSolution> solution_;

    // FR operators for each element type
    std::map<ElementType, FROperators> operators_;

    // Boundary conditions
    std::map<int, std::shared_ptr<BoundaryCondition>> bc_map_;

    // GPU data
    std::unique_ptr<GPUSolverData> gpu_data_;

    // MPI communicator
    std::shared_ptr<MPICommunicator> mpi_comm_;

    // Time and iteration tracking
    Real time_;
    int iter_;

    // Number of solution points per element
    int n_sp_tri_;
    int n_sp_quad_;

    // Stream for async operations
    std::unique_ptr<CudaStream> stream_;
};

}  // namespace zhijian
