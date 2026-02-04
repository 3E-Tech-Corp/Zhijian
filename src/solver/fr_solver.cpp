#include "solver/fr_solver.hpp"
#include "gpu/kernels.hpp"
#include "bc/boundary_condition.hpp"
#include "parallel/mpi_comm.hpp"
#include <iostream>
#include <cmath>

namespace zhijian {

FRSolver::FRSolver()
    : mesh_(nullptr), time_(0.0), iter_(0),
      n_sp_tri_(0), n_sp_quad_(0) {}

FRSolver::~FRSolver() = default;

void FRSolver::initialize(const Mesh& mesh, const SimParams& params) {
    mesh_ = &mesh;
    params_ = params;

    // Initialize FR operators
    initOperators();

    // Calculate number of solution points
    n_sp_quad_ = SolutionPoints::numPointsQuad(params.poly_order);
    n_sp_tri_ = SolutionPoints::numPointsTri(params.poly_order);

    // Initialize solution data
    Index n_elem = mesh.numElements();
    solution_.resize(n_elem);

    for (Index i = 0; i < n_elem; ++i) {
        const Element& elem = mesh.element(i);
        int n_sp = (elem.type == ElementType::Quadrilateral) ? n_sp_quad_ : n_sp_tri_;
        int n_edges = (elem.type == ElementType::Quadrilateral) ? 4 : 3;
        int n_fp = FluxPoints::numPointsPerEdge(params.poly_order);

        solution_[i].sol_pts.resize(n_sp);
        solution_[i].flux_pts.resize(n_edges * n_fp);
        solution_[i].dUdx.resize(n_sp);
        solution_[i].dUdy.resize(n_sp);
    }

    // Allocate GPU memory (must precede initGeometry which writes to gpu_data_)
    allocateGPUMemory();

    // Initialize geometry on GPU
    initGeometry();

    // Create CUDA stream
    stream_ = std::make_unique<CudaStream>();

    std::cout << "FR Solver initialized:" << std::endl;
    std::cout << "  Polynomial order: P" << params.poly_order << std::endl;
    std::cout << "  Elements: " << n_elem << std::endl;
    std::cout << "  Solution points (quad): " << n_sp_quad_ << std::endl;
    std::cout << "  Solution points (tri): " << n_sp_tri_ << std::endl;
    std::cout << "  Flux type: " << static_cast<int>(params.flux_type) << std::endl;
    std::cout << "  Riemann solver: " << static_cast<int>(params.riemann) << std::endl;
}

void FRSolver::initOperators() {
    // Initialize operators for quadrilateral elements
    operators_[ElementType::Quadrilateral].init(
        ElementType::Quadrilateral, params_.poly_order, params_.flux_type);

    // Initialize operators for triangular elements
    operators_[ElementType::Triangle].init(
        ElementType::Triangle, params_.poly_order, params_.flux_type);
}

void FRSolver::initGeometry() {
    // Precompute geometric quantities at solution points for all elements
    Index n_elem = mesh_->numElements();

    // This will be uploaded to GPU
    std::vector<Real> J_host, Jinv_host;

    for (Index i = 0; i < n_elem; ++i) {
        const Element& elem = mesh_->element(i);
        const FROperators& ops = operators_.at(elem.type);

        // Get element nodes
        std::vector<Vec2> nodes(elem.node_ids.size());
        for (size_t j = 0; j < elem.node_ids.size(); ++j) {
            nodes[j] = mesh_->node(elem.node_ids[j]);
        }

        // Compute Jacobian at each solution point
        for (int sp = 0; sp < ops.numSolutionPoints(); ++sp) {
            Real xi = ops.solutionXi()[sp];
            Real eta = ops.solutionEta()[sp];

            Real dxdxi, dxdeta, dydxi, dydeta;
            ElementGeometry::jacobian(nodes, elem.type, elem.order,
                                       xi, eta, dxdxi, dxdeta, dydxi, dydeta);

            Real Jdet = dxdxi * dydeta - dxdeta * dydxi;
            J_host.push_back(Jdet);

            // Inverse Jacobian
            Real invJ = 1.0 / Jdet;
            Jinv_host.push_back(dydeta * invJ);   // dxi/dx
            Jinv_host.push_back(-dxdeta * invJ);  // dxi/dy
            Jinv_host.push_back(-dydxi * invJ);   // deta/dx
            Jinv_host.push_back(dxdxi * invJ);    // deta/dy
        }
    }

    // Upload to GPU
    gpu_data_->J.resize(J_host.size());
    gpu_data_->J.copyToDevice(J_host);

    gpu_data_->Jinv.resize(Jinv_host.size());
    gpu_data_->Jinv.copyToDevice(Jinv_host);
}

void FRSolver::allocateGPUMemory() {
    gpu_data_ = std::make_unique<GPUSolverData>();

    Index n_elem = mesh_->numElements();
    int max_sp = std::max(n_sp_quad_, n_sp_tri_);
    int n_fp_per_edge = FluxPoints::numPointsPerEdge(params_.poly_order);

    // Assuming all quads for simplicity (adjust for mixed meshes)
    size_t sol_size = n_elem * max_sp * N_VARS;
    size_t fp_size = n_elem * 4 * n_fp_per_edge * N_VARS;

    gpu_data_->U.resize(sol_size);
    gpu_data_->dUdt.resize(sol_size);
    gpu_data_->dUdx.resize(sol_size);
    gpu_data_->dUdy.resize(sol_size);

    gpu_data_->U_fp.resize(fp_size);
    gpu_data_->F_fp.resize(fp_size);

    gpu_data_->U0.resize(sol_size);
    gpu_data_->U1.resize(sol_size);

    // Upload FR operators to GPU
    const FROperators& ops = operators_.at(ElementType::Quadrilateral);
    int n_sp = ops.numSolutionPoints();

    std::vector<Real> diff_xi_flat, diff_eta_flat;
    for (int i = 0; i < n_sp; ++i) {
        for (int j = 0; j < n_sp; ++j) {
            diff_xi_flat.push_back(ops.diffXi()[i][j]);
            diff_eta_flat.push_back(ops.diffEta()[i][j]);
        }
    }

    gpu_data_->diff_xi.resize(diff_xi_flat.size());
    gpu_data_->diff_xi.copyToDevice(diff_xi_flat);

    gpu_data_->diff_eta.resize(diff_eta_flat.size());
    gpu_data_->diff_eta.copyToDevice(diff_eta_flat);
}

void FRSolver::setInitialCondition(const std::function<State(Vec2)>& ic_func) {
    Index n_elem = mesh_->numElements();

    for (Index i = 0; i < n_elem; ++i) {
        const Element& elem = mesh_->element(i);
        const FROperators& ops = operators_.at(elem.type);

        // Get element nodes
        std::vector<Vec2> nodes(elem.node_ids.size());
        for (size_t j = 0; j < elem.node_ids.size(); ++j) {
            nodes[j] = mesh_->node(elem.node_ids[j]);
        }

        // Set solution at each solution point
        for (int sp = 0; sp < ops.numSolutionPoints(); ++sp) {
            Real xi = ops.solutionXi()[sp];
            Real eta = ops.solutionEta()[sp];

            // Map to physical coordinates
            Vec2 pos = ElementGeometry::refToPhys(nodes, elem.type, elem.order, xi, eta);

            // Evaluate initial condition
            solution_[i].sol_pts[sp] = ic_func(pos);
        }
    }

    // Copy to GPU
    copyToDevice();

    time_ = 0.0;
    iter_ = 0;
}

void FRSolver::setBoundaryCondition(int bc_tag, std::shared_ptr<BoundaryCondition> bc) {
    bc_map_[bc_tag] = bc;
}

void FRSolver::advance(Real dt) {
    // SSP-RK3 time integration
    int n_elem = static_cast<int>(mesh_->numElements());
    int n_sp = n_sp_quad_;  // Assuming all quads
    int n = n_elem * n_sp * N_VARS;

    // Stage 1: U1 = U0 + dt * R(U0)
    gpu::copySolution(gpu_data_->U0.data(), gpu_data_->U.data(), n, stream_->get());
    computeResidual_GPU();
    gpu::rkStage1(gpu_data_->U.data(), gpu_data_->U0.data(),
                  gpu_data_->dUdt.data(), n, dt, stream_->get());

    // Stage 2: U2 = 0.75*U0 + 0.25*U1 + 0.25*dt*R(U1)
    gpu::copySolution(gpu_data_->U1.data(), gpu_data_->U.data(), n, stream_->get());
    computeResidual_GPU();
    gpu::rkStage2(gpu_data_->U.data(), gpu_data_->U0.data(), gpu_data_->U1.data(),
                  gpu_data_->dUdt.data(), n, dt, stream_->get());

    // Stage 3: U = (1/3)*U0 + (2/3)*U2 + (2/3)*dt*R(U2)
    computeResidual_GPU();
    gpu::rkStage3(gpu_data_->U.data(), gpu_data_->U0.data(), gpu_data_->U.data(),
                  gpu_data_->dUdt.data(), n, dt, stream_->get());

    // Synchronize
    stream_->synchronize();

    time_ += dt;
    iter_++;
}

void FRSolver::computeResidual_GPU() {
    // Exchange halo data for MPI
    if (mpi_comm_) {
        exchangeHalo();
    }

    // Compute inviscid flux
    computeInviscidFlux_GPU();

    // Compute viscous flux if enabled
    if (params_.viscous) {
        computeGradients_GPU();
        computeViscousFlux_GPU();
    }

    // Apply boundary conditions
    applyBoundaryConditions_GPU();
}

void FRSolver::computeInviscidFlux_GPU() {
    int n_elem = static_cast<int>(mesh_->numElements());
    int n_sp = n_sp_quad_;
    int n_fp = FluxPoints::numPointsPerEdge(params_.poly_order);

    // Compute flux at solution points
    gpu::computeInviscidFluxSP(gpu_data_->U.data(),
                                gpu_data_->dUdt.data(),  // Temporary storage for Fx
                                gpu_data_->dUdx.data(),  // Temporary storage for Fy
                                params_.gamma, n_elem, n_sp, stream_->get());

    // Interpolate solution to flux points
    gpu::interpolateToFluxPoints(gpu_data_->U.data(), gpu_data_->U_fp.data(),
                                  gpu_data_->interp_fp.data(),
                                  n_elem, n_sp, 4, n_fp, stream_->get());

    // Compute Riemann flux at interfaces
    int n_faces = static_cast<int>(mesh_->numFaces());
    gpu::computeRiemannFlux(gpu_data_->U_fp.data(),    // Left states
                            gpu_data_->U_fp.data(),    // Right states (needs proper indexing)
                            gpu_data_->F_fp.data(),
                            gpu_data_->face_normals.data(),
                            params_.gamma, params_.riemann,
                            n_faces, n_fp, stream_->get());

    // Compute flux divergence at solution points
    gpu::computeFluxDivergence(gpu_data_->dUdt.data(), gpu_data_->dUdx.data(),
                                gpu_data_->dUdt.data(),
                                gpu_data_->diff_xi.data(), gpu_data_->diff_eta.data(),
                                gpu_data_->Jinv.data(), gpu_data_->J.data(),
                                n_elem, n_sp, stream_->get());

    // Apply FR correction
    gpu::applyFRCorrection(gpu_data_->dUdt.data(), gpu_data_->F_fp.data(),
                            gpu_data_->corr_deriv.data(), gpu_data_->J.data(),
                            n_elem, n_sp, 4, n_fp, stream_->get());
}

void FRSolver::computeGradients_GPU() {
    int n_elem = static_cast<int>(mesh_->numElements());
    int n_sp = n_sp_quad_;

    // Compute reference gradients
    DeviceArray<Real> dUdxi(n_elem * n_sp * N_VARS);
    DeviceArray<Real> dUdeta(n_elem * n_sp * N_VARS);

    gpu::computeReferenceGradients(gpu_data_->U.data(),
                                    dUdxi.data(), dUdeta.data(),
                                    gpu_data_->diff_xi.data(),
                                    gpu_data_->diff_eta.data(),
                                    n_elem, n_sp, stream_->get());

    // Transform to physical gradients
    gpu::transformGradients(dUdxi.data(), dUdeta.data(),
                            gpu_data_->dUdx.data(), gpu_data_->dUdy.data(),
                            gpu_data_->Jinv.data(),
                            n_elem, n_sp, stream_->get());
}

void FRSolver::computeViscousFlux_GPU() {
    int n_elem = static_cast<int>(mesh_->numElements());
    int n_sp = n_sp_quad_;

    // Temporary arrays for viscous flux
    DeviceArray<Real> Fv_x(n_elem * n_sp * N_VARS);
    DeviceArray<Real> Fv_y(n_elem * n_sp * N_VARS);

    gpu::computeViscousFluxSP(gpu_data_->U.data(),
                               gpu_data_->dUdx.data(), gpu_data_->dUdy.data(),
                               Fv_x.data(), Fv_y.data(),
                               params_.gamma, params_.Re, params_.Pr,
                               n_elem, n_sp, stream_->get());

    // Subtract viscous flux divergence from residual
    // (viscous terms have negative sign in NS equations)
}

void FRSolver::applyBoundaryConditions_GPU() {
    // Count boundary faces and prepare data
    std::vector<int> bc_types;
    std::vector<Real> bc_data;
    std::vector<Real> normals;

    int n_fp = FluxPoints::numPointsPerEdge(params_.poly_order);

    for (const auto& face : mesh_->faces()) {
        if (face.is_boundary) {
            bc_types.push_back(static_cast<int>(face.bc_type));

            // Add BC-specific data
            if (face.bc_type == BCType::FarField) {
                // Far-field: store rho_inf, u_inf, v_inf, p_inf
                Real rho = params_.rho_inf;
                Real c = sqrt(params_.gamma * params_.p_inf / rho);
                Real vel = params_.Mach_inf * c;
                Real u = vel * cos(params_.AoA * M_PI / 180.0);
                Real v = vel * sin(params_.AoA * M_PI / 180.0);
                bc_data.push_back(rho);
                bc_data.push_back(u);
                bc_data.push_back(v);
                bc_data.push_back(params_.p_inf);
            } else {
                // Default: zeros
                bc_data.push_back(0);
                bc_data.push_back(0);
                bc_data.push_back(0);
                bc_data.push_back(0);
            }

            // Add face normals for each flux point
            Vec2 n = mesh_->faceNormal(0);  // Would need proper face indexing
            for (int fp = 0; fp < n_fp; ++fp) {
                normals.push_back(n.x);
                normals.push_back(n.y);
            }
        }
    }

    int n_bc_faces = static_cast<int>(bc_types.size());
    if (n_bc_faces == 0) return;

    // Upload to GPU
    DeviceArray<int> d_bc_types(n_bc_faces);
    d_bc_types.copyToDevice(bc_types.data(), n_bc_faces);

    DeviceArray<Real> d_bc_data(bc_data.size());
    d_bc_data.copyToDevice(bc_data);

    DeviceArray<Real> d_normals(normals.size());
    d_normals.copyToDevice(normals);

    // Apply boundary conditions
    gpu::applyBoundaryConditions(gpu_data_->U_fp.data(),  // Interior states
                                  gpu_data_->U_fp.data(),  // Ghost states (output)
                                  d_bc_types.data(),
                                  d_bc_data.data(),
                                  d_normals.data(),
                                  params_.gamma,
                                  n_bc_faces, n_fp, stream_->get());
}

void FRSolver::exchangeHalo() {
    // MPI halo exchange (to be implemented)
    if (mpi_comm_) {
        // Exchange flux point data between partitions
    }
}

Real FRSolver::computeTimeStep() const {
    // Compute time step based on CFL condition
    Real dt_min = 1e20;

    for (Index i = 0; i < mesh_->numElements(); ++i) {
        const Element& elem = mesh_->element(i);
        const FROperators& ops = operators_.at(elem.type);

        // Get minimum element length scale
        Real h = sqrt(mesh_->elementArea(i));

        // Get maximum wave speed in this element
        for (int sp = 0; sp < ops.numSolutionPoints(); ++sp) {
            const State& U = solution_[i].sol_pts[sp];
            IdealGas gas(params_.gamma);
            Real wave_speed = gas.maxWaveSpeed(U);

            // Account for polynomial order
            Real dt_local = params_.CFL * h / (wave_speed * (params_.poly_order + 1));
            dt_min = std::min(dt_min, dt_local);
        }
    }

    return dt_min;
}

Real FRSolver::computeResidual() const {
    Real norm_sq = 0.0;

    for (Index i = 0; i < mesh_->numElements(); ++i) {
        for (const auto& s : solution_[i].sol_pts) {
            for (int v = 0; v < N_VARS; ++v) {
                norm_sq += s[v] * s[v];
            }
        }
    }

    return sqrt(norm_sq);
}

void FRSolver::copyToHost() {
    // Copy solution from GPU to host
    Index n_elem = mesh_->numElements();
    int n_sp = n_sp_quad_;

    std::vector<Real> U_host(n_elem * n_sp * N_VARS);
    gpu_data_->U.copyToHost(U_host);

    // Unpack into solution structure
    int idx = 0;
    for (Index i = 0; i < n_elem; ++i) {
        for (int sp = 0; sp < n_sp; ++sp) {
            for (int v = 0; v < N_VARS; ++v) {
                solution_[i].sol_pts[sp][v] = U_host[idx++];
            }
        }
    }
}

void FRSolver::copyToDevice() {
    // Copy solution from host to GPU
    Index n_elem = mesh_->numElements();
    int n_sp = n_sp_quad_;

    std::vector<Real> U_host(n_elem * n_sp * N_VARS);

    int idx = 0;
    for (Index i = 0; i < n_elem; ++i) {
        for (int sp = 0; sp < n_sp; ++sp) {
            for (int v = 0; v < N_VARS; ++v) {
                U_host[idx++] = solution_[i].sol_pts[sp][v];
            }
        }
    }

    gpu_data_->U.copyToDevice(U_host);
}

void FRSolver::computeDerivedQuantities(std::vector<std::vector<Real>>& pressure,
                                         std::vector<std::vector<Real>>& mach) const {
    Index n_elem = mesh_->numElements();
    pressure.resize(n_elem);
    mach.resize(n_elem);

    IdealGas gas(params_.gamma);

    for (Index i = 0; i < n_elem; ++i) {
        int n_sp = static_cast<int>(solution_[i].sol_pts.size());
        pressure[i].resize(n_sp);
        mach[i].resize(n_sp);

        for (int sp = 0; sp < n_sp; ++sp) {
            const State& U = solution_[i].sol_pts[sp];
            pressure[i][sp] = gas.pressure(U);
            mach[i][sp] = gas.machNumber(U);
        }
    }
}

}  // namespace zhijian
