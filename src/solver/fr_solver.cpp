#include "solver/fr_solver.hpp"
#include "gpu/kernels.hpp"
#include "bc/boundary_condition.hpp"
#include "parallel/mpi_comm.hpp"
#include <iostream>
#include <cmath>

namespace zhijian {

// Check for NaN/Inf in GPU array (production version - no verbose output)
static bool checkNaN(const DeviceArray<Real>& arr, const char* name, size_t size) {
    std::vector<Real> host(size);
    const_cast<DeviceArray<Real>&>(arr).copyToHost(host);
    
    for (size_t i = 0; i < size; ++i) {
        if (std::isnan(host[i]) || std::isinf(host[i])) {
            std::cerr << "ERROR: NaN/Inf in " << name << " at index " << i << std::endl;
            return true;
        }
    }
    return false;
}

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
    gpu_data_->Fx_sp.resize(sol_size);
    gpu_data_->Fy_sp.resize(sol_size);

    gpu_data_->U_fp.resize(fp_size);
    gpu_data_->F_fp.resize(fp_size);

    gpu_data_->U0.resize(sol_size);
    gpu_data_->U1.resize(sol_size);

    // Initialize arrays to zero to avoid NaN from uninitialized memory
    gpu_data_->dUdt.fill(0);
    gpu_data_->dUdx.fill(0);
    gpu_data_->dUdy.fill(0);
    gpu_data_->Fx_sp.fill(0);
    gpu_data_->Fy_sp.fill(0);
    gpu_data_->F_fp.fill(0);

    // Upload FR operators to GPU
    const FROperators& ops = operators_.at(ElementType::Quadrilateral);
    int n_sp = ops.numSolutionPoints();

    std::vector<Real> diff_xi_flat, diff_eta_flat;
    
    // Debug: verify differentiation matrix rows sum to ~0
    Real max_row_sum = 0;
    for (int i = 0; i < n_sp; ++i) {
        Real row_sum = 0;
        for (int j = 0; j < n_sp; ++j) {
            diff_xi_flat.push_back(ops.diffXi()[i][j]);
            diff_eta_flat.push_back(ops.diffEta()[i][j]);
            row_sum += ops.diffXi()[i][j];
        }
        max_row_sum = std::max(max_row_sum, std::abs(row_sum));
    }
    // Debug: Diff matrix row sum check (disable for production)
    // std::cout << "Diff matrix max row sum: " << max_row_sum 
    //           << " (should be ~1e-15 for freestream preservation)" << std::endl;
    (void)max_row_sum;  // suppress unused warning

    gpu_data_->diff_xi.resize(diff_xi_flat.size());
    gpu_data_->diff_xi.copyToDevice(diff_xi_flat);

    gpu_data_->diff_eta.resize(diff_eta_flat.size());
    gpu_data_->diff_eta.copyToDevice(diff_eta_flat);

    // Upload interpolation-to-flux-points matrix
    // For quads: 4 edges, each with n_fp flux points, interpolated from n_sp solution points
    int n_edges = 4;
    size_t interp_size = n_edges * n_fp_per_edge * n_sp;
    std::vector<Real> interp_flat(interp_size, 0.0);
    for (int e = 0; e < n_edges; ++e) {
        const auto& interp_e = ops.interpToFlux(e);
        for (int fp = 0; fp < n_fp_per_edge; ++fp) {
            for (int sp_idx = 0; sp_idx < n_sp; ++sp_idx) {
                interp_flat[e * n_fp_per_edge * n_sp + fp * n_sp + sp_idx] =
                    interp_e[fp][sp_idx];
            }
        }
    }
    gpu_data_->interp_fp.resize(interp_flat.size());
    gpu_data_->interp_fp.copyToDevice(interp_flat);

    // Upload correction function derivatives
    const auto& corr = ops.correctionDeriv();
    size_t corr_size = corr.size() * corr[0].size();
    std::vector<Real> corr_flat(corr_size);
    for (size_t i = 0; i < corr.size(); ++i) {
        for (size_t j = 0; j < corr[i].size(); ++j) {
            corr_flat[i * corr[i].size() + j] = corr[i][j];
        }
    }
    gpu_data_->corr_deriv.resize(corr_flat.size());
    gpu_data_->corr_deriv.copyToDevice(corr_flat);

    // Upload face normals and connectivity
    Index n_faces = mesh_->numFaces();
    std::vector<Real> normals_flat(n_faces * 2);
    std::vector<int> face_left_elem(n_faces);
    std::vector<int> face_left_local(n_faces);
    std::vector<int> face_right_elem(n_faces);
    std::vector<int> face_right_local(n_faces);

    for (Index f = 0; f < n_faces; ++f) {
        Vec2 n = mesh_->faceNormal(f);
        normals_flat[2 * f]     = n.x;
        normals_flat[2 * f + 1] = n.y;

        const Face& face = mesh_->face(f);
        face_left_elem[f] = face.left_elem;
        face_left_local[f] = face.left_face;
        face_right_elem[f] = face.right_elem;  // -1 for boundary
        face_right_local[f] = face.right_face;
    }

    gpu_data_->face_normals.resize(normals_flat.size());
    gpu_data_->face_normals.copyToDevice(normals_flat);

    gpu_data_->face_left_elem.resize(n_faces);
    gpu_data_->face_left_elem.copyToDevice(face_left_elem.data(), n_faces);

    gpu_data_->face_left_local.resize(n_faces);
    gpu_data_->face_left_local.copyToDevice(face_left_local.data(), n_faces);

    gpu_data_->face_right_elem.resize(n_faces);
    gpu_data_->face_right_elem.copyToDevice(face_right_elem.data(), n_faces);

    gpu_data_->face_right_local.resize(n_faces);
    gpu_data_->face_right_local.copyToDevice(face_right_local.data(), n_faces);
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

    // Compute inviscid flux (BCs applied internally before Riemann solver)
    computeInviscidFlux_GPU();

    // Compute viscous flux if enabled
    if (params_.viscous) {
        computeGradients_GPU();
        computeViscousFlux_GPU();
    }
}

void FRSolver::computeInviscidFlux_GPU() {
    int n_elem = static_cast<int>(mesh_->numElements());
    int n_sp = n_sp_quad_;
    int n_fp = FluxPoints::numPointsPerEdge(params_.poly_order);
    int n_edges = 4;  // Currently only support quads

    // Check initial solution
    size_t sol_size = n_elem * n_sp * N_VARS;
    if (checkNaN(gpu_data_->U, "U (initial)", sol_size)) return;

    // Compute flux at solution points
    gpu::computeInviscidFluxSP(gpu_data_->U.data(),
                                gpu_data_->Fx_sp.data(),  // x-flux at solution points
                                gpu_data_->Fy_sp.data(),  // y-flux at solution points
                                params_.gamma, n_elem, n_sp, stream_->get());
    stream_->synchronize();
    if (checkNaN(gpu_data_->Fx_sp, "Fx_sp", sol_size)) return;
    if (checkNaN(gpu_data_->Fy_sp, "Fy_sp", sol_size)) return;
    
    // Debug: Jacobian/metric check (disabled for production)
    // Enable by setting DEBUG_FREESTREAM=1
#if defined(DEBUG_FREESTREAM) && DEBUG_FREESTREAM
    {
        size_t jinv_size = n_elem * n_sp * 4;
        std::vector<Real> h_Jinv(jinv_size);
        gpu_data_->Jinv.copyToHost(h_Jinv);
        
        Real max_metric = 0;
        int max_elem = -1, max_sp = -1;
        for (int e = 0; e < n_elem; ++e) {
            for (int sp = 0; sp < n_sp; ++sp) {
                int idx = e * n_sp * 4 + sp * 4;
                for (int m = 0; m < 4; ++m) {
                    if (std::abs(h_Jinv[idx + m]) > max_metric) {
                        max_metric = std::abs(h_Jinv[idx + m]);
                        max_elem = e;
                        max_sp = sp;
                    }
                }
            }
        }
        std::cout << "Jacobian inverse (metrics) max value: " << max_metric 
                  << " at elem " << max_elem << " sp " << max_sp << std::endl;
        
        if (max_elem >= 0) {
            std::cout << "Metrics at elem " << max_elem << ":" << std::endl;
            for (int sp = 0; sp < std::min(4, n_sp); ++sp) {
                int idx = max_elem * n_sp * 4 + sp * 4;
                std::cout << "  sp" << sp << ": dxi/dx=" << h_Jinv[idx+0] 
                          << " dxi/dy=" << h_Jinv[idx+1]
                          << " deta/dx=" << h_Jinv[idx+2]
                          << " deta/dy=" << h_Jinv[idx+3] << std::endl;
            }
        }
        std::cout << "Expected divergence from roundoff: ~" << 6e-8 * max_metric << std::endl;
    }
#endif

    // Interpolate solution to flux points
    gpu::interpolateToFluxPoints(gpu_data_->U.data(), gpu_data_->U_fp.data(),
                                  gpu_data_->interp_fp.data(),
                                  n_elem, n_sp, 4, n_fp, stream_->get());
    stream_->synchronize();
    size_t fp_size = n_elem * 4 * n_fp * N_VARS;
    if (checkNaN(gpu_data_->U_fp, "U_fp", fp_size)) return;
    
    // Debug: U_fp print disabled (flux is now correct)

    // Build BC data for Riemann flux kernel
    int n_faces = static_cast<int>(mesh_->numFaces());
    std::vector<int> bc_type_per_face(n_faces, static_cast<int>(BCType::Interior));
    std::vector<int> bc_face_map(n_faces, -1);
    std::vector<Real> bc_data;

    // Far-field reference state
    Real rho_inf = params_.rho_inf;
    Real c_inf = sqrt(params_.gamma * params_.p_inf / rho_inf);
    Real vel_inf = params_.Mach_inf * c_inf;
    Real u_inf = vel_inf * cos(params_.AoA * M_PI / 180.0);
    Real v_inf = vel_inf * sin(params_.AoA * M_PI / 180.0);

    int bc_count = 0;
    for (Index f = 0; f < static_cast<Index>(n_faces); ++f) {
        const Face& face = mesh_->face(f);
        if (!face.is_boundary) continue;

        bc_type_per_face[f] = static_cast<int>(face.bc_type);
        bc_face_map[f] = bc_count;

        // Get BC parameters from mesh BCInfo
        Real bc_p_static = params_.p_inf;
        if (face.bc_tag > 0 && mesh_->hasBCInfo(face.bc_tag)) {
            const BCInfo& info = mesh_->getBCInfo(face.bc_tag);
            if (info.p_static > 0) bc_p_static = info.p_static;
        }

        // Pack 8 values per boundary face
        switch (face.bc_type) {
            case BCType::FarField:
                bc_data.push_back(rho_inf);
                bc_data.push_back(u_inf);
                bc_data.push_back(v_inf);
                bc_data.push_back(params_.p_inf);
                bc_data.push_back(0); bc_data.push_back(0);
                bc_data.push_back(0); bc_data.push_back(0);
                break;
            case BCType::Outflow:
                bc_data.push_back(bc_p_static);
                bc_data.push_back(0); bc_data.push_back(0); bc_data.push_back(0);
                bc_data.push_back(0); bc_data.push_back(0);
                bc_data.push_back(0); bc_data.push_back(0);
                break;
            default:
                // Wall, SlipWall, Symmetry, etc.
                for (int i = 0; i < 8; ++i) bc_data.push_back(0);
                break;
        }
        bc_count++;
    }

    // Upload BC data to GPU
    DeviceArray<int> d_bc_type(n_faces);
    d_bc_type.copyToDevice(bc_type_per_face.data(), n_faces);

    DeviceArray<int> d_bc_face_map(n_faces);
    d_bc_face_map.copyToDevice(bc_face_map.data(), n_faces);

    DeviceArray<Real> d_bc_data(bc_data.size() > 0 ? bc_data.size() : 1);
    if (!bc_data.empty()) {
        d_bc_data.copyToDevice(bc_data);
    }

    // Allocate temporary arrays
    DeviceArray<Real> F_face(n_faces * n_fp * N_VARS);   // Common (Riemann) flux
    DeviceArray<Real> F_diff_elem(n_elem * n_edges * n_fp * N_VARS);  // Element-indexed F_diff

    // Debug: check normals
    stream_->synchronize();
    size_t norm_size = n_faces * 2;  // One normal per face
    if (checkNaN(gpu_data_->face_normals, "face_normals", norm_size)) return;

    // Step 1: Compute Riemann (common) flux with BC handling
    gpu::computeRiemannFluxWithBC(
        gpu_data_->U_fp.data(),
        F_face.data(),
        gpu_data_->face_left_elem.data(),
        gpu_data_->face_left_local.data(),
        gpu_data_->face_right_elem.data(),
        gpu_data_->face_right_local.data(),
        d_bc_type.data(),
        d_bc_data.data(),
        d_bc_face_map.data(),
        gpu_data_->face_normals.data(),
        params_.gamma, params_.riemann,
        n_faces, n_fp, n_edges, stream_->get());

    // Step 2: Compute F_diff for FR correction
    // This computes (F_common - F_int) for left elements and (-F_common - F_int) for right elements
    // properly handling the normal direction for each element
    gpu::computeFluxDiffForFR(
        gpu_data_->U_fp.data(),
        F_face.data(),
        F_diff_elem.data(),
        gpu_data_->face_left_elem.data(),
        gpu_data_->face_left_local.data(),
        gpu_data_->face_right_elem.data(),
        gpu_data_->face_right_local.data(),
        gpu_data_->face_normals.data(),
        params_.gamma,
        n_faces, n_fp, n_edges, stream_->get());

    // Step 3: Scatter F_common to element-indexed format (for potential future use)
    gpu::scatterFluxToElements(
        F_face.data(),
        gpu_data_->F_fp.data(),
        gpu_data_->face_left_elem.data(),
        gpu_data_->face_left_local.data(),
        gpu_data_->face_right_elem.data(),
        gpu_data_->face_right_local.data(),
        n_faces, n_fp, n_edges, stream_->get());

    stream_->synchronize();
    if (checkNaN(gpu_data_->F_fp, "F_fp (after scatter)", fp_size)) return;

    // Step 4: Compute flux divergence at solution points (from interior flux)
    // dUdt = -∇·F (negative divergence of flux)
    gpu::computeFluxDivergence(gpu_data_->Fx_sp.data(), gpu_data_->Fy_sp.data(),
                                gpu_data_->dUdt.data(),
                                gpu_data_->diff_xi.data(), gpu_data_->diff_eta.data(),
                                gpu_data_->Jinv.data(), gpu_data_->J.data(),
                                n_elem, n_sp, stream_->get());
    stream_->synchronize();
    if (checkNaN(gpu_data_->dUdt, "dUdt (after divergence)", sol_size)) return;

    // Step 5: Apply FR correction using F_diff = F_common - F_int
    // This adds correction: -g'·F_diff/J to the divergence
    gpu::applyFRCorrection(gpu_data_->dUdt.data(), F_diff_elem.data(),
                            gpu_data_->corr_deriv.data(), gpu_data_->J.data(),
                            n_elem, n_sp, n_edges, n_fp, stream_->get());
    stream_->synchronize();
    if (checkNaN(gpu_data_->dUdt, "dUdt (after FR correction)", sol_size)) return;
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
    std::vector<int> bc_elem_idx;      // Element index for each BC face
    std::vector<int> bc_local_face;    // Local face index (0-3) for each BC face
    std::vector<int> bc_types;
    std::vector<Real> bc_data;
    std::vector<Real> normals;

    int n_fp = FluxPoints::numPointsPerEdge(params_.poly_order);
    int n_edges = 4;  // Currently only support quads

    // Far-field reference state (computed once)
    Real rho_inf = params_.rho_inf;
    Real c_inf = sqrt(params_.gamma * params_.p_inf / rho_inf);
    Real vel_inf = params_.Mach_inf * c_inf;
    Real u_inf = vel_inf * cos(params_.AoA * M_PI / 180.0);
    Real v_inf = vel_inf * sin(params_.AoA * M_PI / 180.0);

    Index n_faces = mesh_->numFaces();
    for (Index face_idx = 0; face_idx < n_faces; ++face_idx) {
        const Face& face = mesh_->face(face_idx);
        if (!face.is_boundary) continue;

        // Store element and local face indices for correct U_fp indexing
        bc_elem_idx.push_back(static_cast<int>(face.left_elem));
        bc_local_face.push_back(static_cast<int>(face.left_face));
        bc_types.push_back(static_cast<int>(face.bc_type));

        // Get BC data from mesh BCInfo if available, otherwise use defaults
        Real bc_p_static = params_.p_inf;
        Real bc_p_total = params_.p_inf * 1.2;
        Real bc_T_total = params_.T_inf;
        Real bc_dir_x = 1.0, bc_dir_y = 0.0;
        Real bc_T_wall = 0.0;  // 0 = adiabatic

        if (face.bc_tag > 0 && mesh_->hasBCInfo(face.bc_tag)) {
            const BCInfo& info = mesh_->getBCInfo(face.bc_tag);
            if (info.p_static > 0) bc_p_static = info.p_static;
            if (info.p_total > 0) bc_p_total = info.p_total;
            if (info.T_total > 0) bc_T_total = info.T_total;
            if (info.flow_direction.norm() > 0) {
                bc_dir_x = info.flow_direction.x;
                bc_dir_y = info.flow_direction.y;
            }
            if (info.T_wall > 0) bc_T_wall = info.T_wall;
        }

        // Pack BC-specific data (8 values per face for flexibility)
        switch (face.bc_type) {
            case BCType::FarField:
                // Far-field: rho_inf, u_inf, v_inf, p_inf, 0, 0, 0, 0
                bc_data.push_back(rho_inf);
                bc_data.push_back(u_inf);
                bc_data.push_back(v_inf);
                bc_data.push_back(params_.p_inf);
                bc_data.push_back(0); bc_data.push_back(0);
                bc_data.push_back(0); bc_data.push_back(0);
                break;
            case BCType::Inflow:
                // Inflow: p_total, T_total, dir_x, dir_y, 0, 0, 0, 0
                bc_data.push_back(bc_p_total);
                bc_data.push_back(bc_T_total);
                bc_data.push_back(bc_dir_x);
                bc_data.push_back(bc_dir_y);
                bc_data.push_back(0); bc_data.push_back(0);
                bc_data.push_back(0); bc_data.push_back(0);
                break;
            case BCType::Outflow:
                // Outflow: p_static, 0, 0, 0, 0, 0, 0, 0
                bc_data.push_back(bc_p_static);
                bc_data.push_back(0); bc_data.push_back(0); bc_data.push_back(0);
                bc_data.push_back(0); bc_data.push_back(0);
                bc_data.push_back(0); bc_data.push_back(0);
                break;
            case BCType::Wall:
                // Wall: T_wall (0 = adiabatic), 0, 0, 0, 0, 0, 0, 0
                bc_data.push_back(bc_T_wall);
                bc_data.push_back(0); bc_data.push_back(0); bc_data.push_back(0);
                bc_data.push_back(0); bc_data.push_back(0);
                bc_data.push_back(0); bc_data.push_back(0);
                break;
            default:
                // SlipWall, Symmetry, etc.: no special data needed
                bc_data.push_back(0); bc_data.push_back(0);
                bc_data.push_back(0); bc_data.push_back(0);
                bc_data.push_back(0); bc_data.push_back(0);
                bc_data.push_back(0); bc_data.push_back(0);
                break;
        }

        // Use the CORRECT face normal for this specific face
        Vec2 n = mesh_->faceNormal(face_idx);
        for (int fp = 0; fp < n_fp; ++fp) {
            normals.push_back(n.x);
            normals.push_back(n.y);
        }
    }

    int n_bc_faces = static_cast<int>(bc_types.size());
    if (n_bc_faces == 0) return;

    // Upload to GPU
    DeviceArray<int> d_elem_idx(n_bc_faces);
    d_elem_idx.copyToDevice(bc_elem_idx.data(), n_bc_faces);

    DeviceArray<int> d_local_face(n_bc_faces);
    d_local_face.copyToDevice(bc_local_face.data(), n_bc_faces);

    DeviceArray<int> d_bc_types(n_bc_faces);
    d_bc_types.copyToDevice(bc_types.data(), n_bc_faces);

    DeviceArray<Real> d_bc_data(bc_data.size());
    d_bc_data.copyToDevice(bc_data);

    DeviceArray<Real> d_normals(normals.size());
    d_normals.copyToDevice(normals);

    // Apply boundary conditions (modifies U_fp in place for boundary faces)
    gpu::applyBoundaryConditions(gpu_data_->U_fp.data(),
                                  d_elem_idx.data(),
                                  d_local_face.data(),
                                  d_bc_types.data(),
                                  d_bc_data.data(),
                                  d_normals.data(),
                                  params_.gamma,
                                  n_bc_faces, n_fp, n_edges, stream_->get());
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
    Real min_h = 1e20;
    Real max_wave_speed = 0.0;

    for (Index i = 0; i < mesh_->numElements(); ++i) {
        const Element& elem = mesh_->element(i);
        const FROperators& ops = operators_.at(elem.type);

        // Get minimum element length scale
        Real h = sqrt(mesh_->elementArea(i));
        min_h = std::min(min_h, h);

        // Get maximum wave speed in this element
        for (int sp = 0; sp < ops.numSolutionPoints(); ++sp) {
            const State& U = solution_[i].sol_pts[sp];
            IdealGas gas(params_.gamma);
            Real wave_speed = gas.maxWaveSpeed(U);
            max_wave_speed = std::max(max_wave_speed, wave_speed);

            // Account for polynomial order
            // For FR/DG, stability requires more restrictive CFL at high orders
            // Using (p+1)^2 for more conservative estimate
            Real dt_local = params_.CFL * h / (wave_speed * (params_.poly_order + 1) * (params_.poly_order + 1));
            dt_min = std::min(dt_min, dt_local);
        }
    }

    // Debug output (first call only)
    static bool first_call = true;
    if (first_call) {
        std::cout << "Time step debug:\n"
                  << "  min_h = " << min_h << "\n"
                  << "  max_wave_speed = " << max_wave_speed << "\n"
                  << "  poly_order = " << params_.poly_order << "\n"
                  << "  CFL = " << params_.CFL << "\n"
                  << "  dt_min = " << dt_min << "\n";
        first_call = false;
    }

    return dt_min;
}

Real FRSolver::computeResidual() const {
    // Compute L2 norm of dU/dt (the actual residual) from GPU
    Index n_elem = mesh_->numElements();
    int n_sp = n_sp_quad_;
    size_t size = n_elem * n_sp * N_VARS;

    std::vector<Real> dUdt_host(size);
    gpu_data_->dUdt.copyToHost(dUdt_host);

    Real norm_sq = 0.0;
    for (size_t i = 0; i < size; ++i) {
        norm_sq += dUdt_host[i] * dUdt_host[i];
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
