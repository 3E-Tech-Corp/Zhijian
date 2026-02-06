#include "gpu/kernels.hpp"
#include "common/cuda_device.cuh"
#include "physics/euler.hpp"
#include "physics/navierstokes.hpp"
// Removed CUB dependency to avoid version conflicts with Thrust/ROCm

namespace zhijian {
namespace gpu {

// ============================================================================
// Helper device functions
// ============================================================================

__device__ inline void loadState(const Real* U, int idx, int n_sp, State& s) {
    for (int v = 0; v < N_VARS; ++v) {
        s[v] = U[idx * N_VARS + v];
    }
}

__device__ inline void storeState(Real* U, int idx, const State& s) {
    for (int v = 0; v < N_VARS; ++v) {
        U[idx * N_VARS + v] = s[v];
    }
}

// ============================================================================
// Interpolation Kernels
// ============================================================================

__global__ void interpolateToFluxPointsKernel(
    const Real* __restrict__ U_sp,
    Real* __restrict__ U_fp,
    const Real* __restrict__ interp,
    int n_elem, int n_sp, int n_edges, int n_fp_per_edge)
{
    int elem = blockIdx.x;
    int fp = threadIdx.x;  // flux point index within element

    if (elem >= n_elem) return;

    int edge = fp / n_fp_per_edge;
    int fp_local = fp % n_fp_per_edge;

    if (edge >= n_edges) return;

    // Interpolate each variable
    for (int v = 0; v < N_VARS; ++v) {
        Real val = 0.0;
        for (int sp = 0; sp < n_sp; ++sp) {
            val += interp[edge * n_fp_per_edge * n_sp + fp_local * n_sp + sp]
                   * U_sp[elem * n_sp * N_VARS + sp * N_VARS + v];
        }
        U_fp[(elem * n_edges + edge) * n_fp_per_edge * N_VARS + fp_local * N_VARS + v] = val;
    }
}

void interpolateToFluxPoints(const Real* U_sp, Real* U_fp,
                             const Real* interp,
                             int n_elem, int n_sp, int n_edges,
                             int n_fp_per_edge,
                             cudaStream_t stream)
{
    int threads = n_edges * n_fp_per_edge;
    int blocks = n_elem;
    interpolateToFluxPointsKernel<<<blocks, threads, 0, stream>>>(
        U_sp, U_fp, interp, n_elem, n_sp, n_edges, n_fp_per_edge);
}

// ============================================================================
// Gradient Computation Kernels
// ============================================================================

__global__ void computeReferenceGradientsKernel(
    const Real* __restrict__ U,
    Real* __restrict__ dUdxi,
    Real* __restrict__ dUdeta,
    const Real* __restrict__ diff_xi,
    const Real* __restrict__ diff_eta,
    int n_elem, int n_sp)
{
    int elem = blockIdx.x;
    int sp = threadIdx.x;
    int var = threadIdx.y;

    if (elem >= n_elem || sp >= n_sp || var >= N_VARS) return;

    Real sum_xi = 0.0, sum_eta = 0.0;
    for (int k = 0; k < n_sp; ++k) {
        Real Uk = U[elem * n_sp * N_VARS + k * N_VARS + var];
        sum_xi += diff_xi[sp * n_sp + k] * Uk;
        sum_eta += diff_eta[sp * n_sp + k] * Uk;
    }

    int idx = elem * n_sp * N_VARS + sp * N_VARS + var;
    dUdxi[idx] = sum_xi;
    dUdeta[idx] = sum_eta;
}

void computeReferenceGradients(const Real* U, Real* dUdxi, Real* dUdeta,
                                const Real* diff_xi, const Real* diff_eta,
                                int n_elem, int n_sp,
                                cudaStream_t stream)
{
    dim3 threads(n_sp, N_VARS);
    int blocks = n_elem;
    computeReferenceGradientsKernel<<<blocks, threads, 0, stream>>>(
        U, dUdxi, dUdeta, diff_xi, diff_eta, n_elem, n_sp);
}

__global__ void transformGradientsKernel(
    const Real* __restrict__ dUdxi,
    const Real* __restrict__ dUdeta,
    Real* __restrict__ dUdx,
    Real* __restrict__ dUdy,
    const Real* __restrict__ Jinv,
    int n_elem, int n_sp)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total = n_elem * n_sp * N_VARS;

    if (idx >= total) return;

    int elem = idx / (n_sp * N_VARS);
    int sp = (idx / N_VARS) % n_sp;
    // int var = idx % N_VARS;  // unused — index into vars handled via sp offset

    int jinv_idx = elem * n_sp * 4 + sp * 4;
    Real dxidx = Jinv[jinv_idx + 0];
    Real dxidy = Jinv[jinv_idx + 1];
    Real detadx = Jinv[jinv_idx + 2];
    Real detady = Jinv[jinv_idx + 3];

    Real dxi = dUdxi[idx];
    Real deta = dUdeta[idx];

    dUdx[idx] = dxi * dxidx + deta * detadx;
    dUdy[idx] = dxi * dxidy + deta * detady;
}

void transformGradients(const Real* dUdxi, const Real* dUdeta,
                        Real* dUdx, Real* dUdy,
                        const Real* Jinv,
                        int n_elem, int n_sp,
                        cudaStream_t stream)
{
    int total = n_elem * n_sp * N_VARS;
    int threads = 256;
    int blocks = (total + threads - 1) / threads;
    transformGradientsKernel<<<blocks, threads, 0, stream>>>(
        dUdxi, dUdeta, dUdx, dUdy, Jinv, n_elem, n_sp);
}

// ============================================================================
// Flux Computation Kernels
// ============================================================================

__global__ void computeInviscidFluxSPKernel(
    const Real* __restrict__ U,
    Real* __restrict__ Fx,
    Real* __restrict__ Fy,
    Real gamma,
    int n_elem, int n_sp)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total = n_elem * n_sp;

    if (idx >= total) return;

    State s;
    for (int v = 0; v < N_VARS; ++v) {
        s[v] = U[idx * N_VARS + v];
    }

    IdealGas gas(gamma);
    State fx = gas.fluxX(s);
    State fy = gas.fluxY(s);

    for (int v = 0; v < N_VARS; ++v) {
        Fx[idx * N_VARS + v] = fx[v];
        Fy[idx * N_VARS + v] = fy[v];
    }
}

void computeInviscidFluxSP(const Real* U, Real* Fx, Real* Fy,
                            Real gamma,
                            int n_elem, int n_sp,
                            cudaStream_t stream)
{
    int total = n_elem * n_sp;
    int threads = 256;
    int blocks = (total + threads - 1) / threads;
    computeInviscidFluxSPKernel<<<blocks, threads, 0, stream>>>(
        U, Fx, Fy, gamma, n_elem, n_sp);
}

__global__ void computeRiemannFluxKernel(
    const Real* __restrict__ U_L,
    const Real* __restrict__ U_R,
    Real* __restrict__ F_common,
    const Real* __restrict__ normals,
    Real gamma,
    RiemannSolver solver_type,
    int n_faces, int n_fp_per_face)
{
    int face = blockIdx.x;
    int fp = threadIdx.x;

    if (face >= n_faces || fp >= n_fp_per_face) return;

    int idx = face * n_fp_per_face + fp;

    // Load states
    State UL, UR;
    for (int v = 0; v < N_VARS; ++v) {
        UL[v] = U_L[idx * N_VARS + v];
        UR[v] = U_R[idx * N_VARS + v];
    }

    // Load normal
    Real nx = normals[idx * 2 + 0];
    Real ny = normals[idx * 2 + 1];

    // Compute Riemann flux
    State F;
    switch (solver_type) {
        case RiemannSolver::Rusanov:
            F = RiemannSolvers::rusanov(UL, UR, nx, ny, gamma);
            break;
        case RiemannSolver::Roe:
            F = RiemannSolvers::roe(UL, UR, nx, ny, gamma);
            break;
        case RiemannSolver::HLLC:
            F = RiemannSolvers::hllc(UL, UR, nx, ny, gamma);
            break;
    }

    // Store result
    for (int v = 0; v < N_VARS; ++v) {
        F_common[idx * N_VARS + v] = F[v];
    }
}

void computeRiemannFlux(const Real* U_L, const Real* U_R,
                         Real* F_common,
                         const Real* normals,
                         Real gamma,
                         RiemannSolver solver_type,
                         int n_faces, int n_fp_per_face,
                         cudaStream_t stream)
{
    int threads = n_fp_per_face;
    int blocks = n_faces;
    computeRiemannFluxKernel<<<blocks, threads, 0, stream>>>(
        U_L, U_R, F_common, normals, gamma, solver_type, n_faces, n_fp_per_face);
}

// ============================================================================
// Viscous Flux Kernels
// ============================================================================

__global__ void computeViscousFluxSPKernel(
    const Real* __restrict__ U,
    const Real* __restrict__ dUdx,
    const Real* __restrict__ dUdy,
    Real* __restrict__ Fv_x,
    Real* __restrict__ Fv_y,
    Real gamma, Real Re, Real Pr,
    int n_elem, int n_sp)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total = n_elem * n_sp;

    if (idx >= total) return;

    // Load state and gradients
    State s, dsdx, dsdy;
    for (int v = 0; v < N_VARS; ++v) {
        s[v] = U[idx * N_VARS + v];
        dsdx[v] = dUdx[idx * N_VARS + v];
        dsdy[v] = dUdy[idx * N_VARS + v];
    }

    // Compute velocity gradients
    Real rho = s.rho();
    Real u = s.rhou() / rho;
    Real v = s.rhov() / rho;

    Real dudx = (dsdx[1] - u * dsdx[0]) / rho;
    Real dudy = (dsdy[1] - u * dsdy[0]) / rho;
    Real dvdx = (dsdx[2] - v * dsdx[0]) / rho;
    Real dvdy = (dsdy[2] - v * dsdy[0]) / rho;

    // Compute temperature gradient
    IdealGas gas(gamma);
    Real p = gas.pressure(s);
    Real T = p / rho;  // Non-dimensional, R = 1

    Real dpdx = (gamma - 1.0) * (dsdx[3] - 0.5 * dsdx[0] * (u*u + v*v)
                                  - rho * (u * dudx + v * dvdx));
    Real dpdy = (gamma - 1.0) * (dsdy[3] - 0.5 * dsdy[0] * (u*u + v*v)
                                  - rho * (u * dudy + v * dvdy));

    Real dTdx = (dpdx - T * dsdx[0]) / rho;
    Real dTdy = (dpdy - T * dsdy[0]) / rho;

    // Viscosity (Sutherland's law simplified)
    Real mu = pow(T, 0.7) / Re;  // Simplified power law
    Real k = mu * gamma / ((gamma - 1.0) * Pr);

    // Stress tensor
    Real div_v = dudx + dvdy;
    Real lambda = -2.0 / 3.0 * mu;
    Real tau_xx = 2.0 * mu * dudx + lambda * div_v;
    Real tau_yy = 2.0 * mu * dvdy + lambda * div_v;
    Real tau_xy = mu * (dudy + dvdx);

    // Viscous fluxes
    State Fx, Fy;
    Fx[0] = 0.0;
    Fx[1] = tau_xx;
    Fx[2] = tau_xy;
    Fx[3] = u * tau_xx + v * tau_xy + k * dTdx;

    Fy[0] = 0.0;
    Fy[1] = tau_xy;
    Fy[2] = tau_yy;
    Fy[3] = u * tau_xy + v * tau_yy + k * dTdy;

    // Store
    for (int var = 0; var < N_VARS; ++var) {
        Fv_x[idx * N_VARS + var] = Fx[var];
        Fv_y[idx * N_VARS + var] = Fy[var];
    }
}

void computeViscousFluxSP(const Real* U,
                           const Real* dUdx, const Real* dUdy,
                           Real* Fv_x, Real* Fv_y,
                           Real gamma, Real Re, Real Pr,
                           int n_elem, int n_sp,
                           cudaStream_t stream)
{
    int total = n_elem * n_sp;
    int threads = 256;
    int blocks = (total + threads - 1) / threads;
    computeViscousFluxSPKernel<<<blocks, threads, 0, stream>>>(
        U, dUdx, dUdy, Fv_x, Fv_y, gamma, Re, Pr, n_elem, n_sp);
}

// ============================================================================
// Divergence and Correction Kernels
// ============================================================================

__global__ void computeFluxDivergenceKernel(
    const Real* __restrict__ Fx,
    const Real* __restrict__ Fy,
    Real* __restrict__ div_F,
    const Real* __restrict__ diff_xi,
    const Real* __restrict__ diff_eta,
    const Real* __restrict__ Jinv,
    const Real* __restrict__ J,
    int n_elem, int n_sp)
{
    int elem = blockIdx.x;
    int sp = threadIdx.x;
    int var = threadIdx.y;

    if (elem >= n_elem || sp >= n_sp || var >= N_VARS) return;

    // Get inverse Jacobian at this solution point
    int jinv_idx = elem * n_sp * 4 + sp * 4;
    Real dxidx = Jinv[jinv_idx + 0];
    Real dxidy = Jinv[jinv_idx + 1];
    Real detadx = Jinv[jinv_idx + 2];
    Real detady = Jinv[jinv_idx + 3];

    // Compute d/dxi and d/deta of transformed flux
    Real dFdxi = 0.0, dFdeta = 0.0;
    Real dGdxi = 0.0, dGdeta = 0.0;

    for (int k = 0; k < n_sp; ++k) {
        int k_idx = elem * n_sp * N_VARS + k * N_VARS + var;
        Real Fx_k = Fx[k_idx];
        Real Fy_k = Fy[k_idx];

        dFdxi += diff_xi[sp * n_sp + k] * Fx_k;
        dFdeta += diff_eta[sp * n_sp + k] * Fx_k;
        dGdxi += diff_xi[sp * n_sp + k] * Fy_k;
        dGdeta += diff_eta[sp * n_sp + k] * Fy_k;
    }

    // Physical divergence: dF/dx + dG/dy
    // Using chain rule: dF/dx = dF/dxi * dxi/dx + dF/deta * deta/dx
    Real div = (dFdxi * dxidx + dFdeta * detadx) + (dGdxi * dxidy + dGdeta * detady);

    // Scale by Jacobian determinant for proper integration
    Real Jdet = J[elem * n_sp + sp];
    div_F[elem * n_sp * N_VARS + sp * N_VARS + var] = div * Jdet;
}

void computeFluxDivergence(const Real* Fx, const Real* Fy,
                            Real* div_F,
                            const Real* diff_xi, const Real* diff_eta,
                            const Real* Jinv, const Real* J,
                            int n_elem, int n_sp,
                            cudaStream_t stream)
{
    dim3 threads(n_sp, N_VARS);
    int blocks = n_elem;
    computeFluxDivergenceKernel<<<blocks, threads, 0, stream>>>(
        Fx, Fy, div_F, diff_xi, diff_eta, Jinv, J, n_elem, n_sp);
}

__global__ void applyFRCorrectionKernel(
    Real* __restrict__ div_F,
    const Real* __restrict__ F_diff,
    const Real* __restrict__ corr_deriv,
    const Real* __restrict__ J,
    int n_elem, int n_sp, int n_edges, int n_fp_per_edge)
{
    int elem = blockIdx.x;
    int sp = threadIdx.x;
    int var = threadIdx.y;

    if (elem >= n_elem || sp >= n_sp || var >= N_VARS) return;

    Real correction = 0.0;

    // Sum contributions from all edges
    for (int edge = 0; edge < n_edges; ++edge) {
        Real edge_corr = 0.0;
        for (int fp = 0; fp < n_fp_per_edge; ++fp) {
            int fp_idx = (elem * n_edges + edge) * n_fp_per_edge + fp;
            Real Fdiff = F_diff[fp_idx * N_VARS + var];

            // Correction function derivative at this solution point for this edge/fp
            // Simplified: use single correction value per edge
            Real g_deriv = corr_deriv[edge * n_sp + sp];

            edge_corr += Fdiff * g_deriv / n_fp_per_edge;  // Average over flux points
        }
        correction += edge_corr;
    }

    // Add correction to divergence
    Real Jdet = J[elem * n_sp + sp];
    int idx = elem * n_sp * N_VARS + sp * N_VARS + var;
    div_F[idx] += correction * Jdet;
}

void applyFRCorrection(Real* div_F,
                        const Real* F_diff,
                        const Real* corr_deriv,
                        const Real* J,
                        int n_elem, int n_sp, int n_edges, int n_fp_per_edge,
                        cudaStream_t stream)
{
    dim3 threads(n_sp, N_VARS);
    int blocks = n_elem;
    applyFRCorrectionKernel<<<blocks, threads, 0, stream>>>(
        div_F, F_diff, corr_deriv, J, n_elem, n_sp, n_edges, n_fp_per_edge);
}

// ============================================================================
// Boundary Condition Kernels
// ============================================================================

__global__ void applyBoundaryConditionsKernel(
    const Real* __restrict__ U_int,
    Real* __restrict__ U_ghost,
    const int* __restrict__ bc_type,
    const Real* __restrict__ bc_data,
    const Real* __restrict__ normals,
    Real gamma,
    int n_bc_faces, int n_fp_per_face)
{
    int face = blockIdx.x;
    int fp = threadIdx.x;

    if (face >= n_bc_faces || fp >= n_fp_per_face) return;

    int idx = face * n_fp_per_face + fp;
    BCType type = static_cast<BCType>(bc_type[face]);

    // Load interior state
    State U;
    for (int v = 0; v < N_VARS; ++v) {
        U[v] = U_int[idx * N_VARS + v];
    }

    // Load normal
    Real nx = normals[idx * 2 + 0];
    Real ny = normals[idx * 2 + 1];

    IdealGas gas(gamma);
    State ghost;

    switch (type) {
        case BCType::Wall:
        case BCType::SlipWall:
        case BCType::Symmetry: {
            // Reflect velocity
            Real rho = U.rho();
            Real u = U.rhou() / rho;
            Real v = U.rhov() / rho;
            Real p = gas.pressure(U);

            Real vn = u * nx + v * ny;
            Real u_ghost = u - 2.0 * vn * nx;
            Real v_ghost = v - 2.0 * vn * ny;

            ghost = gas.primToConserv(rho, u_ghost, v_ghost, p);
            break;
        }

        case BCType::FarField: {
            // Far-field BC data: [rho_inf, u_inf, v_inf, p_inf, ...]
            int bc_offset = face * 8;  // 8 values per face
            Real rho_inf = bc_data[bc_offset + 0];
            Real u_inf = bc_data[bc_offset + 1];
            Real v_inf = bc_data[bc_offset + 2];
            Real p_inf = bc_data[bc_offset + 3];

            State U_inf = gas.primToConserv(rho_inf, u_inf, v_inf, p_inf);
            Real c_inf = gas.soundSpeed(U_inf);

            // Interior state
            Real rho_i = U.rho();
            Real u_i = U.rhou() / rho_i;
            Real v_i = U.rhov() / rho_i;
            Real c_i = gas.soundSpeed(U);

            Real vn_i = u_i * nx + v_i * ny;
            Real vn_inf = u_inf * nx + v_inf * ny;

            // Riemann invariants
            Real R_plus = vn_i + 2.0 * c_i / (gamma - 1.0);
            Real R_minus = vn_inf - 2.0 * c_inf / (gamma - 1.0);

            Real vn_b = 0.5 * (R_plus + R_minus);
            Real c_b = 0.25 * (gamma - 1.0) * (R_plus - R_minus);

            Real rho_b, u_b, v_b, p_b;

            if (vn_b >= 0) {
                // Outflow
                Real vt_i = u_i * (-ny) + v_i * nx;
                Real s_i = gas.pressure(U) / pow(rho_i, gamma);
                rho_b = pow(c_b * c_b / (gamma * s_i), 1.0 / (gamma - 1.0));
                p_b = rho_b * c_b * c_b / gamma;
                u_b = vn_b * nx - vt_i * ny;
                v_b = vn_b * ny + vt_i * nx;
            } else {
                // Inflow
                Real vt_inf = u_inf * (-ny) + v_inf * nx;
                Real s_inf = p_inf / pow(rho_inf, gamma);
                rho_b = pow(c_b * c_b / (gamma * s_inf), 1.0 / (gamma - 1.0));
                p_b = rho_b * c_b * c_b / gamma;
                u_b = vn_b * nx - vt_inf * ny;
                v_b = vn_b * ny + vt_inf * nx;
            }

            ghost = gas.primToConserv(rho_b, u_b, v_b, p_b);
            break;
        }

        case BCType::Outflow: {
            // Extrapolate density and velocity, use specified pressure
            Real rho = U.rho();
            Real u = U.rhou() / rho;
            Real v = U.rhov() / rho;
            Real p_out = bc_data[face * 8 + 0];  // Specified outlet pressure (8 values per face)
            ghost = gas.primToConserv(rho, u, v, p_out);
            break;
        }

        case BCType::Inflow: {
            // Subsonic inflow: specify total conditions, extrapolate pressure from interior
            int bc_offset = face * 8;
            Real p_total = bc_data[bc_offset + 0];
            Real T_total = bc_data[bc_offset + 1];
            Real dir_x = bc_data[bc_offset + 2];
            Real dir_y = bc_data[bc_offset + 3];

            // Normalize direction
            Real dir_mag = sqrt(dir_x * dir_x + dir_y * dir_y);
            if (dir_mag > 1e-10) {
                dir_x /= dir_mag;
                dir_y /= dir_mag;
            } else {
                dir_x = -nx;  // Default: opposite to outward normal
                dir_y = -ny;
            }

            // Extrapolate static pressure from interior
            Real p_static = gas.pressure(U);

            // Compute Mach number from isentropic relation: p_total/p = (1 + (gamma-1)/2 * M^2)^(gamma/(gamma-1))
            Real pr = p_total / p_static;
            Real exp_inv = (gamma - 1.0) / gamma;
            Real M_sq = 2.0 / (gamma - 1.0) * (pow(pr, exp_inv) - 1.0);
            M_sq = fmax(M_sq, 0.0);  // Ensure non-negative
            Real M = sqrt(M_sq);

            // Compute static temperature: T_total/T = 1 + (gamma-1)/2 * M^2
            Real T_static = T_total / (1.0 + 0.5 * (gamma - 1.0) * M_sq);

            // Compute density and velocity (assuming R = p/(rho*T), with R absorbed into units)
            Real rho_b = p_static / T_static;  // Assumes p = rho * T (R=1 non-dimensional)
            Real c_b = sqrt(gamma * p_static / rho_b);
            Real vel_mag = M * c_b;

            Real u_b = vel_mag * dir_x;
            Real v_b = vel_mag * dir_y;

            ghost = gas.primToConserv(rho_b, u_b, v_b, p_static);
            break;
        }

        default:
            ghost = U;  // Copy interior state
            break;
    }

    // Store ghost state
    for (int var = 0; var < N_VARS; ++var) {
        U_ghost[idx * N_VARS + var] = ghost[var];
    }
}

void applyBoundaryConditions(const Real* U_int, Real* U_ghost,
                              const int* bc_type,
                              const Real* bc_data,
                              const Real* normals,
                              Real gamma,
                              int n_bc_faces, int n_fp_per_face,
                              cudaStream_t stream)
{
    int threads = n_fp_per_face;
    int blocks = n_bc_faces;
    applyBoundaryConditionsKernel<<<blocks, threads, 0, stream>>>(
        U_int, U_ghost, bc_type, bc_data, normals, gamma, n_bc_faces, n_fp_per_face);
}

// ============================================================================
// Time Integration Kernels (SSP-RK3)
// ============================================================================

__global__ void rkStage1Kernel(Real* U1, const Real* U0, const Real* R,
                                int n, Real dt)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    U1[idx] = U0[idx] + dt * R[idx];
}

void rkStage1(Real* U1, const Real* U0, const Real* R,
              int n, Real dt, cudaStream_t stream)
{
    int threads = 256;
    int blocks = (n + threads - 1) / threads;
    rkStage1Kernel<<<blocks, threads, 0, stream>>>(U1, U0, R, n, dt);
}

__global__ void rkStage2Kernel(Real* U2, const Real* U0, const Real* U1,
                                const Real* R, int n, Real dt)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    U2[idx] = 0.75 * U0[idx] + 0.25 * U1[idx] + 0.25 * dt * R[idx];
}

void rkStage2(Real* U2, const Real* U0, const Real* U1, const Real* R,
              int n, Real dt, cudaStream_t stream)
{
    int threads = 256;
    int blocks = (n + threads - 1) / threads;
    rkStage2Kernel<<<blocks, threads, 0, stream>>>(U2, U0, U1, R, n, dt);
}

__global__ void rkStage3Kernel(Real* U, const Real* U0, const Real* U2,
                                const Real* R, int n, Real dt)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    U[idx] = (1.0/3.0) * U0[idx] + (2.0/3.0) * U2[idx] + (2.0/3.0) * dt * R[idx];
}

void rkStage3(Real* U, const Real* U0, const Real* U2, const Real* R,
              int n, Real dt, cudaStream_t stream)
{
    int threads = 256;
    int blocks = (n + threads - 1) / threads;
    rkStage3Kernel<<<blocks, threads, 0, stream>>>(U, U0, U2, R, n, dt);
}

// ============================================================================
// Utility Kernels
// ============================================================================

__global__ void computeTimeStepKernel(
    const Real* __restrict__ U,
    const Real* __restrict__ J,
    const Real* __restrict__ h_min,
    Real* __restrict__ dt_local,
    Real CFL, Real gamma,
    int n_elem, int n_sp)
{
    extern __shared__ Real shared_dt[];

    int elem = blockIdx.x;
    int sp = threadIdx.x;

    if (elem >= n_elem) return;

    Real local_dt = 1e20;

    if (sp < n_sp) {
        int idx = elem * n_sp + sp;

        // Load state
        State s;
        for (int v = 0; v < N_VARS; ++v) {
            s[v] = U[idx * N_VARS + v];
        }

        IdealGas gas(gamma);
        Real wave_speed = gas.maxWaveSpeed(s);
        Real h = h_min[elem];
        Real Jdet = fabs(J[idx]);

        local_dt = CFL * h / (wave_speed * sqrt(Jdet));
    }

    shared_dt[sp] = local_dt;
    __syncthreads();

    // Reduction within block
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (sp < s && sp + s < n_sp) {
            shared_dt[sp] = fmin(shared_dt[sp], shared_dt[sp + s]);
        }
        __syncthreads();
    }

    if (sp == 0) {
        dt_local[elem] = shared_dt[0];
    }
}

// Kernel to find minimum value using parallel reduction
__global__ void reduceMinKernel(const Real* __restrict__ input, Real* __restrict__ output, int n)
{
    extern __shared__ Real shared[];

    int tid = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // Load data into shared memory (use large value for out-of-bounds)
    shared[tid] = (idx < n) ? input[idx] : 1e20;
    __syncthreads();

    // Parallel reduction in shared memory
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            shared[tid] = fmin(shared[tid], shared[tid + s]);
        }
        __syncthreads();
    }

    // Write result for this block
    if (tid == 0) {
        output[blockIdx.x] = shared[0];
    }
}

Real computeTimeStep(const Real* U, const Real* J,
                      const Real* h_min,
                      Real CFL, Real gamma,
                      int n_elem, int n_sp,
                      cudaStream_t stream)
{
    // Compute local dt for each element
    DeviceArray<Real> dt_local(n_elem);

    int threads = n_sp;
    int blocks = n_elem;
    size_t shared_size = threads * sizeof(Real);

    computeTimeStepKernel<<<blocks, threads, shared_size, stream>>>(
        U, J, h_min, dt_local.data(), CFL, gamma, n_elem, n_sp);

    // Find global minimum using manual parallel reduction
    // Phase 1: Reduce n_elem values to fewer blocks
    const int reduce_threads = 256;
    int reduce_blocks = (n_elem + reduce_threads - 1) / reduce_threads;

    DeviceArray<Real> partial_min(reduce_blocks);
    reduceMinKernel<<<reduce_blocks, reduce_threads, reduce_threads * sizeof(Real), stream>>>(
        dt_local.data(), partial_min.data(), n_elem);

    // Phase 2: Further reduce if needed (for very large meshes)
    while (reduce_blocks > 1) {
        int new_blocks = (reduce_blocks + reduce_threads - 1) / reduce_threads;
        DeviceArray<Real> temp_min(new_blocks);
        reduceMinKernel<<<new_blocks, reduce_threads, reduce_threads * sizeof(Real), stream>>>(
            partial_min.data(), temp_min.data(), reduce_blocks);

        // Copy result back
        cudaMemcpyAsync(partial_min.data(), temp_min.data(), new_blocks * sizeof(Real),
                        cudaMemcpyDeviceToDevice, stream);
        reduce_blocks = new_blocks;
    }

    // Copy final minimum to host
    Real min_dt;
    cudaMemcpy(&min_dt, partial_min.data(), sizeof(Real), cudaMemcpyDeviceToHost);

    return min_dt;
}

__global__ void computeResidualNormKernel(
    const Real* __restrict__ R,
    Real* __restrict__ norm_sq,
    int n)
{
    extern __shared__ Real shared[];

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int tid = threadIdx.x;

    Real val = 0.0;
    if (idx < n) {
        val = R[idx] * R[idx];
    }
    shared[tid] = val;
    __syncthreads();

    // Reduction
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            shared[tid] += shared[tid + s];
        }
        __syncthreads();
    }

    if (tid == 0) {
        atomicAdd_double(norm_sq, shared[0]);
    }
}

Real computeResidualNorm(const Real* R, int n_elem, int n_sp, cudaStream_t stream)
{
    int n = n_elem * n_sp * N_VARS;
    int threads = 256;
    int blocks = (n + threads - 1) / threads;

    Real* d_norm_sq;
    cudaMalloc(&d_norm_sq, sizeof(Real));
    cudaMemset(d_norm_sq, 0, sizeof(Real));

    computeResidualNormKernel<<<blocks, threads, threads * sizeof(Real), stream>>>(
        R, d_norm_sq, n);

    Real norm_sq;
    cudaMemcpy(&norm_sq, d_norm_sq, sizeof(Real), cudaMemcpyDeviceToHost);
    cudaFree(d_norm_sq);

    return sqrt(norm_sq);
}

__global__ void setInitialConditionKernel(
    Real* U, Real rho, Real u, Real v, Real p, Real gamma, int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;

    // Convert to conservative variables
    Real rhoE = p / (gamma - 1.0) + 0.5 * rho * (u * u + v * v);

    U[idx * N_VARS + 0] = rho;
    U[idx * N_VARS + 1] = rho * u;
    U[idx * N_VARS + 2] = rho * v;
    U[idx * N_VARS + 3] = rhoE;
}

void setInitialCondition(Real* U, Real rho, Real u, Real v, Real p,
                          Real gamma, int n_elem, int n_sp, cudaStream_t stream)
{
    int n = n_elem * n_sp;
    int threads = 256;
    int blocks = (n + threads - 1) / threads;
    setInitialConditionKernel<<<blocks, threads, 0, stream>>>(
        U, rho, u, v, p, gamma, n);
}

void copySolution(Real* dst, const Real* src, int n, cudaStream_t stream)
{
    cudaMemcpyAsync(dst, src, n * sizeof(Real), cudaMemcpyDeviceToDevice, stream);
}

__global__ void scaleArrayKernel(Real* dst, const Real* src, Real alpha, int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    dst[idx] = alpha * src[idx];
}

void scaleArray(Real* dst, const Real* src, Real alpha, int n, cudaStream_t stream)
{
    int threads = 256;
    int blocks = (n + threads - 1) / threads;
    scaleArrayKernel<<<blocks, threads, 0, stream>>>(dst, src, alpha, n);
}

__global__ void axpyKernel(Real* y, const Real* x, Real alpha, int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    y[idx] += alpha * x[idx];
}

void axpy(Real* y, const Real* x, Real alpha, int n, cudaStream_t stream)
{
    int threads = 256;
    int blocks = (n + threads - 1) / threads;
    axpyKernel<<<blocks, threads, 0, stream>>>(y, x, alpha, n);
}

}  // namespace gpu
}  // namespace zhijian
