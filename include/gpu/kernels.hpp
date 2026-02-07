#pragma once

#include "common/types.hpp"
#include <cuda_runtime.h>

namespace zhijian {
namespace gpu {

// ============================================================================
// Solution Point Operations
// ============================================================================

// Interpolate solution from solution points to flux points
// U_sp: solution at solution points [n_elem, n_sp, N_VARS]
// U_fp: solution at flux points [n_elem, n_edges, n_fp_per_edge, N_VARS]
// interp: interpolation matrix [n_fp_per_edge, n_sp]
void interpolateToFluxPoints(const Real* U_sp, Real* U_fp,
                             const Real* interp,
                             int n_elem, int n_sp, int n_edges,
                             int n_fp_per_edge,
                             cudaStream_t stream = 0);

// Compute gradients using differentiation matrices
// U: solution at solution points [n_elem, n_sp, N_VARS]
// dUdxi, dUdeta: gradients in reference coords [n_elem, n_sp, N_VARS]
// diff_xi, diff_eta: differentiation matrices [n_sp, n_sp]
void computeReferenceGradients(const Real* U, Real* dUdxi, Real* dUdeta,
                                const Real* diff_xi, const Real* diff_eta,
                                int n_elem, int n_sp,
                                cudaStream_t stream = 0);

// Transform gradients from reference to physical coordinates
// dUdxi, dUdeta: reference gradients
// dUdx, dUdy: physical gradients
// Jinv: inverse Jacobian [n_elem, n_sp, 4] (dxi/dx, dxi/dy, deta/dx, deta/dy)
void transformGradients(const Real* dUdxi, const Real* dUdeta,
                        Real* dUdx, Real* dUdy,
                        const Real* Jinv,
                        int n_elem, int n_sp,
                        cudaStream_t stream = 0);

// ============================================================================
// Flux Computation
// ============================================================================

// Compute inviscid flux at solution points
// U: conservative variables
// Fx, Fy: fluxes in x and y directions
void computeInviscidFluxSP(const Real* U, Real* Fx, Real* Fy,
                            Real gamma,
                            int n_elem, int n_sp,
                            cudaStream_t stream = 0);

// Compute common (Riemann) flux at flux points
// U_L, U_R: left and right states at flux points
// F_common: common flux
// normals: face normals [n_faces, n_fp, 2]
void computeRiemannFlux(const Real* U_L, const Real* U_R,
                         Real* F_common,
                         const Real* normals,
                         Real gamma,
                         RiemannSolver solver_type,
                         int n_faces, int n_fp_per_face,
                         cudaStream_t stream = 0);

// New: Compute Riemann flux using element-indexed U_fp with integrated BC handling
// U_fp: element-indexed flux point states [n_elem][4][n_fp][N_VARS]
// F_common: output common flux (face-indexed)
// face_*: face connectivity arrays
// bc_type: BC type for each face (Interior=5 for internal faces)
// bc_data: BC parameters (8 values per boundary face)
// bc_face_map: maps global face index to BC data index (-1 for interior)
void computeRiemannFluxWithBC(
    const Real* U_fp,
    Real* F_common,
    const int* face_left_elem,
    const int* face_left_local,
    const int* face_right_elem,
    const int* face_right_local,
    const int* bc_type,
    const Real* bc_data,
    const int* bc_face_map,
    const Real* normals,
    Real gamma,
    RiemannSolver solver_type,
    int n_faces, int n_fp_per_face,
    cudaStream_t stream = 0);

// Compute viscous flux at solution points
// U: conservative variables
// dUdx, dUdy: gradients
// Fv_x, Fv_y: viscous fluxes
void computeViscousFluxSP(const Real* U,
                           const Real* dUdx, const Real* dUdy,
                           Real* Fv_x, Real* Fv_y,
                           Real gamma, Real Re, Real Pr,
                           int n_elem, int n_sp,
                           cudaStream_t stream = 0);

// Compute interior normal flux at flux points (element-indexed)
void computeInviscidFluxAtFP(const Real* U_fp, Real* F_fp,
                              const Real* face_normals,
                              const int* face_left_elem,
                              const int* face_left_local,
                              Real gamma,
                              int n_elem, int n_edges, int n_fp_per_edge,
                              cudaStream_t stream = 0);

// Compute flux difference: F_diff = F_common - F_int
void computeFluxDifference(const Real* F_common, const Real* F_int,
                            Real* F_diff, int n, cudaStream_t stream = 0);

// Compute common viscous flux at flux points (BR2 scheme)
void computeViscousFluxFP(const Real* U_L, const Real* U_R,
                           const Real* dUdx_L, const Real* dUdy_L,
                           const Real* dUdx_R, const Real* dUdy_R,
                           Real* Fv_common,
                           const Real* normals,
                           Real gamma, Real Re, Real Pr,
                           int n_faces, int n_fp_per_face,
                           cudaStream_t stream = 0);

// ============================================================================
// Divergence and Correction
// ============================================================================

// Compute flux divergence at solution points
// Fx, Fy: fluxes at solution points
// dFdx: flux divergence
// diff_xi, diff_eta: differentiation matrices
// Jinv: inverse Jacobian
void computeFluxDivergence(const Real* Fx, const Real* Fy,
                            Real* div_F,
                            const Real* diff_xi, const Real* diff_eta,
                            const Real* Jinv, const Real* J,
                            int n_elem, int n_sp,
                            cudaStream_t stream = 0);

// Scatter face-indexed common flux to element-indexed format
void scatterFluxToElements(
    const Real* F_face,                // Face-indexed common flux
    Real* F_elem,                      // Element-indexed flux (output)
    const int* face_left_elem,
    const int* face_left_local,
    const int* face_right_elem,
    const int* face_right_local,
    int n_faces, int n_fp_per_face,
    cudaStream_t stream = 0);

// Apply FR correction at solution points
// div_F: current divergence (modified in place)
// F_diff: flux difference at flux points (F_common - F_interior)
// corr_deriv: correction function derivatives at solution points
void applyFRCorrection(Real* div_F,
                        const Real* F_diff,
                        const Real* corr_deriv,
                        const Real* J,
                        int n_elem, int n_sp, int n_edges, int n_fp_per_edge,
                        cudaStream_t stream = 0);

// ============================================================================
// Boundary Conditions
// ============================================================================

// Apply boundary conditions at flux points
// U_int: interior solution at boundary flux points
// U_fp: flux point states (element-indexed, modified in place for boundary faces)
// elem_idx: element index for each boundary face
// local_face: local face index (0-3) for each boundary face
// bc_type: boundary condition type for each face
// bc_data: boundary condition parameters
// normals: face outward normals
void applyBoundaryConditions(Real* U_fp,
                              const int* elem_idx,
                              const int* local_face,
                              const int* bc_type,
                              const Real* bc_data,
                              const Real* normals,
                              Real gamma,
                              int n_bc_faces, int n_fp_per_face,
                              cudaStream_t stream = 0);

// ============================================================================
// Time Integration (SSP-RK3)
// ============================================================================

// RK stage 1: U1 = U0 + dt * R
void rkStage1(Real* U1, const Real* U0, const Real* R,
              int n, Real dt,
              cudaStream_t stream = 0);

// RK stage 2: U2 = 0.75*U0 + 0.25*U1 + 0.25*dt*R
void rkStage2(Real* U2, const Real* U0, const Real* U1, const Real* R,
              int n, Real dt,
              cudaStream_t stream = 0);

// RK stage 3: U = (1/3)*U0 + (2/3)*U2 + (2/3)*dt*R
void rkStage3(Real* U, const Real* U0, const Real* U2, const Real* R,
              int n, Real dt,
              cudaStream_t stream = 0);

// Positivity-preserving limiter: enforces rho > 0, p > 0, bounded velocity
void enforcePositivity(Real* U, int n_states, cudaStream_t stream = 0);

// ============================================================================
// Utility Operations
// ============================================================================

// Compute time step based on CFL condition
// Returns minimum dt across all elements
Real computeTimeStep(const Real* U, const Real* J,
                      const Real* h_min,  // minimum element length scale
                      Real CFL, Real gamma,
                      int n_elem, int n_sp,
                      cudaStream_t stream = 0);

// Compute residual norm (L2)
Real computeResidualNorm(const Real* R,
                          int n_elem, int n_sp,
                          cudaStream_t stream = 0);

// Set initial condition
void setInitialCondition(Real* U,
                          Real rho, Real u, Real v, Real p,
                          Real gamma,
                          int n_elem, int n_sp,
                          cudaStream_t stream = 0);

// Copy solution between arrays
void copySolution(Real* dst, const Real* src, int n,
                  cudaStream_t stream = 0);

// Scale array: dst = alpha * src
void scaleArray(Real* dst, const Real* src, Real alpha, int n,
                cudaStream_t stream = 0);

// AXPY: y = alpha * x + y
void axpy(Real* y, const Real* x, Real alpha, int n,
          cudaStream_t stream = 0);

}  // namespace gpu
}  // namespace zhijian
