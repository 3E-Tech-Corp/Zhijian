#include "basis/basis.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace zhijian {

// ============================================================================
// Gauss-Legendre Quadrature
// ============================================================================

void GaussLegendre::getQuadrature(int n, std::vector<Real>& points, std::vector<Real>& weights) {
    points.resize(n);
    weights.resize(n);

    // Precomputed values for common orders
    switch (n) {
        case 1:
            points[0] = 0.0;
            weights[0] = 2.0;
            break;

        case 2:
            points[0] = -1.0 / sqrt(3.0);
            points[1] = 1.0 / sqrt(3.0);
            weights[0] = weights[1] = 1.0;
            break;

        case 3:
            points[0] = -sqrt(3.0/5.0);
            points[1] = 0.0;
            points[2] = sqrt(3.0/5.0);
            weights[0] = weights[2] = 5.0/9.0;
            weights[1] = 8.0/9.0;
            break;

        case 4:
            points[0] = -sqrt(3.0/7.0 + 2.0/7.0*sqrt(6.0/5.0));
            points[1] = -sqrt(3.0/7.0 - 2.0/7.0*sqrt(6.0/5.0));
            points[2] = sqrt(3.0/7.0 - 2.0/7.0*sqrt(6.0/5.0));
            points[3] = sqrt(3.0/7.0 + 2.0/7.0*sqrt(6.0/5.0));
            weights[0] = weights[3] = (18.0 - sqrt(30.0)) / 36.0;
            weights[1] = weights[2] = (18.0 + sqrt(30.0)) / 36.0;
            break;

        case 5:
            points[0] = -sqrt(5.0 + 2.0*sqrt(10.0/7.0)) / 3.0;
            points[1] = -sqrt(5.0 - 2.0*sqrt(10.0/7.0)) / 3.0;
            points[2] = 0.0;
            points[3] = sqrt(5.0 - 2.0*sqrt(10.0/7.0)) / 3.0;
            points[4] = sqrt(5.0 + 2.0*sqrt(10.0/7.0)) / 3.0;
            weights[0] = weights[4] = (322.0 - 13.0*sqrt(70.0)) / 900.0;
            weights[1] = weights[3] = (322.0 + 13.0*sqrt(70.0)) / 900.0;
            weights[2] = 128.0/225.0;
            break;

        case 6:
            points[0] = -0.932469514203152;
            points[1] = -0.661209386466265;
            points[2] = -0.238619186083197;
            points[3] = 0.238619186083197;
            points[4] = 0.661209386466265;
            points[5] = 0.932469514203152;
            weights[0] = weights[5] = 0.171324492379170;
            weights[1] = weights[4] = 0.360761573048139;
            weights[2] = weights[3] = 0.467913934572691;
            break;

        default:
            // Use Newton-Raphson iteration for higher orders
            for (int i = 0; i < n; ++i) {
                // Initial guess
                Real x = cos(M_PI * (i + 0.75) / (n + 0.5));

                // Newton iteration
                for (int iter = 0; iter < 20; ++iter) {
                    Real P0 = 1.0, P1 = x;
                    for (int j = 2; j <= n; ++j) {
                        Real Pn = ((2*j - 1) * x * P1 - (j - 1) * P0) / j;
                        P0 = P1;
                        P1 = Pn;
                    }
                    // P1 = P_n(x), derivative
                    Real dP = n * (x * P1 - P0) / (x * x - 1.0);
                    Real dx = P1 / dP;
                    x -= dx;
                    if (std::abs(dx) < 1e-15) break;
                }

                points[i] = x;

                // Weight
                Real P0 = 1.0, P1 = x;
                for (int j = 2; j <= n; ++j) {
                    Real Pn = ((2*j - 1) * x * P1 - (j - 1) * P0) / j;
                    P0 = P1;
                    P1 = Pn;
                }
                Real dP = n * (x * P1 - P0) / (x * x - 1.0);
                weights[i] = 2.0 / ((1.0 - x * x) * dP * dP);
            }
            break;
    }
}

void GaussLegendre::getGLL(int n, std::vector<Real>& points, std::vector<Real>& weights) {
    points.resize(n);
    weights.resize(n);

    if (n < 2) {
        throw std::invalid_argument("GLL requires at least 2 points");
    }

    // Endpoints are always -1 and 1
    points[0] = -1.0;
    points[n-1] = 1.0;

    if (n == 2) {
        weights[0] = weights[1] = 1.0;
        return;
    }

    // Interior points are roots of P'_{n-1}(x)
    // Use Newton iteration
    for (int i = 1; i < n - 1; ++i) {
        // Initial guess
        Real x = cos(M_PI * i / (n - 1));

        for (int iter = 0; iter < 20; ++iter) {
            // Evaluate P_{n-1} and P'_{n-1}
            Real P0 = 1.0, P1 = x;
            for (int j = 2; j <= n - 1; ++j) {
                Real Pn = ((2*j - 1) * x * P1 - (j - 1) * P0) / j;
                P0 = P1;
                P1 = Pn;
            }
            Real dP = (n - 1) * (x * P1 - P0) / (x * x - 1.0);

            // d^2P/dx^2
            Real d2P = (2.0 * x * dP - (n - 1) * n * P1) / (1.0 - x * x);

            Real dx = dP / d2P;
            x -= dx;
            if (std::abs(dx) < 1e-15) break;
        }
        points[i] = x;
    }

    // Sort points
    std::sort(points.begin(), points.end());

    // Compute weights
    for (int i = 0; i < n; ++i) {
        Real x = points[i];
        Real P0 = 1.0, P1 = x;
        for (int j = 2; j <= n - 1; ++j) {
            Real Pn = ((2*j - 1) * x * P1 - (j - 1) * P0) / j;
            P0 = P1;
            P1 = Pn;
        }
        weights[i] = 2.0 / (n * (n - 1) * P1 * P1);
    }
}

// ============================================================================
// Lagrange Basis Functions
// ============================================================================

Real Lagrange1D::eval(const std::vector<Real>& nodes, int j, Real x) {
    int n = static_cast<int>(nodes.size());
    Real result = 1.0;
    for (int i = 0; i < n; ++i) {
        if (i != j) {
            result *= (x - nodes[i]) / (nodes[j] - nodes[i]);
        }
    }
    return result;
}

Real Lagrange1D::evalDeriv(const std::vector<Real>& nodes, int j, Real x) {
    int n = static_cast<int>(nodes.size());
    Real result = 0.0;

    for (int k = 0; k < n; ++k) {
        if (k == j) continue;

        Real term = 1.0 / (nodes[j] - nodes[k]);
        for (int i = 0; i < n; ++i) {
            if (i != j && i != k) {
                term *= (x - nodes[i]) / (nodes[j] - nodes[i]);
            }
        }
        result += term;
    }
    return result;
}

void Lagrange1D::evalAll(const std::vector<Real>& nodes, Real x, std::vector<Real>& values) {
    int n = static_cast<int>(nodes.size());
    values.resize(n);
    for (int j = 0; j < n; ++j) {
        values[j] = eval(nodes, j, x);
    }
}

void Lagrange1D::evalAllDeriv(const std::vector<Real>& nodes, Real x, std::vector<Real>& derivs) {
    int n = static_cast<int>(nodes.size());
    derivs.resize(n);
    for (int j = 0; j < n; ++j) {
        derivs[j] = evalDeriv(nodes, j, x);
    }
}

void Lagrange1D::buildDiffMatrix(const std::vector<Real>& nodes,
                                  std::vector<std::vector<Real>>& D) {
    int n = static_cast<int>(nodes.size());
    D.resize(n, std::vector<Real>(n));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            D[i][j] = evalDeriv(nodes, j, nodes[i]);
        }
    }
}

// ============================================================================
// Correction Functions
// ============================================================================

Real CorrectionFunction::getC(FluxType type, int order) {
    // c parameter for different schemes (Vincent et al. 2011)
    int p = order;

    switch (type) {
        case FluxType::DG:
            // DG: c = 0
            return 0.0;

        case FluxType::SD:
            // Spectral Difference: c_SD
            return 2.0 * (p + 1) / ((2*p + 1) * p * std::tgamma(p + 1) * std::tgamma(p + 1));

        case FluxType::HU:
            // Huynh's g2 (energy-stable): c = 2*(p+1) / (p * (2p+1) * factorial(p)^2)
            {
                Real fact_p = std::tgamma(p + 1);
                return 2.0 * (p + 1) / (p * (2*p + 1) * fact_p * fact_p);
            }

        case FluxType::GA:
            // Gauss: c_+ from Jameson
            {
                Real fact_p = std::tgamma(p + 1);
                Real fact_2p = std::tgamma(2*p + 1);
                return 2.0 * (p + 1) * fact_p * fact_p / (fact_2p * (2*p + 1));
            }

        default:
            return 0.0;
    }
}

Real CorrectionFunction::evalLeft(FluxType type, int order, Real xi) {
    // g_L(xi) = (-1)^p / 2 * (P_p(xi) - P_{p+1}(xi))
    // where P_n is the Legendre polynomial
    // This is the correction function for the left boundary (xi = -1)

    int p = order;

    // Evaluate Legendre polynomials
    auto legendre = [](int n, Real x) -> Real {
        if (n == 0) return 1.0;
        if (n == 1) return x;
        Real P0 = 1.0, P1 = x;
        for (int i = 2; i <= n; ++i) {
            Real Pn = ((2*i - 1) * x * P1 - (i - 1) * P0) / i;
            P0 = P1;
            P1 = Pn;
        }
        return P1;
    };

    Real Pp = legendre(p, xi);
    Real Pp1 = legendre(p + 1, xi);

    Real sign = (p % 2 == 0) ? 1.0 : -1.0;
    return sign * 0.5 * (Pp - Pp1);
}

Real CorrectionFunction::evalRight(FluxType type, int order, Real xi) {
    // g_R(xi) = g_L(-xi) with sign adjustment
    // g_R(xi) = 1/2 * (P_p(xi) + P_{p+1}(xi))

    int p = order;

    auto legendre = [](int n, Real x) -> Real {
        if (n == 0) return 1.0;
        if (n == 1) return x;
        Real P0 = 1.0, P1 = x;
        for (int i = 2; i <= n; ++i) {
            Real Pn = ((2*i - 1) * x * P1 - (i - 1) * P0) / i;
            P0 = P1;
            P1 = Pn;
        }
        return P1;
    };

    Real Pp = legendre(p, xi);
    Real Pp1 = legendre(p + 1, xi);

    return 0.5 * (Pp + Pp1);
}

Real CorrectionFunction::evalLeftDeriv(FluxType type, int order, Real xi) {
    // Derivative of g_L

    int p = order;

    // Legendre polynomial derivative: P'_n = n/(x^2-1) * (x*P_n - P_{n-1})
    auto legendre = [](int n, Real x) -> Real {
        if (n == 0) return 1.0;
        if (n == 1) return x;
        Real P0 = 1.0, P1 = x;
        for (int i = 2; i <= n; ++i) {
            Real Pn = ((2*i - 1) * x * P1 - (i - 1) * P0) / i;
            P0 = P1;
            P1 = Pn;
        }
        return P1;
    };

    auto legendreDeriv = [&legendre](int n, Real x) -> Real {
        if (n == 0) return 0.0;
        Real Pn = legendre(n, x);
        Real Pn1 = legendre(n - 1, x);
        if (std::abs(x * x - 1.0) < 1e-14) {
            // At endpoints, use recurrence
            return n * (n + 1) / 2.0 * ((x > 0) ? 1.0 : ((n % 2 == 0) ? 1.0 : -1.0));
        }
        return n / (x * x - 1.0) * (x * Pn - Pn1);
    };

    Real dPp = legendreDeriv(p, xi);
    Real dPp1 = legendreDeriv(p + 1, xi);

    Real sign = (p % 2 == 0) ? 1.0 : -1.0;
    return sign * 0.5 * (dPp - dPp1);
}

Real CorrectionFunction::evalRightDeriv(FluxType type, int order, Real xi) {
    int p = order;

    auto legendre = [](int n, Real x) -> Real {
        if (n == 0) return 1.0;
        if (n == 1) return x;
        Real P0 = 1.0, P1 = x;
        for (int i = 2; i <= n; ++i) {
            Real Pn = ((2*i - 1) * x * P1 - (i - 1) * P0) / i;
            P0 = P1;
            P1 = Pn;
        }
        return P1;
    };

    auto legendreDeriv = [&legendre](int n, Real x) -> Real {
        if (n == 0) return 0.0;
        Real Pn = legendre(n, x);
        Real Pn1 = legendre(n - 1, x);
        if (std::abs(x * x - 1.0) < 1e-14) {
            return n * (n + 1) / 2.0 * ((x > 0) ? 1.0 : ((n % 2 == 0) ? 1.0 : -1.0));
        }
        return n / (x * x - 1.0) * (x * Pn - Pn1);
    };

    Real dPp = legendreDeriv(p, xi);
    Real dPp1 = legendreDeriv(p + 1, xi);

    return 0.5 * (dPp + dPp1);
}

// ============================================================================
// Solution Points
// ============================================================================

void SolutionPoints::get1D(int order, std::vector<Real>& points) {
    // Use Gauss-Legendre points for solution points
    std::vector<Real> weights;
    GaussLegendre::getQuadrature(order + 1, points, weights);
}

void SolutionPoints::getTriangle(int order, std::vector<Real>& xi, std::vector<Real>& eta) {
    // Use collapsed coordinate system for triangle
    // Williams-Shunn-Jameson points or similar
    int n_pts = numPointsTri(order);
    xi.resize(n_pts);
    eta.resize(n_pts);

    // Get 1D Gauss points
    std::vector<Real> pts1d, wts1d;
    GaussLegendre::getQuadrature(order + 1, pts1d, wts1d);

    // Map to triangle using collapsed coordinates
    // Triangle: (0,0), (1,0), (0,1)
    int idx = 0;
    for (int j = 0; j <= order; ++j) {
        for (int i = 0; i <= order - j; ++i) {
            // Collapsed coordinate mapping
            Real s = 0.5 * (1.0 + pts1d[i]);
            Real t = 0.5 * (1.0 + pts1d[j]);

            // Map to triangle
            xi[idx] = s * (1.0 - t);
            eta[idx] = t;
            idx++;
        }
    }
}

void SolutionPoints::getQuad(int order, std::vector<Real>& xi, std::vector<Real>& eta) {
    // Tensor product of 1D Gauss points
    std::vector<Real> pts1d, wts1d;
    GaussLegendre::getQuadrature(order + 1, pts1d, wts1d);

    int n = order + 1;
    xi.resize(n * n);
    eta.resize(n * n);

    int idx = 0;
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            xi[idx] = pts1d[i];
            eta[idx] = pts1d[j];
            idx++;
        }
    }
}

// ============================================================================
// Flux Points
// ============================================================================

void FluxPoints::get1D(int order, std::vector<Real>& points) {
    // Use Gauss-Legendre points for flux points on edges
    std::vector<Real> weights;
    GaussLegendre::getQuadrature(order + 1, points, weights);
}

// ============================================================================
// FR Operators
// ============================================================================

void FROperators::init(ElementType elem_type, int poly_order, FluxType flux_type) {
    elem_type_ = elem_type;
    poly_order_ = poly_order;
    flux_type_ = flux_type;

    if (elem_type == ElementType::Quadrilateral) {
        initQuad();
    } else {
        initTriangle();
    }
}

void FROperators::initQuad() {
    int p = poly_order_;
    num_sol_pts_ = (p + 1) * (p + 1);
    num_flux_pts_per_edge_ = p + 1;
    num_edges_ = 4;

    // Get 1D solution points
    std::vector<Real> pts1d, wts1d;
    GaussLegendre::getQuadrature(p + 1, pts1d, wts1d);

    // Build 2D solution points (tensor product)
    sol_xi_.resize(num_sol_pts_);
    sol_eta_.resize(num_sol_pts_);
    quad_weights_.resize(num_sol_pts_);

    int idx = 0;
    for (int j = 0; j <= p; ++j) {
        for (int i = 0; i <= p; ++i) {
            sol_xi_[idx] = pts1d[i];
            sol_eta_[idx] = pts1d[j];
            quad_weights_[idx] = wts1d[i] * wts1d[j];
            idx++;
        }
    }

    // Get 1D flux points
    std::vector<Real> fp1d, fp_wts;
    FluxPoints::get1D(p, fp1d);

    // Flux point locations for each edge
    // Edge 0: bottom (eta = -1)
    // Edge 1: right (xi = 1)
    // Edge 2: top (eta = 1)
    // Edge 3: left (xi = -1)
    flux_pt_loc_.resize(num_edges_);
    for (int e = 0; e < num_edges_; ++e) {
        flux_pt_loc_[e] = fp1d;
    }

    // Build 1D differentiation matrix
    std::vector<std::vector<Real>> D1d;
    Lagrange1D::buildDiffMatrix(pts1d, D1d);

    // Build 2D differentiation matrices (tensor product)
    diff_xi_.resize(num_sol_pts_, std::vector<Real>(num_sol_pts_, 0.0));
    diff_eta_.resize(num_sol_pts_, std::vector<Real>(num_sol_pts_, 0.0));

    for (int j = 0; j <= p; ++j) {
        for (int i = 0; i <= p; ++i) {
            int row = j * (p + 1) + i;
            for (int k = 0; k <= p; ++k) {
                // d/dxi: derivative in i direction
                diff_xi_[row][j * (p + 1) + k] = D1d[i][k];
                // d/deta: derivative in j direction
                diff_eta_[row][k * (p + 1) + i] = D1d[j][k];
            }
        }
    }

    // Build interpolation to flux points
    interp_to_flux_.resize(num_edges_);
    extrap_to_fp_.resize(num_edges_);

    for (int e = 0; e < num_edges_; ++e) {
        interp_to_flux_[e].resize(num_flux_pts_per_edge_,
                                   std::vector<Real>(num_sol_pts_, 0.0));
        extrap_to_fp_[e].resize(num_flux_pts_per_edge_,
                                 std::vector<Real>(num_sol_pts_, 0.0));

        for (int fp = 0; fp < num_flux_pts_per_edge_; ++fp) {
            std::vector<Real> L1d(p + 1);

            // Location on this edge
            Real s = fp1d[fp];  // Parameter along edge

            Real xi_fp, eta_fp;
            if (e == 0) {       // Bottom: eta = -1
                xi_fp = s;
                eta_fp = -1.0;
            } else if (e == 1) { // Right: xi = 1
                xi_fp = 1.0;
                eta_fp = s;
            } else if (e == 2) { // Top: eta = 1
                xi_fp = -s;  // Reverse direction
                eta_fp = 1.0;
            } else {            // Left: xi = -1
                xi_fp = -1.0;
                eta_fp = -s;  // Reverse direction
            }

            // Evaluate 1D basis functions
            std::vector<Real> Lxi(p + 1), Leta(p + 1);
            Lagrange1D::evalAll(pts1d, xi_fp, Lxi);
            Lagrange1D::evalAll(pts1d, eta_fp, Leta);

            // Tensor product for 2D interpolation
            for (int j = 0; j <= p; ++j) {
                for (int i = 0; i <= p; ++i) {
                    int sp = j * (p + 1) + i;
                    Real val = Lxi[i] * Leta[j];
                    interp_to_flux_[e][fp][sp] = val;
                    extrap_to_fp_[e][fp][sp] = val;
                }
            }
        }
    }

    // Build correction function derivatives at solution points
    correction_deriv_.resize(num_edges_);

    for (int e = 0; e < num_edges_; ++e) {
        correction_deriv_[e].resize(num_sol_pts_, 0.0);

        for (int j = 0; j <= p; ++j) {
            for (int i = 0; i <= p; ++i) {
                int sp = j * (p + 1) + i;
                Real xi = sol_xi_[sp];
                Real eta = sol_eta_[sp];

                // Correction function derivative depends on which edge
                if (e == 0 || e == 2) {
                    // Bottom/top edges: correction in eta direction
                    std::vector<Real> Leta(p + 1);
                    Lagrange1D::evalAll(pts1d, eta, Leta);

                    if (e == 0) {
                        // Bottom edge (eta = -1): g_L correction
                        correction_deriv_[e][sp] = CorrectionFunction::evalLeftDeriv(flux_type_, p, eta)
                                                   * Lagrange1D::eval(pts1d, i, xi);
                    } else {
                        // Top edge (eta = 1): g_R correction
                        correction_deriv_[e][sp] = CorrectionFunction::evalRightDeriv(flux_type_, p, eta)
                                                   * Lagrange1D::eval(pts1d, i, xi);
                    }
                } else {
                    // Left/right edges: correction in xi direction
                    if (e == 3) {
                        // Left edge (xi = -1): g_L correction
                        correction_deriv_[e][sp] = CorrectionFunction::evalLeftDeriv(flux_type_, p, xi)
                                                   * Lagrange1D::eval(pts1d, j, eta);
                    } else {
                        // Right edge (xi = 1): g_R correction
                        correction_deriv_[e][sp] = CorrectionFunction::evalRightDeriv(flux_type_, p, xi)
                                                   * Lagrange1D::eval(pts1d, j, eta);
                    }
                }
            }
        }
    }
}

void FROperators::initTriangle() {
    int p = poly_order_;
    num_sol_pts_ = (p + 1) * (p + 2) / 2;
    num_flux_pts_per_edge_ = p + 1;
    num_edges_ = 3;

    // Get solution points for triangle
    SolutionPoints::getTriangle(p, sol_xi_, sol_eta_);

    // Quadrature weights for triangle (approximate)
    quad_weights_.resize(num_sol_pts_);
    Real total_area = 0.5;  // Reference triangle area
    for (int i = 0; i < num_sol_pts_; ++i) {
        quad_weights_[i] = total_area / num_sol_pts_;  // Simple equal weights
    }

    // Get 1D flux points
    std::vector<Real> fp1d, fp_wts;
    FluxPoints::get1D(p, fp1d);

    // Flux point locations for each edge
    // Edge 0: (0,0) to (1,0)
    // Edge 1: (1,0) to (0,1)
    // Edge 2: (0,1) to (0,0)
    flux_pt_loc_.resize(num_edges_);
    for (int e = 0; e < num_edges_; ++e) {
        flux_pt_loc_[e].resize(num_flux_pts_per_edge_);
        for (int fp = 0; fp < num_flux_pts_per_edge_; ++fp) {
            flux_pt_loc_[e][fp] = 0.5 * (1.0 + fp1d[fp]);  // Map [-1,1] to [0,1]
        }
    }

    // Build interpolation to flux points and differentiation matrices
    // (More complex for triangles - use modal basis or Warp & Blend nodes)
    // Simplified implementation here

    interp_to_flux_.resize(num_edges_);
    extrap_to_fp_.resize(num_edges_);
    correction_deriv_.resize(num_edges_);

    for (int e = 0; e < num_edges_; ++e) {
        interp_to_flux_[e].resize(num_flux_pts_per_edge_,
                                   std::vector<Real>(num_sol_pts_, 0.0));
        extrap_to_fp_[e].resize(num_flux_pts_per_edge_,
                                 std::vector<Real>(num_sol_pts_, 0.0));
        correction_deriv_[e].resize(num_sol_pts_, 0.0);

        // Initialize with simple nearest-neighbor interpolation for now
        // Full implementation would use proper triangle basis functions
    }

    // Differentiation matrices for triangle (simplified)
    diff_xi_.resize(num_sol_pts_, std::vector<Real>(num_sol_pts_, 0.0));
    diff_eta_.resize(num_sol_pts_, std::vector<Real>(num_sol_pts_, 0.0));
}

// ============================================================================
// Element Geometry
// ============================================================================

void ElementGeometry::shapeFunc(ElementType type, int order, Real xi, Real eta,
                                 std::vector<Real>& N) {
    if (type == ElementType::Triangle) {
        if (order == 1) {
            // Linear triangle: N = [1-xi-eta, xi, eta]
            N.resize(3);
            N[0] = 1.0 - xi - eta;
            N[1] = xi;
            N[2] = eta;
        } else {
            // Quadratic triangle (6 nodes)
            N.resize(6);
            Real L0 = 1.0 - xi - eta;
            Real L1 = xi;
            Real L2 = eta;
            N[0] = L0 * (2*L0 - 1);
            N[1] = L1 * (2*L1 - 1);
            N[2] = L2 * (2*L2 - 1);
            N[3] = 4 * L0 * L1;
            N[4] = 4 * L1 * L2;
            N[5] = 4 * L2 * L0;
        }
    } else {
        // Quadrilateral
        if (order == 1) {
            // Bilinear quad: N = 0.25*(1+xi_i*xi)*(1+eta_i*eta)
            N.resize(4);
            N[0] = 0.25 * (1 - xi) * (1 - eta);
            N[1] = 0.25 * (1 + xi) * (1 - eta);
            N[2] = 0.25 * (1 + xi) * (1 + eta);
            N[3] = 0.25 * (1 - xi) * (1 + eta);
        } else {
            // Biquadratic quad (9 nodes)
            N.resize(9);
            // Corner nodes
            N[0] = 0.25 * xi * (xi - 1) * eta * (eta - 1);
            N[1] = 0.25 * xi * (xi + 1) * eta * (eta - 1);
            N[2] = 0.25 * xi * (xi + 1) * eta * (eta + 1);
            N[3] = 0.25 * xi * (xi - 1) * eta * (eta + 1);
            // Edge nodes
            N[4] = 0.5 * (1 - xi*xi) * eta * (eta - 1);
            N[5] = 0.5 * xi * (xi + 1) * (1 - eta*eta);
            N[6] = 0.5 * (1 - xi*xi) * eta * (eta + 1);
            N[7] = 0.5 * xi * (xi - 1) * (1 - eta*eta);
            // Center node
            N[8] = (1 - xi*xi) * (1 - eta*eta);
        }
    }
}

void ElementGeometry::shapeFuncDeriv(ElementType type, int order, Real xi, Real eta,
                                      std::vector<Real>& dNdxi, std::vector<Real>& dNdeta) {
    if (type == ElementType::Triangle) {
        if (order == 1) {
            dNdxi.resize(3);
            dNdeta.resize(3);
            dNdxi[0] = -1.0;  dNdeta[0] = -1.0;
            dNdxi[1] = 1.0;   dNdeta[1] = 0.0;
            dNdxi[2] = 0.0;   dNdeta[2] = 1.0;
        } else {
            // Quadratic triangle
            dNdxi.resize(6);
            dNdeta.resize(6);
            Real L0 = 1.0 - xi - eta;
            Real L1 = xi;
            Real L2 = eta;
            dNdxi[0] = -(4*L0 - 1);  dNdeta[0] = -(4*L0 - 1);
            dNdxi[1] = 4*L1 - 1;     dNdeta[1] = 0;
            dNdxi[2] = 0;            dNdeta[2] = 4*L2 - 1;
            dNdxi[3] = 4*(L0 - L1);  dNdeta[3] = -4*L1;
            dNdxi[4] = 4*L2;         dNdeta[4] = 4*L1;
            dNdxi[5] = -4*L2;        dNdeta[5] = 4*(L0 - L2);
        }
    } else {
        // Quadrilateral
        if (order == 1) {
            dNdxi.resize(4);
            dNdeta.resize(4);
            dNdxi[0] = -0.25 * (1 - eta);  dNdeta[0] = -0.25 * (1 - xi);
            dNdxi[1] = 0.25 * (1 - eta);   dNdeta[1] = -0.25 * (1 + xi);
            dNdxi[2] = 0.25 * (1 + eta);   dNdeta[2] = 0.25 * (1 + xi);
            dNdxi[3] = -0.25 * (1 + eta);  dNdeta[3] = 0.25 * (1 - xi);
        } else {
            // Biquadratic quad (9 nodes)
            dNdxi.resize(9);
            dNdeta.resize(9);
            // Corner nodes
            dNdxi[0] = 0.25 * (2*xi - 1) * eta * (eta - 1);
            dNdeta[0] = 0.25 * xi * (xi - 1) * (2*eta - 1);
            dNdxi[1] = 0.25 * (2*xi + 1) * eta * (eta - 1);
            dNdeta[1] = 0.25 * xi * (xi + 1) * (2*eta - 1);
            dNdxi[2] = 0.25 * (2*xi + 1) * eta * (eta + 1);
            dNdeta[2] = 0.25 * xi * (xi + 1) * (2*eta + 1);
            dNdxi[3] = 0.25 * (2*xi - 1) * eta * (eta + 1);
            dNdeta[3] = 0.25 * xi * (xi - 1) * (2*eta + 1);
            // Edge midside nodes
            dNdxi[4] = -xi * eta * (eta - 1);
            dNdeta[4] = 0.5 * (1 - xi*xi) * (2*eta - 1);
            dNdxi[5] = 0.5 * (2*xi + 1) * (1 - eta*eta);
            dNdeta[5] = -xi * (xi + 1) * eta;
            dNdxi[6] = -xi * eta * (eta + 1);
            dNdeta[6] = 0.5 * (1 - xi*xi) * (2*eta + 1);
            dNdxi[7] = 0.5 * (2*xi - 1) * (1 - eta*eta);
            dNdeta[7] = -xi * (xi - 1) * eta;
            // Center node
            dNdxi[8] = -2*xi * (1 - eta*eta);
            dNdeta[8] = -2*eta * (1 - xi*xi);
        }
    }
}

void ElementGeometry::jacobian(const std::vector<Vec2>& nodes, ElementType type, int order,
                                Real xi, Real eta,
                                Real& dxdxi, Real& dxdeta, Real& dydxi, Real& dydeta) {
    std::vector<Real> dNdxi, dNdeta;
    shapeFuncDeriv(type, order, xi, eta, dNdxi, dNdeta);

    dxdxi = dxdeta = dydxi = dydeta = 0.0;
    int n = static_cast<int>(dNdxi.size());
    for (int i = 0; i < n; ++i) {
        dxdxi += dNdxi[i] * nodes[i].x;
        dxdeta += dNdeta[i] * nodes[i].x;
        dydxi += dNdxi[i] * nodes[i].y;
        dydeta += dNdeta[i] * nodes[i].y;
    }
}

Real ElementGeometry::jacobianDet(const std::vector<Vec2>& nodes, ElementType type, int order,
                                   Real xi, Real eta) {
    Real dxdxi, dxdeta, dydxi, dydeta;
    jacobian(nodes, type, order, xi, eta, dxdxi, dxdeta, dydxi, dydeta);
    return dxdxi * dydeta - dxdeta * dydxi;
}

void ElementGeometry::jacobianInverse(Real dxdxi, Real dxdeta, Real dydxi, Real dydeta,
                                       Real& dxidx, Real& dxidy, Real& detadx, Real& detady) {
    Real J = dxdxi * dydeta - dxdeta * dydxi;
    Real invJ = 1.0 / J;
    dxidx = dydeta * invJ;
    dxidy = -dxdeta * invJ;
    detadx = -dydxi * invJ;
    detady = dxdxi * invJ;
}

Vec2 ElementGeometry::refToPhys(const std::vector<Vec2>& nodes, ElementType type, int order,
                                 Real xi, Real eta) {
    std::vector<Real> N;
    shapeFunc(type, order, xi, eta, N);

    Vec2 p(0, 0);
    for (size_t i = 0; i < N.size(); ++i) {
        p.x += N[i] * nodes[i].x;
        p.y += N[i] * nodes[i].y;
    }
    return p;
}

Vec2 ElementGeometry::refFaceNormal(ElementType type, int face_id) {
    if (type == ElementType::Triangle) {
        // Triangle face normals (outward)
        switch (face_id) {
            case 0: return Vec2(0, -1);   // Bottom: (0,0)-(1,0)
            case 1: return Vec2(1, 1) / sqrt(2.0);  // Hypotenuse: (1,0)-(0,1)
            case 2: return Vec2(-1, 0);   // Left: (0,1)-(0,0)
            default: return Vec2(0, 0);
        }
    } else {
        // Quad face normals
        switch (face_id) {
            case 0: return Vec2(0, -1);   // Bottom
            case 1: return Vec2(1, 0);    // Right
            case 2: return Vec2(0, 1);    // Top
            case 3: return Vec2(-1, 0);   // Left
            default: return Vec2(0, 0);
        }
    }
}

void ElementGeometry::faceParamToRef(ElementType type, int face_id, Real s,
                                      Real& xi, Real& eta) {
    // s in [0, 1] parametrizes the face
    if (type == ElementType::Triangle) {
        switch (face_id) {
            case 0: xi = s; eta = 0; break;        // Bottom
            case 1: xi = 1-s; eta = s; break;      // Hypotenuse
            case 2: xi = 0; eta = 1-s; break;      // Left
        }
    } else {
        // Quad: s in [-1, 1]
        switch (face_id) {
            case 0: xi = s; eta = -1; break;       // Bottom
            case 1: xi = 1; eta = s; break;        // Right
            case 2: xi = -s; eta = 1; break;       // Top
            case 3: xi = -1; eta = -s; break;      // Left
        }
    }
}

}  // namespace zhijian
