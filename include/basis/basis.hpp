#pragma once

#include "common/types.hpp"
#include <vector>
#include <array>

namespace zhijian {

// Gauss-Legendre quadrature points and weights
class GaussLegendre {
public:
    // Get quadrature points and weights for n-point rule
    static void getQuadrature(int n, std::vector<Real>& points, std::vector<Real>& weights);

    // Get Gauss-Lobatto-Legendre points and weights
    static void getGLL(int n, std::vector<Real>& points, std::vector<Real>& weights);
};

// 1D Lagrange basis functions
class Lagrange1D {
public:
    // Evaluate Lagrange basis function j at point x
    // nodes: interpolation nodes in [-1, 1]
    static Real eval(const std::vector<Real>& nodes, int j, Real x);

    // Evaluate derivative of Lagrange basis function j at point x
    static Real evalDeriv(const std::vector<Real>& nodes, int j, Real x);

    // Evaluate all basis functions at a point
    static void evalAll(const std::vector<Real>& nodes, Real x, std::vector<Real>& values);

    // Evaluate all derivatives at a point
    static void evalAllDeriv(const std::vector<Real>& nodes, Real x, std::vector<Real>& derivs);

    // Build Vandermonde matrix (nodes to polynomial coefficients)
    static void buildVandermonde(const std::vector<Real>& nodes,
                                  std::vector<std::vector<Real>>& V);

    // Build differentiation matrix
    static void buildDiffMatrix(const std::vector<Real>& nodes,
                                 std::vector<std::vector<Real>>& D);
};

// FR correction function types
// These define the correction functions g_L and g_R at the element boundaries
class CorrectionFunction {
public:
    // Get correction function based on flux type
    // Returns the value of g(xi) where xi in [-1, 1]
    // For left boundary: g_L(-1) = 1, g_L(1) = 0
    // For right boundary: g_R(-1) = 0, g_R(1) = 1
    static Real evalLeft(FluxType type, int order, Real xi);
    static Real evalRight(FluxType type, int order, Real xi);

    // Evaluate derivative of correction function
    static Real evalLeftDeriv(FluxType type, int order, Real xi);
    static Real evalRightDeriv(FluxType type, int order, Real xi);

    // Get the c parameter for different flux types (Vincent et al.)
    static Real getC(FluxType type, int order);
};

// Solution points for FR method
class SolutionPoints {
public:
    // Get solution point locations in 1D reference element [-1, 1]
    // For FR, typically use Gauss or Gauss-Lobatto points
    static void get1D(int order, std::vector<Real>& points);

    // Get solution points for triangle in reference coordinates
    // Reference triangle: (0,0), (1,0), (0,1)
    static void getTriangle(int order, std::vector<Real>& xi, std::vector<Real>& eta);

    // Get solution points for quadrilateral in reference coordinates
    // Reference quad: [-1,1] x [-1,1]
    static void getQuad(int order, std::vector<Real>& xi, std::vector<Real>& eta);

    // Number of solution points
    static int numPoints1D(int order) { return order + 1; }
    static int numPointsTri(int order) { return (order + 1) * (order + 2) / 2; }
    static int numPointsQuad(int order) { return (order + 1) * (order + 1); }
};

// Flux points (interface points) for FR method
class FluxPoints {
public:
    // Get flux points on a 1D edge (in [-1, 1])
    static void get1D(int order, std::vector<Real>& points);

    // Number of flux points per edge
    static int numPointsPerEdge(int order) { return order + 1; }
};

// FR Operator matrices for an element
class FROperators {
public:
    FROperators() = default;

    // Initialize operators for given element type and polynomial order
    void init(ElementType elem_type, int poly_order, FluxType flux_type);

    // Interpolation from solution points to flux points
    // For each edge of the element
    const std::vector<std::vector<Real>>& interpToFlux(int edge) const {
        return interp_to_flux_[edge];
    }

    // Extrapolation from solution points to a single flux point
    const std::vector<Real>& extrapToFluxPoint(int edge, int fp) const {
        return extrap_to_fp_[edge][fp];
    }

    // Differentiation matrices (d/dxi, d/deta)
    const std::vector<std::vector<Real>>& diffXi() const { return diff_xi_; }
    const std::vector<std::vector<Real>>& diffEta() const { return diff_eta_; }

    // Correction function derivative at solution points
    // correction_deriv_[edge][sp] = dg/dxi or dg/deta at solution point sp
    // for the correction function associated with edge
    const std::vector<std::vector<Real>>& correctionDeriv() const {
        return correction_deriv_;
    }

    // Quadrature weights for integration
    const std::vector<Real>& quadWeights() const { return quad_weights_; }

    // Solution point locations
    const std::vector<Real>& solutionXi() const { return sol_xi_; }
    const std::vector<Real>& solutionEta() const { return sol_eta_; }

    // Flux point locations per edge
    const std::vector<std::vector<Real>>& fluxPointLoc() const { return flux_pt_loc_; }

    // Number of solution points
    int numSolutionPoints() const { return num_sol_pts_; }

    // Number of flux points per edge
    int numFluxPointsPerEdge() const { return num_flux_pts_per_edge_; }

    // Number of edges
    int numEdges() const { return num_edges_; }

private:
    ElementType elem_type_;
    int poly_order_;
    FluxType flux_type_;

    int num_sol_pts_;
    int num_flux_pts_per_edge_;
    int num_edges_;

    // Solution point locations
    std::vector<Real> sol_xi_;
    std::vector<Real> sol_eta_;

    // Flux point locations for each edge
    std::vector<std::vector<Real>> flux_pt_loc_;

    // Interpolation to flux points (per edge)
    std::vector<std::vector<std::vector<Real>>> interp_to_flux_;

    // Extrapolation to individual flux points
    std::vector<std::vector<std::vector<Real>>> extrap_to_fp_;

    // Differentiation matrices
    std::vector<std::vector<Real>> diff_xi_;
    std::vector<std::vector<Real>> diff_eta_;

    // Correction function derivatives at solution points (per edge)
    std::vector<std::vector<Real>> correction_deriv_;

    // Quadrature weights
    std::vector<Real> quad_weights_;

    void initQuad();
    void initTriangle();
};

// Geometry transformation for elements
class ElementGeometry {
public:
    // Compute Jacobian matrix at a point in reference coordinates
    // Returns [dx/dxi, dx/deta; dy/dxi, dy/deta]
    static void jacobian(const std::vector<Vec2>& nodes, ElementType type, int order,
                         Real xi, Real eta,
                         Real& dxdxi, Real& dxdeta, Real& dydxi, Real& dydeta);

    // Compute Jacobian determinant
    static Real jacobianDet(const std::vector<Vec2>& nodes, ElementType type, int order,
                            Real xi, Real eta);

    // Compute inverse Jacobian
    static void jacobianInverse(Real dxdxi, Real dxdeta, Real dydxi, Real dydeta,
                                 Real& dxidx, Real& dxidy, Real& detadx, Real& detady);

    // Map from reference to physical coordinates
    static Vec2 refToPhys(const std::vector<Vec2>& nodes, ElementType type, int order,
                          Real xi, Real eta);

    // Shape functions for geometric mapping
    static void shapeFunc(ElementType type, int order, Real xi, Real eta,
                          std::vector<Real>& N);

    // Shape function derivatives
    static void shapeFuncDeriv(ElementType type, int order, Real xi, Real eta,
                                std::vector<Real>& dNdxi, std::vector<Real>& dNdeta);

    // Face outward normal in reference coordinates
    static Vec2 refFaceNormal(ElementType type, int face_id);

    // Face parametric coordinate to reference coordinate
    static void faceParamToRef(ElementType type, int face_id, Real s,
                                Real& xi, Real& eta);
};

}  // namespace zhijian
