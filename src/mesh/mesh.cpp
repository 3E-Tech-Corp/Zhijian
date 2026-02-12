#include "mesh/mesh.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <map>

namespace zhijian {

Index Mesh::numBoundaryFaces() const {
    Index count = 0;
    for (const auto& face : faces_) {
        if (face.is_boundary) ++count;
    }
    return count;
}

void Mesh::buildFaces() {
    faces_.clear();

    // Map from edge (sorted node pair) to face info
    // For internal edges, we'll have two elements sharing the edge
    std::map<std::pair<Index, Index>, std::vector<std::pair<LocalIndex, LocalIndex>>> edge_to_elem;

    // Iterate over all elements and their edges
    for (Index elem_id = 0; elem_id < static_cast<Index>(elements_.size()); ++elem_id) {
        const Element& elem = elements_[elem_id];
        int n_vertices = elem.numVertices();

        for (int face_id = 0; face_id < n_vertices; ++face_id) {
            // Get the two vertices of this edge
            Index v0 = elem.node_ids[face_id];
            Index v1 = elem.node_ids[(face_id + 1) % n_vertices];

            // Create sorted edge key
            auto edge_key = std::make_pair(std::min(v0, v1), std::max(v0, v1));

            // Store element and local face index
            edge_to_elem[edge_key].push_back({static_cast<LocalIndex>(elem_id),
                                               static_cast<LocalIndex>(face_id)});
        }
    }

    // Create faces from edge map
    for (const auto& [edge_key, elems] : edge_to_elem) {
        Face face;

        if (elems.size() == 1) {
            // Boundary face
            face.left_elem = elems[0].first;
            face.left_face = elems[0].second;
            face.right_elem = -1;
            face.right_face = -1;
            face.is_boundary = true;
            face.bc_type = BCType::FarField;  // Default to FarField (safer for external boundaries)
            face.bc_tag = 0;
            face.mpi_rank = -1;
        } else if (elems.size() == 2) {
            // Interior face
            face.left_elem = elems[0].first;
            face.left_face = elems[0].second;
            face.right_elem = elems[1].first;
            face.right_face = elems[1].second;
            face.is_boundary = false;
            face.bc_type = BCType::Interior;
            face.bc_tag = -1;
            face.mpi_rank = -1;
        } else {
            std::cerr << "Warning: Edge shared by " << elems.size()
                      << " elements (non-manifold mesh)" << std::endl;
            continue;
        }

        LocalIndex face_idx = static_cast<LocalIndex>(faces_.size());
        faces_.push_back(face);

        // Update element face IDs
        elements_[face.left_elem].face_ids.push_back(face_idx);
        if (!face.is_boundary) {
            elements_[face.right_elem].face_ids.push_back(face_idx);
        }
    }
}

void Mesh::computeGeometry() {
    Index n_elem = numElements();
    Index n_face = numFaces();

    elem_centroids_.resize(n_elem);
    elem_areas_.resize(n_elem);
    face_lengths_.resize(n_face);
    face_normals_.resize(n_face);
    face_centroids_.resize(n_face);

    // Compute element centroids and areas
    for (Index i = 0; i < n_elem; ++i) {
        elem_centroids_[i] = elementCentroid(i);
        elem_areas_[i] = elementArea(i);
    }

    // Compute face geometry
    for (Index i = 0; i < n_face; ++i) {
        face_lengths_[i] = faceLength(i);
        face_normals_[i] = faceNormal(i);
        face_centroids_[i] = faceCentroid(i);
    }
}

Vec2 Mesh::elementCentroid(Index elem_id) const {
    const Element& elem = elements_[elem_id];
    Vec2 centroid(0, 0);

    int n_vertices = elem.numVertices();
    for (int i = 0; i < n_vertices; ++i) {
        centroid = centroid + nodes_[elem.node_ids[i]];
    }
    centroid = centroid / static_cast<Real>(n_vertices);

    return centroid;
}

Real Mesh::elementArea(Index elem_id) const {
    const Element& elem = elements_[elem_id];

    if (elem.type == ElementType::Triangle) {
        // Shoelace formula for triangle
        const Vec2& p0 = nodes_[elem.node_ids[0]];
        const Vec2& p1 = nodes_[elem.node_ids[1]];
        const Vec2& p2 = nodes_[elem.node_ids[2]];

        return 0.5 * std::abs((p1.x - p0.x) * (p2.y - p0.y) -
                               (p2.x - p0.x) * (p1.y - p0.y));
    } else {
        // Quadrilateral: split into two triangles
        const Vec2& p0 = nodes_[elem.node_ids[0]];
        const Vec2& p1 = nodes_[elem.node_ids[1]];
        const Vec2& p2 = nodes_[elem.node_ids[2]];
        const Vec2& p3 = nodes_[elem.node_ids[3]];

        Real area1 = 0.5 * std::abs((p1.x - p0.x) * (p2.y - p0.y) -
                                     (p2.x - p0.x) * (p1.y - p0.y));
        Real area2 = 0.5 * std::abs((p2.x - p0.x) * (p3.y - p0.y) -
                                     (p3.x - p0.x) * (p2.y - p0.y));
        return area1 + area2;
    }
}

Real Mesh::faceLength(Index face_id) const {
    const Face& face = faces_[face_id];
    const Element& elem = elements_[face.left_elem];
    int n_vertices = elem.numVertices();

    Index v0 = elem.node_ids[face.left_face];
    Index v1 = elem.node_ids[(face.left_face + 1) % n_vertices];

    return (nodes_[v1] - nodes_[v0]).norm();
}

Vec2 Mesh::faceNormal(Index face_id) const {
    const Face& face = faces_[face_id];
    const Element& elem = elements_[face.left_elem];
    int n_vertices = elem.numVertices();

    Index v0 = elem.node_ids[face.left_face];
    Index v1 = elem.node_ids[(face.left_face + 1) % n_vertices];

    Vec2 tangent = nodes_[v1] - nodes_[v0];
    Real len = tangent.norm();

    // Normal pointing outward from left element (90 degree rotation)
    Vec2 normal(tangent.y / len, -tangent.x / len);

    // Check orientation: normal should point from left element centroid outward
    Vec2 elem_center = elementCentroid(face.left_elem);
    Vec2 face_center = faceCentroid(face_id);
    Vec2 to_face = face_center - elem_center;

    if (normal.dot(to_face) < 0) {
        normal = normal * (-1.0);
    }

    return normal;
}

Vec2 Mesh::faceCentroid(Index face_id) const {
    const Face& face = faces_[face_id];
    const Element& elem = elements_[face.left_elem];
    int n_vertices = elem.numVertices();

    Index v0 = elem.node_ids[face.left_face];
    Index v1 = elem.node_ids[(face.left_face + 1) % n_vertices];

    return (nodes_[v0] + nodes_[v1]) * 0.5;
}

void Mesh::printStats() const {
    std::cout << "=== Mesh Statistics ===" << std::endl;
    std::cout << "Nodes: " << numNodes() << std::endl;
    std::cout << "Elements: " << numElements() << std::endl;
    std::cout << "Faces: " << numFaces() << std::endl;
    std::cout << "Boundary faces: " << numBoundaryFaces() << std::endl;

    // Count element types
    int n_tri = 0, n_quad = 0;
    int n_q1 = 0, n_q2 = 0, n_q4 = 0;
    for (const auto& elem : elements_) {
        if (elem.type == ElementType::Triangle) ++n_tri;
        else ++n_quad;
        if (elem.order == 1) ++n_q1;
        else if (elem.order == 2) ++n_q2;
        else if (elem.order == 4) ++n_q4;
    }

    std::cout << "Triangles: " << n_tri << std::endl;
    std::cout << "Quadrilaterals: " << n_quad << std::endl;
    if (n_q1) std::cout << "Q1 elements: " << n_q1 << std::endl;
    if (n_q2) std::cout << "Q2 elements: " << n_q2 << std::endl;
    if (n_q4) std::cout << "Q4 elements: " << n_q4 << std::endl;

    // Bounding box
    if (!nodes_.empty()) {
        Real xmin = nodes_[0].x, xmax = nodes_[0].x;
        Real ymin = nodes_[0].y, ymax = nodes_[0].y;
        for (const auto& node : nodes_) {
            xmin = std::min(xmin, node.x);
            xmax = std::max(xmax, node.x);
            ymin = std::min(ymin, node.y);
            ymax = std::max(ymax, node.y);
        }
        std::cout << "Bounding box: [" << xmin << ", " << xmax << "] x ["
                  << ymin << ", " << ymax << "]" << std::endl;
    }

    // Boundary condition summary
    std::cout << "Boundary conditions: " << bc_info_.size() << " defined" << std::endl;
    for (const auto& [tag, info] : bc_info_) {
        std::cout << "  Tag " << tag << ": " << info.name << " (type "
                  << static_cast<int>(info.type) << ")" << std::endl;
    }

    std::cout << "=======================" << std::endl;
}

}  // namespace zhijian
