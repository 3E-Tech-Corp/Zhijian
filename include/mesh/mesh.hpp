#pragma once

#include "common/types.hpp"
#include <vector>
#include <map>
#include <set>
#include <string>

namespace zhijian {

// Face data structure
struct Face {
    LocalIndex left_elem;       // Left element index
    LocalIndex right_elem;      // Right element index (-1 for boundary)
    LocalIndex left_face;       // Local face index in left element
    LocalIndex right_face;      // Local face index in right element
    BCType bc_type;             // Boundary condition type
    int bc_tag;                 // Boundary condition tag/ID
    bool is_boundary;           // Is this a boundary face?
    int mpi_rank;               // MPI rank of neighbor (-1 if local or boundary)
};

// Element data structure
struct Element {
    ElementType type;           // Element type (tri/quad)
    int order;                  // Element order (1 for Q1, 2 for Q2)
    std::vector<Index> node_ids;  // Global node IDs
    std::vector<LocalIndex> face_ids;  // Face IDs
    int partition;              // Partition ID (MPI rank)

    // Number of vertices for this element type
    int numVertices() const {
        return (type == ElementType::Triangle) ? 3 : 4;
    }

    // Number of nodes (including high-order nodes)
    int numNodes() const {
        if (type == ElementType::Triangle) {
            return (order == 1) ? 3 : 6;  // Q1: 3 nodes, Q2: 6 nodes
        } else {
            return (order == 1) ? 4 : 9;  // Q1: 4 nodes, Q2: 9 nodes
        }
    }
};

// Boundary condition info
struct BCInfo {
    BCType type;
    std::string name;
    int tag;
    // BC-specific data
    Real p_total;               // Total pressure (inflow)
    Real T_total;               // Total temperature (inflow)
    Vec2 flow_direction;        // Flow direction
    Real p_static;              // Static pressure (outflow)
    State far_field_state;      // Far-field state
};

// Mesh class
class Mesh {
public:
    Mesh() = default;
    ~Mesh() = default;

    // Accessors
    Index numNodes() const { return nodes_.size(); }
    Index numElements() const { return elements_.size(); }
    Index numFaces() const { return faces_.size(); }
    Index numBoundaryFaces() const;

    const Vec2& node(Index i) const { return nodes_[i]; }
    Vec2& node(Index i) { return nodes_[i]; }

    const Element& element(Index i) const { return elements_[i]; }
    Element& element(Index i) { return elements_[i]; }

    const Face& face(Index i) const { return faces_[i]; }
    Face& face(Index i) { return faces_[i]; }

    const std::vector<Vec2>& nodes() const { return nodes_; }
    const std::vector<Element>& elements() const { return elements_; }
    const std::vector<Face>& faces() const { return faces_; }

    std::vector<Vec2>& nodes() { return nodes_; }
    std::vector<Element>& elements() { return elements_; }
    std::vector<Face>& faces() { return faces_; }

    // Boundary conditions
    void addBCInfo(int tag, const BCInfo& info) { bc_info_[tag] = info; }
    const BCInfo& getBCInfo(int tag) const { return bc_info_.at(tag); }
    bool hasBCInfo(int tag) const { return bc_info_.count(tag) > 0; }
    const std::map<int, BCInfo>& allBCInfo() const { return bc_info_; }

    // Build face connectivity from elements
    void buildFaces();

    // Compute geometric quantities
    void computeGeometry();

    // Get element centroid
    Vec2 elementCentroid(Index elem_id) const;

    // Get element area
    Real elementArea(Index elem_id) const;

    // Get face length
    Real faceLength(Index face_id) const;

    // Get face normal (pointing from left to right element)
    Vec2 faceNormal(Index face_id) const;

    // Get face centroid
    Vec2 faceCentroid(Index face_id) const;

    // Print mesh statistics
    void printStats() const;

private:
    std::vector<Vec2> nodes_;
    std::vector<Element> elements_;
    std::vector<Face> faces_;
    std::map<int, BCInfo> bc_info_;

    // Cached geometry data
    std::vector<Vec2> elem_centroids_;
    std::vector<Real> elem_areas_;
    std::vector<Real> face_lengths_;
    std::vector<Vec2> face_normals_;
    std::vector<Vec2> face_centroids_;
};

// CGNS mesh reader
class CGNSReader {
public:
    CGNSReader() = default;

    // Read mesh from CGNS file
    bool read(const std::string& filename, Mesh& mesh);

    // Get error message
    const std::string& errorMessage() const { return error_msg_; }

private:
    std::string error_msg_;

    bool readNodes(int file_id, int base_id, int zone_id, Mesh& mesh);
    bool readElements(int file_id, int base_id, int zone_id, Mesh& mesh);
    bool readBoundaryConditions(int file_id, int base_id, int zone_id, Mesh& mesh);
};

// Gmsh MSH mesh reader (v2.2 and v4.1 ASCII)
class GmshReader {
public:
    GmshReader() = default;

    // Read mesh from Gmsh .msh file
    // Populates nodes, elements, faces, and boundary conditions.
    // Calls mesh.buildFaces() and mesh.computeGeometry() internally.
    bool read(const std::string& filename, Mesh& mesh);

    // Apply stored boundary tags to faces (called internally after buildFaces)
    void applyBoundaryTags(Mesh& mesh);

    // Get error message
    const std::string& errorMessage() const { return error_msg_; }

private:
    std::string error_msg_;

    // Gmsh node-tag -> 0-based index mapping (may be non-contiguous in file)
    std::map<Index, Index> node_id_map_;

    // Boundary edge lookup: sorted vertex pair -> physical tag
    std::map<std::pair<Index, Index>, int> boundary_edge_tags_;

    // Temporary boundary edge storage used during parsing
    struct BoundaryEdge {
        std::vector<Index> nodes;   // 0-based node IDs (2 for Q1 line, 3 for Q2)
        int physical_tag;
    };

    // Version-specific section readers
    bool readNodesV2(std::ifstream& file, Mesh& mesh);
    bool readNodesV4(std::ifstream& file, Mesh& mesh);
    bool readElementsV2(std::ifstream& file, Mesh& mesh,
                        std::vector<BoundaryEdge>& boundary_edges);
    bool readElementsV4(std::ifstream& file, Mesh& mesh,
                        std::vector<BoundaryEdge>& boundary_edges,
                        const std::map<int, int>& entity_phys_map);

    // Build BCInfo entries and the boundary_edge_tags_ lookup
    void mapBoundaryConditions(Mesh& mesh,
                               const std::map<std::pair<int,int>, std::string>& physical_names,
                               const std::vector<BoundaryEdge>& boundary_edges);

    // Infer BCType from a physical group name (case-insensitive heuristics)
    static BCType inferBCType(const std::string& name);

    // Number of nodes for a given Gmsh element type (-1 if unsupported)
    static int gmshElemNumNodes(int elm_type);
};

// Mesh partitioner using METIS
class MeshPartitioner {
public:
    MeshPartitioner() = default;

    // Partition mesh into n_parts partitions
    bool partition(Mesh& mesh, int n_parts);

    // Extract local mesh for a given partition
    void extractLocalMesh(const Mesh& global_mesh, int partition_id,
                          Mesh& local_mesh,
                          std::vector<Index>& local_to_global_elem,
                          std::vector<Index>& local_to_global_node);

    // Get error message
    const std::string& errorMessage() const { return error_msg_; }

private:
    std::string error_msg_;
};

}  // namespace zhijian
