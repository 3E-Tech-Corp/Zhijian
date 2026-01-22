#include "mesh/mesh.hpp"
#include <metis.h>
#include <iostream>
#include <algorithm>
#include <set>
#include <map>

namespace zhijian {

bool MeshPartitioner::partition(Mesh& mesh, int n_parts) {
    if (n_parts <= 1) {
        // No partitioning needed, assign all elements to partition 0
        for (auto& elem : mesh.elements()) {
            elem.partition = 0;
        }
        return true;
    }

    Index n_elem = mesh.numElements();
    if (n_elem < static_cast<Index>(n_parts)) {
        error_msg_ = "Number of elements less than number of partitions";
        return false;
    }

    // Build dual graph (element-to-element connectivity)
    // Two elements are connected if they share a face

    // Count neighbors per element
    std::vector<idx_t> xadj(n_elem + 1, 0);
    std::vector<std::vector<idx_t>> neighbors(n_elem);

    for (const auto& face : mesh.faces()) {
        if (!face.is_boundary) {
            neighbors[face.left_elem].push_back(face.right_elem);
            neighbors[face.right_elem].push_back(face.left_elem);
        }
    }

    // Build CSR format for METIS
    xadj[0] = 0;
    for (Index i = 0; i < n_elem; ++i) {
        xadj[i + 1] = xadj[i] + static_cast<idx_t>(neighbors[i].size());
    }

    std::vector<idx_t> adjncy;
    adjncy.reserve(xadj[n_elem]);
    for (Index i = 0; i < n_elem; ++i) {
        for (idx_t neighbor : neighbors[i]) {
            adjncy.push_back(neighbor);
        }
    }

    // METIS parameters
    idx_t nvtxs = static_cast<idx_t>(n_elem);
    idx_t ncon = 1;  // Number of balancing constraints
    idx_t nparts = static_cast<idx_t>(n_parts);
    idx_t objval;  // Edge-cut or communication volume
    std::vector<idx_t> part(n_elem);

    // Call METIS
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_NUMBERING] = 0;  // C-style 0-based numbering

    int ret = METIS_PartGraphKway(&nvtxs, &ncon, xadj.data(), adjncy.data(),
                                   nullptr, nullptr, nullptr,
                                   &nparts, nullptr, nullptr, options,
                                   &objval, part.data());

    if (ret != METIS_OK) {
        error_msg_ = "METIS partitioning failed";
        return false;
    }

    // Store partition info in elements
    for (Index i = 0; i < n_elem; ++i) {
        mesh.elements()[i].partition = part[i];
    }

    // Print partitioning statistics
    std::vector<int> part_count(n_parts, 0);
    for (Index i = 0; i < n_elem; ++i) {
        part_count[part[i]]++;
    }

    std::cout << "Mesh partitioned into " << n_parts << " parts" << std::endl;
    std::cout << "Edge cut: " << objval << std::endl;
    for (int p = 0; p < n_parts; ++p) {
        std::cout << "  Partition " << p << ": " << part_count[p] << " elements" << std::endl;
    }

    return true;
}

void MeshPartitioner::extractLocalMesh(const Mesh& global_mesh, int partition_id,
                                        Mesh& local_mesh,
                                        std::vector<Index>& local_to_global_elem,
                                        std::vector<Index>& local_to_global_node) {
    local_mesh = Mesh();
    local_to_global_elem.clear();
    local_to_global_node.clear();

    // Find elements belonging to this partition
    std::vector<Index> local_elems;
    for (Index i = 0; i < global_mesh.numElements(); ++i) {
        if (global_mesh.element(i).partition == partition_id) {
            local_elems.push_back(i);
        }
    }

    // Also include ghost elements (neighbors in other partitions)
    std::set<Index> ghost_elems;
    for (Index local_elem : local_elems) {
        for (const auto& face : global_mesh.faces()) {
            if (face.is_boundary) continue;

            if (face.left_elem == local_elem &&
                global_mesh.element(face.right_elem).partition != partition_id) {
                ghost_elems.insert(face.right_elem);
            }
            if (face.right_elem == local_elem &&
                global_mesh.element(face.left_elem).partition != partition_id) {
                ghost_elems.insert(face.left_elem);
            }
        }
    }

    // Build local-to-global element mapping
    local_to_global_elem = local_elems;
    for (Index ghost : ghost_elems) {
        local_to_global_elem.push_back(ghost);
    }

    // Build global-to-local element mapping
    std::map<Index, Index> global_to_local_elem;
    for (size_t i = 0; i < local_to_global_elem.size(); ++i) {
        global_to_local_elem[local_to_global_elem[i]] = static_cast<Index>(i);
    }

    // Find all nodes used by local and ghost elements
    std::set<Index> used_nodes;
    for (Index global_elem : local_to_global_elem) {
        for (Index node : global_mesh.element(global_elem).node_ids) {
            used_nodes.insert(node);
        }
    }

    // Build local-to-global node mapping
    local_to_global_node.assign(used_nodes.begin(), used_nodes.end());

    // Build global-to-local node mapping
    std::map<Index, Index> global_to_local_node;
    for (size_t i = 0; i < local_to_global_node.size(); ++i) {
        global_to_local_node[local_to_global_node[i]] = static_cast<Index>(i);
    }

    // Copy nodes
    local_mesh.nodes().resize(local_to_global_node.size());
    for (size_t i = 0; i < local_to_global_node.size(); ++i) {
        local_mesh.nodes()[i] = global_mesh.node(local_to_global_node[i]);
    }

    // Copy elements with remapped node indices
    local_mesh.elements().resize(local_to_global_elem.size());
    for (size_t i = 0; i < local_to_global_elem.size(); ++i) {
        const Element& global_elem = global_mesh.element(local_to_global_elem[i]);
        Element& local_elem = local_mesh.elements()[i];

        local_elem.type = global_elem.type;
        local_elem.order = global_elem.order;
        local_elem.partition = global_elem.partition;
        local_elem.node_ids.resize(global_elem.node_ids.size());

        for (size_t j = 0; j < global_elem.node_ids.size(); ++j) {
            local_elem.node_ids[j] = global_to_local_node.at(global_elem.node_ids[j]);
        }
    }

    // Build local faces
    local_mesh.buildFaces();

    // Update face connectivity for inter-partition faces
    for (auto& face : local_mesh.faces()) {
        if (!face.is_boundary) {
            Index global_left = local_to_global_elem[face.left_elem];
            Index global_right = local_to_global_elem[face.right_elem];

            int left_part = global_mesh.element(global_left).partition;
            int right_part = global_mesh.element(global_right).partition;

            if (left_part != partition_id) {
                face.mpi_rank = left_part;
            } else if (right_part != partition_id) {
                face.mpi_rank = right_part;
            }
        }
    }

    // Copy boundary condition info
    for (const auto& [tag, info] : global_mesh.allBCInfo()) {
        local_mesh.addBCInfo(tag, info);
    }

    // Compute geometry
    local_mesh.computeGeometry();
}

}  // namespace zhijian
