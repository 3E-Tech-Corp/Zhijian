#include "mesh/mesh.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <set>
#include <cctype>

namespace zhijian {

// ============================================================================
// Static helpers
// ============================================================================

int GmshReader::gmshElemNumNodes(int elm_type) {
    switch (elm_type) {
        case 1:  return 2;   // 2-node line
        case 2:  return 3;   // 3-node triangle
        case 3:  return 4;   // 4-node quadrilateral
        case 4:  return 4;   // 4-node tetrahedron (3D — detected but rejected)
        case 5:  return 8;   // 8-node hexahedron (3D — detected but rejected)
        case 6:  return 6;   // 6-node prism (3D — detected but rejected)
        case 7:  return 5;   // 5-node pyramid (3D — detected but rejected)
        case 8:  return 3;   // 3-node second-order line
        case 9:  return 6;   // 6-node second-order triangle
        case 10: return 9;   // 9-node second-order quadrilateral
        case 11: return 10;  // 10-node second-order tetrahedron (3D)
        case 12: return 27;  // 27-node second-order hexahedron (3D)
        case 13: return 18;  // 18-node second-order prism (3D)
        case 14: return 14;  // 14-node second-order pyramid (3D)
        case 15: return 1;   // 1-node point
        case 16: return 8;   // 8-node second-order quadrilateral (serendipity)
        case 20: return 9;   // 9-node third-order triangle
        case 21: return 10;  // 10-node third-order triangle
        case 36: return 16;  // 16-node third-order quadrilateral
        case 37: return 25;  // 25-node fourth-order quadrilateral
        default: return -1;  // unsupported
    }
}

const char* GmshReader::gmshElemTypeName(int elm_type) {
    switch (elm_type) {
        case 1:  return "2-node line";
        case 2:  return "3-node triangle";
        case 3:  return "4-node quadrilateral";
        case 4:  return "4-node tetrahedron";
        case 5:  return "8-node hexahedron";
        case 6:  return "6-node prism";
        case 7:  return "5-node pyramid";
        case 8:  return "3-node line (order 2)";
        case 9:  return "6-node triangle (order 2)";
        case 10: return "9-node quadrilateral (order 2)";
        case 11: return "10-node tetrahedron (order 2)";
        case 12: return "27-node hexahedron (order 2)";
        case 13: return "18-node prism (order 2)";
        case 14: return "14-node pyramid (order 2)";
        case 15: return "1-node point";
        case 16: return "8-node quadrilateral (serendipity)";
        case 20: return "9-node triangle (order 3)";
        case 21: return "10-node triangle (order 3)";
        case 36: return "16-node quadrilateral (order 3)";
        case 37: return "25-node quadrilateral (order 4)";
        default: return "unknown";
    }
}

bool GmshReader::is3DElement(int elm_type) {
    return elm_type == 4 || elm_type == 5 || elm_type == 6 || elm_type == 7 ||
           elm_type == 11 || elm_type == 12 || elm_type == 13 || elm_type == 14;
}

bool GmshReader::is2DElement(int elm_type) {
    return elm_type == 2 || elm_type == 3 || elm_type == 9 || elm_type == 10 ||
           elm_type == 16 || elm_type == 20 || elm_type == 21 || elm_type == 36 ||
           elm_type == 37;
}

Index GmshReader::addSerendipityCenterNode(Mesh& mesh, const std::vector<Index>& node_ids) {
    // Synthesize the center node for an 8-node serendipity quad by converting
    // to a 9-node biquadratic quad.  For a serendipity quad the Gmsh node
    // ordering is:
    //
    //   3 -- 6 -- 2
    //   |         |
    //   7    ?    5
    //   |         |
    //   0 -- 4 -- 1
    //
    // The proper conversion uses:
    //   center = 0.5*(n4+n5+n6+n7) - 0.25*(n0+n1+n2+n3)
    // which is exact for bilinear geometry.

    const Vec2& n0 = mesh.node(node_ids[0]);
    const Vec2& n1 = mesh.node(node_ids[1]);
    const Vec2& n2 = mesh.node(node_ids[2]);
    const Vec2& n3 = mesh.node(node_ids[3]);
    const Vec2& n4 = mesh.node(node_ids[4]);
    const Vec2& n5 = mesh.node(node_ids[5]);
    const Vec2& n6 = mesh.node(node_ids[6]);
    const Vec2& n7 = mesh.node(node_ids[7]);

    Vec2 center;
    center.x = 0.5 * (n4.x + n5.x + n6.x + n7.x) - 0.25 * (n0.x + n1.x + n2.x + n3.x);
    center.y = 0.5 * (n4.y + n5.y + n6.y + n7.y) - 0.25 * (n0.y + n1.y + n2.y + n3.y);

    Index new_idx = static_cast<Index>(mesh.nodes().size());
    mesh.nodes().push_back(center);
    return new_idx;
}

BCType GmshReader::inferBCType(const std::string& name) {
    std::string lower = name;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    if (lower.find("wall") != std::string::npos) {
        if (lower.find("slip") != std::string::npos ||
            lower.find("inviscid") != std::string::npos)
            return BCType::SlipWall;
        return BCType::Wall;
    }
    if (lower.find("airfoil") != std::string::npos ||
        lower.find("wing") != std::string::npos ||
        lower.find("body") != std::string::npos ||
        lower.find("surface") != std::string::npos)
        return BCType::Wall;
    if (lower.find("farfield") != std::string::npos ||
        lower.find("far_field") != std::string::npos ||
        lower.find("far-field") != std::string::npos ||
        lower.find("freestream") != std::string::npos)
        return BCType::FarField;
    if (lower.find("inflow") != std::string::npos ||
        lower.find("inlet") != std::string::npos)
        return BCType::Inflow;
    if (lower.find("outflow") != std::string::npos ||
        lower.find("outlet") != std::string::npos ||
        lower.find("exit") != std::string::npos ||
        lower == "out")
        return BCType::Outflow;
    if (lower.find("symmetry") != std::string::npos ||
        lower.find("sym") != std::string::npos)
        return BCType::Symmetry;
    if (lower.find("periodic") != std::string::npos)
        return BCType::Periodic;

    std::cerr << "Warning: Could not infer BC type for '" << name
              << "', defaulting to Wall\n";
    return BCType::Wall;
}

// ============================================================================
// Top-level read — detect format version, then dispatch
// ============================================================================

bool GmshReader::read(const std::string& filename, Mesh& mesh) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        error_msg_ = "Failed to open Gmsh file: " + filename;
        return false;
    }

    // ---- Detect format version from $MeshFormat ----
    double version = 0.0;
    int file_type = 0, data_size = 0;
    std::string line;

    while (std::getline(file, line)) {
        if (line.find("$MeshFormat") != std::string::npos) {
            std::getline(file, line);
            std::istringstream iss(line);
            iss >> version >> file_type >> data_size;
            // consume $EndMeshFormat
            std::getline(file, line);
            break;
        }
    }

    if (version == 0.0) {
        error_msg_ = "No $MeshFormat section found in: " + filename;
        return false;
    }
    if (file_type != 0) {
        error_msg_ = "Only ASCII Gmsh files are supported (file_type must be 0)";
        return false;
    }

    // ---- Rewind and parse all sections ----
    file.clear();
    file.seekg(0);

    // Physical names: (dimension, tag) -> name
    std::map<std::pair<int, int>, std::string> physical_names;

    // v4 entity-tag -> physical-tag mapping (curves only, for boundary edges)
    std::map<int, int> entity_phys_map;

    // Collected boundary edges
    std::vector<BoundaryEdge> boundary_edges;

    // Track element types encountered for diagnostic reporting
    std::map<int, int> element_type_counts;

    // Clear any state from a previous read
    node_id_map_.clear();
    boundary_edge_tags_.clear();

    bool nodes_read = false;
    bool elements_read = false;

    while (std::getline(file, line)) {
        // ---- $MeshFormat (skip, already parsed) ----
        if (line.find("$MeshFormat") != std::string::npos) {
            while (std::getline(file, line)) {
                if (line.find("$EndMeshFormat") != std::string::npos) break;
            }
        }

        // ---- $PhysicalNames ----
        else if (line.find("$PhysicalNames") != std::string::npos) {
            int num_names;
            std::getline(file, line);
            std::istringstream hdr(line);
            hdr >> num_names;

            for (int i = 0; i < num_names; ++i) {
                std::getline(file, line);
                std::istringstream lss(line);
                int dim, tag;
                std::string name;
                lss >> dim >> tag;
                // Read the rest of the line which contains the quoted name
                std::getline(lss, name);
                // Trim leading whitespace
                auto start = name.find_first_not_of(" \t");
                if (start != std::string::npos) name = name.substr(start);
                // Strip quotes
                if (!name.empty() && name.front() == '"') name = name.substr(1);
                if (!name.empty() && name.back() == '"') name.pop_back();

                physical_names[{dim, tag}] = name;
            }
            // consume $EndPhysicalNames
            while (std::getline(file, line)) {
                if (line.find("$EndPhysicalNames") != std::string::npos) break;
            }
        }

        // ---- $Entities (v4 only) ----
        else if (line.find("$Entities") != std::string::npos) {
            std::getline(file, line);
            std::istringstream hdr(line);
            int num_points, num_curves, num_surfaces, num_volumes;
            hdr >> num_points >> num_curves >> num_surfaces >> num_volumes;

            // Skip point entities
            for (int i = 0; i < num_points; ++i) {
                std::getline(file, line);
            }

            // Parse curve entities: extract entity_tag -> first physical tag
            for (int i = 0; i < num_curves; ++i) {
                std::getline(file, line);
                std::istringstream iss(line);
                int entity_tag;
                double xmin, ymin, zmin, xmax, ymax, zmax;
                int num_phys;
                iss >> entity_tag >> xmin >> ymin >> zmin >> xmax >> ymax >> zmax >> num_phys;
                if (num_phys > 0) {
                    int phys_tag;
                    iss >> phys_tag;
                    entity_phys_map[entity_tag] = phys_tag;
                }
            }

            // Skip surface and volume entities
            for (int i = 0; i < num_surfaces + num_volumes; ++i) {
                std::getline(file, line);
            }

            while (std::getline(file, line)) {
                if (line.find("$EndEntities") != std::string::npos) break;
            }
        }

        // ---- $Nodes ----
        else if (line.find("$Nodes") != std::string::npos) {
            if (version >= 4.0) {
                nodes_read = readNodesV4(file, mesh);
            } else {
                nodes_read = readNodesV2(file, mesh);
            }
            if (!nodes_read) return false;
        }

        // ---- $Elements ----
        else if (line.find("$Elements") != std::string::npos) {
            if (version >= 4.0) {
                elements_read = readElementsV4(file, mesh, boundary_edges, entity_phys_map, element_type_counts);
            } else {
                elements_read = readElementsV2(file, mesh, boundary_edges, element_type_counts);
            }
            if (!elements_read) return false;
        }

        // All other sections — skip to matching $End tag
        else if (!line.empty() && line[0] == '$' &&
                 line.find("$End") == std::string::npos) {
            std::string section_name = line.substr(1);
            std::string end_marker = "$End" + section_name;
            while (std::getline(file, line)) {
                if (line.find(end_marker) != std::string::npos) break;
            }
        }
    }

    if (!nodes_read) {
        error_msg_ = "No $Nodes section found";
        return false;
    }
    if (!elements_read) {
        error_msg_ = "No $Elements section found";
        return false;
    }
    if (mesh.numElements() == 0) {
        std::string diag = "No 2D elements (triangles/quads) found in Gmsh file.";
        if (element_type_counts.empty()) {
            diag += " The $Elements section contained no elements at all.";
        } else {
            diag += "\nElement types found in file:";
            bool has_3d = false;
            for (const auto& kv : element_type_counts) {
                diag += "\n  Type " + std::to_string(kv.first) + " ("
                     + gmshElemTypeName(kv.first) + "): " + std::to_string(kv.second);
                if (is3DElement(kv.first)) has_3d = true;
            }
            if (has_3d) {
                diag += "\n\nThis appears to be a 3D mesh. Zhijian is a 2D solver and "
                        "requires triangles (type 2/9) or quadrilaterals (type 3/10/16/37). "
                        "Please regenerate the mesh as 2D, or export only the surface elements.";
            } else {
                diag += "\n\nNone of the element types above are supported 2D elements. "
                        "Zhijian requires triangles (type 2/9) or quadrilaterals (type 3/10/16/37).";
            }
        }
        error_msg_ = diag;
        return false;
    }

    // ---- Map physical names to boundary conditions ----
    mapBoundaryConditions(mesh, physical_names, boundary_edges);

    // ---- Build face connectivity, then stamp boundary tags ----
    mesh.buildFaces();
    applyBoundaryTags(mesh);
    mesh.computeGeometry();

    std::cout << "Gmsh mesh loaded: "
              << mesh.numNodes() << " nodes, "
              << mesh.numElements() << " elements, "
              << boundary_edges.size() << " boundary edges, "
              << mesh.allBCInfo().size() << " boundary groups\n";

    return true;
}

// ============================================================================
// v2.2 Nodes: node-number x y z
// ============================================================================

bool GmshReader::readNodesV2(std::ifstream& file, Mesh& mesh) {
    std::string line;
    std::getline(file, line);
    int num_nodes = std::stoi(line);

    mesh.nodes().resize(num_nodes);

    for (int i = 0; i < num_nodes; ++i) {
        std::getline(file, line);
        std::istringstream iss(line);
        int gmsh_id;
        double x, y, z;
        iss >> gmsh_id >> x >> y >> z;
        node_id_map_[gmsh_id] = i;
        mesh.nodes()[i] = Vec2(x, y);
    }

    // consume $EndNodes
    std::getline(file, line);
    return true;
}

// ============================================================================
// v4.1 Nodes: entity blocks with separate tag listing and coordinate listing
// ============================================================================

bool GmshReader::readNodesV4(std::ifstream& file, Mesh& mesh) {
    std::string line;
    std::getline(file, line);
    std::istringstream hdr(line);
    int num_entity_blocks, num_nodes, min_tag, max_tag;
    hdr >> num_entity_blocks >> num_nodes >> min_tag >> max_tag;

    mesh.nodes().resize(num_nodes);

    int global_idx = 0;
    for (int b = 0; b < num_entity_blocks; ++b) {
        std::getline(file, line);
        std::istringstream bss(line);
        int entity_dim, entity_tag, parametric, nodes_in_block;
        bss >> entity_dim >> entity_tag >> parametric >> nodes_in_block;

        // Read node tags first (one per line)
        std::vector<int> block_tags(nodes_in_block);
        for (int i = 0; i < nodes_in_block; ++i) {
            std::getline(file, line);
            block_tags[i] = std::stoi(line);
            node_id_map_[block_tags[i]] = global_idx + i;
        }

        // Read coordinates (one per line: x y z)
        for (int i = 0; i < nodes_in_block; ++i) {
            std::getline(file, line);
            std::istringstream css(line);
            double x, y, z;
            css >> x >> y >> z;
            mesh.nodes()[global_idx + i] = Vec2(x, y);
        }

        global_idx += nodes_in_block;
    }

    // consume $EndNodes
    std::getline(file, line);
    return true;
}

// ============================================================================
// v2.2 Elements: elm-number elm-type num-tags tag... node...
// ============================================================================

bool GmshReader::readElementsV2(std::ifstream& file, Mesh& mesh,
                                 std::vector<BoundaryEdge>& boundary_edges,
                                 std::map<int, int>& element_type_counts) {
    std::string line;
    std::getline(file, line);
    int num_elements = std::stoi(line);

    for (int i = 0; i < num_elements; ++i) {
        std::getline(file, line);
        std::istringstream iss(line);

        int elm_number, elm_type, num_tags;
        iss >> elm_number >> elm_type >> num_tags;
        element_type_counts[elm_type]++;

        // First tag is the physical group; second is the elementary entity
        int physical_tag = 0;
        for (int t = 0; t < num_tags; ++t) {
            int tag;
            iss >> tag;
            if (t == 0) physical_tag = tag;
        }

        int nn = gmshElemNumNodes(elm_type);
        if (nn < 0) continue;  // skip unsupported types

        // Read and remap node IDs
        std::vector<Index> node_ids(nn);
        for (int j = 0; j < nn; ++j) {
            int gmsh_nid;
            iss >> gmsh_nid;
            auto it = node_id_map_.find(gmsh_nid);
            if (it == node_id_map_.end()) {
                error_msg_ = "Element references unknown node tag " + std::to_string(gmsh_nid);
                return false;
            }
            node_ids[j] = it->second;
        }

        // Classify by element type
        if (elm_type == 1 || elm_type == 8) {
            // 1D boundary line
            BoundaryEdge edge;
            edge.nodes = node_ids;
            edge.physical_tag = physical_tag;
            boundary_edges.push_back(std::move(edge));
        }
        else if (elm_type == 2) {
            Element elem;
            elem.type = ElementType::Triangle;
            elem.order = 1;
            elem.node_ids = node_ids;
            elem.partition = 0;
            mesh.elements().push_back(std::move(elem));
        }
        else if (elm_type == 9) {
            Element elem;
            elem.type = ElementType::Triangle;
            elem.order = 2;
            elem.node_ids = node_ids;
            elem.partition = 0;
            mesh.elements().push_back(std::move(elem));
        }
        else if (elm_type == 3) {
            Element elem;
            elem.type = ElementType::Quadrilateral;
            elem.order = 1;
            elem.node_ids = node_ids;
            elem.partition = 0;
            mesh.elements().push_back(std::move(elem));
        }
        else if (elm_type == 10) {
            Element elem;
            elem.type = ElementType::Quadrilateral;
            elem.order = 2;
            elem.node_ids = node_ids;
            elem.partition = 0;
            mesh.elements().push_back(std::move(elem));
        }
        else if (elm_type == 16) {
            // 8-node serendipity quad → convert to 9-node biquadratic quad
            // by synthesizing the center node
            Index center = addSerendipityCenterNode(mesh, node_ids);
            node_ids.push_back(center);
            Element elem;
            elem.type = ElementType::Quadrilateral;
            elem.order = 2;
            elem.node_ids = node_ids;
            elem.partition = 0;
            mesh.elements().push_back(std::move(elem));
        }
        else if (elm_type == 37) {
            // 25-node fourth-order (quartic) quadrilateral
            Element elem;
            elem.type = ElementType::Quadrilateral;
            elem.order = 4;
            elem.node_ids = node_ids;
            elem.partition = 0;
            mesh.elements().push_back(std::move(elem));
        }
        // else: skip (e.g. type 15 = point, 3D elements)
    }

    // consume $EndElements
    std::getline(file, line);
    return true;
}

// ============================================================================
// v4.1 Elements: entity blocks with (entityDim entityTag elemType count)
// ============================================================================

bool GmshReader::readElementsV4(std::ifstream& file, Mesh& mesh,
                                 std::vector<BoundaryEdge>& boundary_edges,
                                 const std::map<int, int>& entity_phys_map,
                                 std::map<int, int>& element_type_counts) {
    std::string line;
    std::getline(file, line);
    std::istringstream hdr(line);
    int num_entity_blocks, num_elements, min_tag, max_tag;
    hdr >> num_entity_blocks >> num_elements >> min_tag >> max_tag;

    for (int b = 0; b < num_entity_blocks; ++b) {
        std::getline(file, line);
        std::istringstream bss(line);
        int entity_dim, entity_tag, elm_type, elems_in_block;
        bss >> entity_dim >> entity_tag >> elm_type >> elems_in_block;
        element_type_counts[elm_type] += elems_in_block;

        int nn = gmshElemNumNodes(elm_type);

        // Resolve physical tag for this entity via the entity->physical map
        int physical_tag = 0;
        if (entity_dim == 1) {
            auto it = entity_phys_map.find(entity_tag);
            if (it != entity_phys_map.end()) physical_tag = it->second;
        }

        for (int e = 0; e < elems_in_block; ++e) {
            std::getline(file, line);

            if (nn < 0 || elm_type == 15) continue;  // skip unsupported/point

            std::istringstream ess(line);
            int elm_tag;
            ess >> elm_tag;

            std::vector<Index> node_ids(nn);
            for (int j = 0; j < nn; ++j) {
                int gmsh_nid;
                ess >> gmsh_nid;
                auto it = node_id_map_.find(gmsh_nid);
                if (it == node_id_map_.end()) {
                    error_msg_ = "Element references unknown node tag " + std::to_string(gmsh_nid);
                    return false;
                }
                node_ids[j] = it->second;
            }

            if (elm_type == 1 || elm_type == 8) {
                BoundaryEdge edge;
                edge.nodes = node_ids;
                edge.physical_tag = physical_tag;
                boundary_edges.push_back(std::move(edge));
            }
            else if (elm_type == 2) {
                Element elem;
                elem.type = ElementType::Triangle;
                elem.order = 1;
                elem.node_ids = node_ids;
                elem.partition = 0;
                mesh.elements().push_back(std::move(elem));
            }
            else if (elm_type == 9) {
                Element elem;
                elem.type = ElementType::Triangle;
                elem.order = 2;
                elem.node_ids = node_ids;
                elem.partition = 0;
                mesh.elements().push_back(std::move(elem));
            }
            else if (elm_type == 3) {
                Element elem;
                elem.type = ElementType::Quadrilateral;
                elem.order = 1;
                elem.node_ids = node_ids;
                elem.partition = 0;
                mesh.elements().push_back(std::move(elem));
            }
            else if (elm_type == 10) {
                Element elem;
                elem.type = ElementType::Quadrilateral;
                elem.order = 2;
                elem.node_ids = node_ids;
                elem.partition = 0;
                mesh.elements().push_back(std::move(elem));
            }
            else if (elm_type == 16) {
                // 8-node serendipity quad → convert to 9-node biquadratic quad
                // by synthesizing the center node
                Index center = addSerendipityCenterNode(mesh, node_ids);
                node_ids.push_back(center);
                Element elem;
                elem.type = ElementType::Quadrilateral;
                elem.order = 2;
                elem.node_ids = node_ids;
                elem.partition = 0;
                mesh.elements().push_back(std::move(elem));
            }
            else if (elm_type == 37) {
                // 25-node fourth-order (quartic) quadrilateral
                Element elem;
                elem.type = ElementType::Quadrilateral;
                elem.order = 4;
                elem.node_ids = node_ids;
                elem.partition = 0;
                mesh.elements().push_back(std::move(elem));
            }
        }
    }

    // consume $EndElements
    std::getline(file, line);
    return true;
}

// ============================================================================
// Map physical group names to BCInfo entries and build edge lookup
// ============================================================================

void GmshReader::mapBoundaryConditions(
        Mesh& mesh,
        const std::map<std::pair<int, int>, std::string>& physical_names,
        const std::vector<BoundaryEdge>& boundary_edges) {

    // Collect unique physical tags from boundary edges
    std::set<int> bc_tags;
    for (const auto& edge : boundary_edges) {
        if (edge.physical_tag > 0) {
            bc_tags.insert(edge.physical_tag);
        }
    }

    // Register BCInfo for each tag
    for (int tag : bc_tags) {
        BCInfo info;
        info.tag = tag;

        // Look up name from $PhysicalNames (dim=1 for boundary lines)
        auto it = physical_names.find({1, tag});
        if (it != physical_names.end()) {
            info.name = it->second;
        } else {
            info.name = "BC_" + std::to_string(tag);
        }

        info.type = inferBCType(info.name);

        // Zero-init BC-specific data
        info.p_total = 0;
        info.T_total = 0;
        info.flow_direction = Vec2(1, 0);
        info.p_static = 0;
        info.far_field_state = State();

        mesh.addBCInfo(tag, info);
    }

    // Build sorted vertex-pair -> physical_tag lookup for applyBoundaryTags()
    boundary_edge_tags_.clear();
    for (const auto& edge : boundary_edges) {
        if (edge.nodes.size() >= 2 && edge.physical_tag > 0) {
            Index n0 = edge.nodes[0];
            Index n1 = edge.nodes[1];
            if (n0 > n1) std::swap(n0, n1);
            boundary_edge_tags_[{n0, n1}] = edge.physical_tag;
        }
    }
}

// ============================================================================
// Stamp bc_tag / bc_type on boundary faces after buildFaces()
// ============================================================================

void GmshReader::applyBoundaryTags(Mesh& mesh) {
    for (Index i = 0; i < mesh.numFaces(); ++i) {
        Face& face = mesh.face(i);
        if (!face.is_boundary) continue;

        // Extract the two vertex node IDs for this face from its parent element
        const Element& elem = mesh.element(face.left_elem);
        int nv = elem.numVertices();
        int lf = face.left_face;

        Index v0 = elem.node_ids[lf % nv];
        Index v1 = elem.node_ids[(lf + 1) % nv];
        if (v0 > v1) std::swap(v0, v1);

        auto it = boundary_edge_tags_.find({v0, v1});
        if (it != boundary_edge_tags_.end()) {
            face.bc_tag = it->second;
            if (mesh.hasBCInfo(it->second)) {
                face.bc_type = mesh.getBCInfo(it->second).type;
            }
        }
    }
}

}  // namespace zhijian
