#include "mesh/mesh.hpp"
#include <cgnslib.h>
#include <iostream>
#include <cstring>

namespace zhijian {

bool CGNSReader::read(const std::string& filename, Mesh& mesh) {
    int file_id;
    int ierr = cg_open(filename.c_str(), CG_MODE_READ, &file_id);
    if (ierr != CG_OK) {
        error_msg_ = "Failed to open CGNS file: " + filename + " - " + cg_get_error();
        return false;
    }

    // Get number of bases
    int n_bases;
    ierr = cg_nbases(file_id, &n_bases);
    if (ierr != CG_OK || n_bases < 1) {
        error_msg_ = "No bases found in CGNS file";
        cg_close(file_id);
        return false;
    }

    // Read from first base
    int base_id = 1;
    char base_name[33];
    int cell_dim, phys_dim;
    ierr = cg_base_read(file_id, base_id, base_name, &cell_dim, &phys_dim);
    if (ierr != CG_OK) {
        error_msg_ = "Failed to read base: " + std::string(cg_get_error());
        cg_close(file_id);
        return false;
    }

    if (cell_dim != 2) {
        error_msg_ = "Expected 2D mesh (cell_dim=2), got cell_dim=" + std::to_string(cell_dim);
        cg_close(file_id);
        return false;
    }

    // Get number of zones
    int n_zones;
    ierr = cg_nzones(file_id, base_id, &n_zones);
    if (ierr != CG_OK || n_zones < 1) {
        error_msg_ = "No zones found in CGNS file";
        cg_close(file_id);
        return false;
    }

    // Read from first zone (assuming single-zone mesh for now)
    int zone_id = 1;

    // Read nodes
    if (!readNodes(file_id, base_id, zone_id, mesh)) {
        cg_close(file_id);
        return false;
    }

    // Read elements
    if (!readElements(file_id, base_id, zone_id, mesh)) {
        cg_close(file_id);
        return false;
    }

    // Read boundary conditions
    if (!readBoundaryConditions(file_id, base_id, zone_id, mesh)) {
        // Warning only, not fatal
        std::cerr << "Warning: " << error_msg_ << std::endl;
    }

    cg_close(file_id);

    // Build face connectivity
    mesh.buildFaces();

    // Compute geometry
    mesh.computeGeometry();

    return true;
}

bool CGNSReader::readNodes(int file_id, int base_id, int zone_id, Mesh& mesh) {
    char zone_name[33];
    cgsize_t sizes[3];
    ZoneType_t zone_type;

    int ierr = cg_zone_read(file_id, base_id, zone_id, zone_name, sizes);
    if (ierr != CG_OK) {
        error_msg_ = "Failed to read zone: " + std::string(cg_get_error());
        return false;
    }

    ierr = cg_zone_type(file_id, base_id, zone_id, &zone_type);
    if (ierr != CG_OK || zone_type != Unstructured) {
        error_msg_ = "Only unstructured zones supported";
        return false;
    }

    cgsize_t n_nodes = sizes[0];
    cgsize_t n_cells = sizes[1];

    // Read coordinates
    std::vector<double> x(n_nodes), y(n_nodes);

    cgsize_t range_min = 1, range_max = n_nodes;

    ierr = cg_coord_read(file_id, base_id, zone_id, "CoordinateX",
                         RealDouble, &range_min, &range_max, x.data());
    if (ierr != CG_OK) {
        error_msg_ = "Failed to read X coordinates: " + std::string(cg_get_error());
        return false;
    }

    ierr = cg_coord_read(file_id, base_id, zone_id, "CoordinateY",
                         RealDouble, &range_min, &range_max, y.data());
    if (ierr != CG_OK) {
        error_msg_ = "Failed to read Y coordinates: " + std::string(cg_get_error());
        return false;
    }

    // Store nodes in mesh
    mesh.nodes().resize(n_nodes);
    for (cgsize_t i = 0; i < n_nodes; ++i) {
        mesh.nodes()[i] = Vec2(x[i], y[i]);
    }

    return true;
}

bool CGNSReader::readElements(int file_id, int base_id, int zone_id, Mesh& mesh) {
    // Get number of element sections
    int n_sections;
    int ierr = cg_nsections(file_id, base_id, zone_id, &n_sections);
    if (ierr != CG_OK) {
        error_msg_ = "Failed to get number of sections: " + std::string(cg_get_error());
        return false;
    }

    mesh.elements().clear();

    for (int sect_id = 1; sect_id <= n_sections; ++sect_id) {
        char section_name[33];
        ElementType_t cgns_elem_type;
        cgsize_t elem_start, elem_end;
        int n_bndry, parent_flag;

        ierr = cg_section_read(file_id, base_id, zone_id, sect_id,
                               section_name, &cgns_elem_type,
                               &elem_start, &elem_end, &n_bndry, &parent_flag);
        if (ierr != CG_OK) {
            error_msg_ = "Failed to read section: " + std::string(cg_get_error());
            return false;
        }

        // Skip boundary element sections (1D elements)
        if (cgns_elem_type == BAR_2 || cgns_elem_type == BAR_3) {
            continue;
        }

        // Determine element type and order
        ElementType elem_type;
        int elem_order;
        int nodes_per_elem;

        switch (cgns_elem_type) {
            case TRI_3:
                elem_type = ElementType::Triangle;
                elem_order = 1;
                nodes_per_elem = 3;
                break;
            case TRI_6:
                elem_type = ElementType::Triangle;
                elem_order = 2;
                nodes_per_elem = 6;
                break;
            case QUAD_4:
                elem_type = ElementType::Quadrilateral;
                elem_order = 1;
                nodes_per_elem = 4;
                break;
            case QUAD_9:
                elem_type = ElementType::Quadrilateral;
                elem_order = 2;
                nodes_per_elem = 9;
                break;
            case MIXED:
                // Handle mixed elements
                break;
            default:
                std::cerr << "Skipping unsupported element type: " << cgns_elem_type << std::endl;
                continue;
        }

        cgsize_t n_elems = elem_end - elem_start + 1;

        if (cgns_elem_type == MIXED) {
            // Read mixed element connectivity
            cgsize_t conn_size;
            ierr = cg_ElementDataSize(file_id, base_id, zone_id, sect_id, &conn_size);
            if (ierr != CG_OK) {
                error_msg_ = "Failed to get element data size: " + std::string(cg_get_error());
                return false;
            }

            std::vector<cgsize_t> conn(conn_size);
            ierr = cg_elements_read(file_id, base_id, zone_id, sect_id,
                                    conn.data(), nullptr);
            if (ierr != CG_OK) {
                error_msg_ = "Failed to read elements: " + std::string(cg_get_error());
                return false;
            }

            // Parse mixed connectivity
            cgsize_t idx = 0;
            while (idx < conn_size) {
                ElementType_t mixed_type = static_cast<ElementType_t>(conn[idx++]);

                switch (mixed_type) {
                    case TRI_3:
                        elem_type = ElementType::Triangle;
                        elem_order = 1;
                        nodes_per_elem = 3;
                        break;
                    case TRI_6:
                        elem_type = ElementType::Triangle;
                        elem_order = 2;
                        nodes_per_elem = 6;
                        break;
                    case QUAD_4:
                        elem_type = ElementType::Quadrilateral;
                        elem_order = 1;
                        nodes_per_elem = 4;
                        break;
                    case QUAD_9:
                        elem_type = ElementType::Quadrilateral;
                        elem_order = 2;
                        nodes_per_elem = 9;
                        break;
                    default:
                        std::cerr << "Skipping unsupported mixed element type: "
                                  << mixed_type << std::endl;
                        // Skip unknown element (assume 3 nodes minimum)
                        idx += 3;
                        continue;
                }

                Element elem;
                elem.type = elem_type;
                elem.order = elem_order;
                elem.node_ids.resize(nodes_per_elem);
                for (int i = 0; i < nodes_per_elem; ++i) {
                    elem.node_ids[i] = conn[idx++] - 1;  // CGNS uses 1-based indexing
                }
                mesh.elements().push_back(elem);
            }
        } else {
            // Read connectivity for homogeneous element section
            std::vector<cgsize_t> conn(n_elems * nodes_per_elem);
            ierr = cg_elements_read(file_id, base_id, zone_id, sect_id,
                                    conn.data(), nullptr);
            if (ierr != CG_OK) {
                error_msg_ = "Failed to read elements: " + std::string(cg_get_error());
                return false;
            }

            // Create elements
            for (cgsize_t i = 0; i < n_elems; ++i) {
                Element elem;
                elem.type = elem_type;
                elem.order = elem_order;
                elem.node_ids.resize(nodes_per_elem);
                for (int j = 0; j < nodes_per_elem; ++j) {
                    elem.node_ids[j] = conn[i * nodes_per_elem + j] - 1;
                }
                mesh.elements().push_back(elem);
            }
        }
    }

    return true;
}

bool CGNSReader::readBoundaryConditions(int file_id, int base_id, int zone_id, Mesh& mesh) {
    int n_bcs;
    int ierr = cg_nbocos(file_id, base_id, zone_id, &n_bcs);
    if (ierr != CG_OK) {
        error_msg_ = "Failed to get number of BCs: " + std::string(cg_get_error());
        return false;
    }

    for (int bc_id = 1; bc_id <= n_bcs; ++bc_id) {
        char bc_name[33];
        BCType_t cgns_bc_type;
        PointSetType_t point_set_type;
        cgsize_t n_pts;
        int normal_index;
        cgsize_t normal_list_size;
        DataType_t normal_data_type;
        int n_datasets;

        ierr = cg_boco_info(file_id, base_id, zone_id, bc_id,
                            bc_name, &cgns_bc_type, &point_set_type,
                            &n_pts, nullptr, &normal_list_size,
                            &normal_data_type, &n_datasets);
        if (ierr != CG_OK) {
            error_msg_ = "Failed to read BC info: " + std::string(cg_get_error());
            continue;
        }

        // Map CGNS BC type to our BC type
        BCInfo bc_info;
        bc_info.name = bc_name;
        bc_info.tag = bc_id;

        switch (cgns_bc_type) {
            case BCWall:
            case BCWallViscous:
            case BCWallViscousHeatFlux:
            case BCWallViscousIsothermal:
                bc_info.type = BCType::Wall;
                break;
            case BCWallInviscid:
                bc_info.type = BCType::SlipWall;
                break;
            case BCSymmetryPlane:
                bc_info.type = BCType::Symmetry;
                break;
            case BCFarfield:
                bc_info.type = BCType::FarField;
                break;
            case BCInflow:
            case BCInflowSubsonic:
                bc_info.type = BCType::Inflow;
                break;
            case BCOutflow:
            case BCOutflowSubsonic:
                bc_info.type = BCType::Outflow;
                break;
            default:
                bc_info.type = BCType::Wall;  // Default
                std::cerr << "Unknown BC type " << cgns_bc_type
                          << " for BC '" << bc_name << "', defaulting to Wall" << std::endl;
        }

        mesh.addBCInfo(bc_id, bc_info);

        // Read point list for this BC
        std::vector<cgsize_t> points(n_pts);
        ierr = cg_boco_read(file_id, base_id, zone_id, bc_id,
                            points.data(), nullptr);
        if (ierr != CG_OK) {
            error_msg_ = "Failed to read BC points: " + std::string(cg_get_error());
            continue;
        }

        // Mark faces with this boundary condition
        // Note: Points are typically node indices on the boundary
        // We need to find the corresponding faces
        // For now, we'll handle this after face building
    }

    return true;
}

}  // namespace zhijian
