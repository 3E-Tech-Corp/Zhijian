#include "io/visualization.hpp"
#include "basis/basis.hpp"
#include <iomanip>
#include <sstream>

namespace zhijian {

bool VTKWriter::writeLegacy(const std::string& filename,
                             const Mesh& mesh,
                             const std::vector<ElementSolution>& solution,
                             const SimParams& params) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        error_msg_ = "Failed to open file: " + filename;
        return false;
    }

    // Collect all points and cells from subdivided elements
    std::vector<Vec2> all_points;
    std::vector<std::array<int, 3>> all_triangles;
    std::vector<State> all_states;

    for (Index i = 0; i < mesh.numElements(); ++i) {
        std::vector<Vec2> elem_points;
        std::vector<std::array<int, 3>> elem_tris;
        std::vector<State> elem_states;

        subdivideElement(mesh, i, solution[i], params, elem_points, elem_tris, elem_states);

        int offset = static_cast<int>(all_points.size());
        for (const auto& p : elem_points) {
            all_points.push_back(p);
        }
        for (const auto& s : elem_states) {
            all_states.push_back(s);
        }
        for (auto tri : elem_tris) {
            tri[0] += offset;
            tri[1] += offset;
            tri[2] += offset;
            all_triangles.push_back(tri);
        }
    }

    // Write VTK header
    file << "# vtk DataFile Version 3.0\n";
    file << "Zhijian FR Solver Output\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    // Write points
    file << "POINTS " << all_points.size() << " double\n";
    for (const auto& p : all_points) {
        file << std::scientific << std::setprecision(12);
        file << p.x << " " << p.y << " 0.0\n";
    }

    // Write cells
    int n_cells = static_cast<int>(all_triangles.size());
    file << "CELLS " << n_cells << " " << n_cells * 4 << "\n";
    for (const auto& tri : all_triangles) {
        file << "3 " << tri[0] << " " << tri[1] << " " << tri[2] << "\n";
    }

    // Cell types (5 = VTK_TRIANGLE)
    file << "CELL_TYPES " << n_cells << "\n";
    for (int i = 0; i < n_cells; ++i) {
        file << "5\n";
    }

    // Write point data
    file << "POINT_DATA " << all_points.size() << "\n";

    IdealGas gas(params.gamma);

    // Density
    file << "SCALARS Density double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& s : all_states) {
        file << s.rho() << "\n";
    }

    // Pressure
    file << "SCALARS Pressure double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& s : all_states) {
        file << gas.pressure(s) << "\n";
    }

    // Mach number
    file << "SCALARS Mach double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& s : all_states) {
        file << gas.machNumber(s) << "\n";
    }

    // Velocity
    file << "VECTORS Velocity double\n";
    for (const auto& s : all_states) {
        Real u = s.rhou() / s.rho();
        Real v = s.rhov() / s.rho();
        file << u << " " << v << " 0.0\n";
    }

    file.close();
    return true;
}

bool VTKWriter::writeVTU(const std::string& filename,
                          const Mesh& mesh,
                          const std::vector<ElementSolution>& solution,
                          const SimParams& params) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        error_msg_ = "Failed to open file: " + filename;
        return false;
    }

    // Collect all data
    std::vector<Vec2> all_points;
    std::vector<std::array<int, 3>> all_triangles;
    std::vector<State> all_states;

    for (Index i = 0; i < mesh.numElements(); ++i) {
        std::vector<Vec2> elem_points;
        std::vector<std::array<int, 3>> elem_tris;
        std::vector<State> elem_states;

        subdivideElement(mesh, i, solution[i], params, elem_points, elem_tris, elem_states);

        int offset = static_cast<int>(all_points.size());
        for (const auto& p : elem_points) all_points.push_back(p);
        for (const auto& s : elem_states) all_states.push_back(s);
        for (auto tri : elem_tris) {
            tri[0] += offset;
            tri[1] += offset;
            tri[2] += offset;
            all_triangles.push_back(tri);
        }
    }

    int n_points = static_cast<int>(all_points.size());
    int n_cells = static_cast<int>(all_triangles.size());

    // VTU XML format
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";
    file << "    <Piece NumberOfPoints=\"" << n_points << "\" NumberOfCells=\"" << n_cells << "\">\n";

    // Points
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& p : all_points) {
        file << "          " << p.x << " " << p.y << " 0.0\n";
    }
    file << "        </DataArray>\n";
    file << "      </Points>\n";

    // Cells
    file << "      <Cells>\n";
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& tri : all_triangles) {
        file << "          " << tri[0] << " " << tri[1] << " " << tri[2] << "\n";
    }
    file << "        </DataArray>\n";
    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 0; i < n_cells; ++i) {
        file << "          " << (i + 1) * 3 << "\n";
    }
    file << "        </DataArray>\n";
    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < n_cells; ++i) {
        file << "          5\n";  // VTK_TRIANGLE
    }
    file << "        </DataArray>\n";
    file << "      </Cells>\n";

    // Point data
    IdealGas gas(params.gamma);

    file << "      <PointData Scalars=\"Density\" Vectors=\"Velocity\">\n";

    // Density
    file << "        <DataArray type=\"Float64\" Name=\"Density\" format=\"ascii\">\n";
    for (const auto& s : all_states) {
        file << "          " << s.rho() << "\n";
    }
    file << "        </DataArray>\n";

    // Pressure
    file << "        <DataArray type=\"Float64\" Name=\"Pressure\" format=\"ascii\">\n";
    for (const auto& s : all_states) {
        file << "          " << gas.pressure(s) << "\n";
    }
    file << "        </DataArray>\n";

    // Mach
    file << "        <DataArray type=\"Float64\" Name=\"Mach\" format=\"ascii\">\n";
    for (const auto& s : all_states) {
        file << "          " << gas.machNumber(s) << "\n";
    }
    file << "        </DataArray>\n";

    // Velocity
    file << "        <DataArray type=\"Float64\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& s : all_states) {
        Real u = s.rhou() / s.rho();
        Real v = s.rhov() / s.rho();
        file << "          " << u << " " << v << " 0.0\n";
    }
    file << "        </DataArray>\n";

    file << "      </PointData>\n";

    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";

    file.close();
    return true;
}

bool VTKWriter::writeParallelVTU(const std::string& base_filename,
                                  int rank, int num_procs,
                                  const Mesh& local_mesh,
                                  const std::vector<ElementSolution>& solution,
                                  const SimParams& params) {
    // Write local VTU file
    std::string local_filename = base_filename + "_" + std::to_string(rank) + ".vtu";
    if (!writeVTU(local_filename, local_mesh, solution, params)) {
        return false;
    }

    // Rank 0 writes PVTU file
    if (rank == 0) {
        std::vector<std::string> piece_files;
        for (int p = 0; p < num_procs; ++p) {
            piece_files.push_back(base_filename + "_" + std::to_string(p) + ".vtu");
        }
        return writePVTU(base_filename + ".pvtu", num_procs, piece_files);
    }

    return true;
}

bool VTKWriter::writePVTU(const std::string& filename,
                           int num_procs,
                           const std::vector<std::string>& piece_filenames) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        error_msg_ = "Failed to open file: " + filename;
        return false;
    }

    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <PUnstructuredGrid GhostLevel=\"0\">\n";

    // Point data
    file << "    <PPointData Scalars=\"Density\" Vectors=\"Velocity\">\n";
    file << "      <PDataArray type=\"Float64\" Name=\"Density\"/>\n";
    file << "      <PDataArray type=\"Float64\" Name=\"Pressure\"/>\n";
    file << "      <PDataArray type=\"Float64\" Name=\"Mach\"/>\n";
    file << "      <PDataArray type=\"Float64\" Name=\"Velocity\" NumberOfComponents=\"3\"/>\n";
    file << "    </PPointData>\n";

    // Points
    file << "    <PPoints>\n";
    file << "      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>\n";
    file << "    </PPoints>\n";

    // Pieces
    for (const auto& pf : piece_filenames) {
        file << "    <Piece Source=\"" << pf << "\"/>\n";
    }

    file << "  </PUnstructuredGrid>\n";
    file << "</VTKFile>\n";

    file.close();
    return true;
}

void VTKWriter::subdivideElement(const Mesh& mesh, Index elem_id,
                                  const ElementSolution& sol,
                                  const SimParams& params,
                                  std::vector<Vec2>& points,
                                  std::vector<std::array<int, 3>>& triangles,
                                  std::vector<State>& states) {
    const Element& elem = mesh.element(elem_id);
    int n = subdivisions_;

    // Get element nodes
    std::vector<Vec2> nodes(elem.node_ids.size());
    for (size_t j = 0; j < elem.node_ids.size(); ++j) {
        nodes[j] = mesh.node(elem.node_ids[j]);
    }

    if (elem.type == ElementType::Quadrilateral) {
        // Create (n+1) x (n+1) grid of points
        int n_pts = (n + 1) * (n + 1);
        points.resize(n_pts);
        states.resize(n_pts);

        // Get FR operators for interpolation
        FROperators ops;
        ops.init(elem.type, params.poly_order, params.flux_type);

        for (int j = 0; j <= n; ++j) {
            for (int i = 0; i <= n; ++i) {
                int idx = j * (n + 1) + i;

                // Reference coordinates
                Real xi = -1.0 + 2.0 * i / n;
                Real eta = -1.0 + 2.0 * j / n;

                // Physical coordinates
                points[idx] = ElementGeometry::refToPhys(nodes, elem.type, elem.order, xi, eta);

                // Interpolate solution
                states[idx] = State();
                std::vector<Real> basis_vals;
                // Evaluate basis functions and interpolate (simplified)
                // For now, just use cell-average
                if (!sol.sol_pts.empty()) {
                    states[idx] = sol.sol_pts[0];  // Simplified: use first SP
                }
            }
        }

        // Create triangles (2 per quad cell)
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                int v00 = j * (n + 1) + i;
                int v10 = j * (n + 1) + i + 1;
                int v01 = (j + 1) * (n + 1) + i;
                int v11 = (j + 1) * (n + 1) + i + 1;

                triangles.push_back({v00, v10, v11});
                triangles.push_back({v00, v11, v01});
            }
        }
    } else {
        // Triangle element
        // Create barycentric subdivision
        int n_pts = (n + 1) * (n + 2) / 2;
        points.resize(n_pts);
        states.resize(n_pts);

        int idx = 0;
        for (int j = 0; j <= n; ++j) {
            for (int i = 0; i <= n - j; ++i) {
                Real xi = static_cast<Real>(i) / n;
                Real eta = static_cast<Real>(j) / n;

                points[idx] = ElementGeometry::refToPhys(nodes, elem.type, elem.order, xi, eta);

                if (!sol.sol_pts.empty()) {
                    states[idx] = sol.sol_pts[0];
                }
                idx++;
            }
        }

        // Create triangles
        idx = 0;
        for (int j = 0; j < n; ++j) {
            int row_start = j * (n + 1) - j * (j - 1) / 2;
            int next_row_start = (j + 1) * (n + 1) - (j + 1) * j / 2;

            for (int i = 0; i < n - j; ++i) {
                int v0 = row_start + i;
                int v1 = row_start + i + 1;
                int v2 = next_row_start + i;

                triangles.push_back({v0, v1, v2});

                if (i < n - j - 1) {
                    int v3 = next_row_start + i + 1;
                    triangles.push_back({v1, v3, v2});
                }
            }
        }
    }
}

// ============================================================================
// Residual Writer
// ============================================================================

bool ResidualWriter::open(const std::string& filename) {
    file_.open(filename);
    if (!file_.is_open()) {
        return false;
    }

    file_ << "# Iteration, Time, Residual, Rho_res, RhoU_res, RhoV_res, RhoE_res\n";
    return true;
}

void ResidualWriter::write(int iteration, Real time, Real residual,
                            const std::array<Real, N_VARS>& component_residuals) {
    if (file_.is_open()) {
        file_ << iteration << ", " << std::scientific << std::setprecision(8)
              << time << ", " << residual;
        for (int v = 0; v < N_VARS; ++v) {
            file_ << ", " << component_residuals[v];
        }
        file_ << "\n";
        file_.flush();
    }
}

void ResidualWriter::close() {
    if (file_.is_open()) {
        file_.close();
    }
}

}  // namespace zhijian
