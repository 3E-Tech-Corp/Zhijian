#pragma once

#include "common/types.hpp"
#include "mesh/mesh.hpp"
#include "solver/fr_solver.hpp"
#include <string>
#include <vector>
#include <fstream>

namespace zhijian {

// VTK file writer for visualization
// Supports both legacy VTK and VTU (XML) formats
class VTKWriter {
public:
    VTKWriter() = default;

    // Write solution to VTK file (legacy format)
    bool writeLegacy(const std::string& filename,
                     const Mesh& mesh,
                     const std::vector<ElementSolution>& solution,
                     const SimParams& params);

    // Write solution to VTU file (XML format)
    bool writeVTU(const std::string& filename,
                  const Mesh& mesh,
                  const std::vector<ElementSolution>& solution,
                  const SimParams& params);

    // Write parallel VTU files with PVTU metadata
    bool writeParallelVTU(const std::string& base_filename,
                          int rank, int num_procs,
                          const Mesh& local_mesh,
                          const std::vector<ElementSolution>& solution,
                          const SimParams& params);

    // Write PVTU metadata file (rank 0 only)
    bool writePVTU(const std::string& filename,
                   int num_procs,
                   const std::vector<std::string>& piece_filenames);

    // Set output options
    void setBinary(bool binary) { binary_ = binary; }
    void setHighOrderOutput(bool high_order) { high_order_ = high_order; }
    void setSubdivisions(int n) { subdivisions_ = n; }

    // Get error message
    const std::string& errorMessage() const { return error_msg_; }

private:
    std::string error_msg_;
    bool binary_ = false;
    bool high_order_ = true;
    int subdivisions_ = 4;  // Number of subdivisions per element for visualization

    // Internal helper functions
    void writeHeader(std::ofstream& file, const std::string& title);
    void writePoints(std::ofstream& file, const Mesh& mesh,
                     const std::vector<ElementSolution>& solution,
                     const SimParams& params);
    void writeCells(std::ofstream& file, const Mesh& mesh,
                    const std::vector<ElementSolution>& solution);
    void writePointData(std::ofstream& file, const Mesh& mesh,
                        const std::vector<ElementSolution>& solution,
                        const SimParams& params);

    // Subdivide element for high-resolution output
    void subdivideElement(const Mesh& mesh, Index elem_id,
                          const ElementSolution& sol,
                          const SimParams& params,
                          std::vector<Vec2>& points,
                          std::vector<std::array<int, 3>>& triangles,
                          std::vector<State>& states);
};

// Tecplot file writer
class TecplotWriter {
public:
    TecplotWriter() = default;

    // Write solution to Tecplot ASCII file
    bool writeASCII(const std::string& filename,
                    const Mesh& mesh,
                    const std::vector<ElementSolution>& solution,
                    const SimParams& params);

    // Write solution to Tecplot binary file (PLT format)
    bool writeBinary(const std::string& filename,
                     const Mesh& mesh,
                     const std::vector<ElementSolution>& solution,
                     const SimParams& params);

    // Get error message
    const std::string& errorMessage() const { return error_msg_; }

private:
    std::string error_msg_;
};

// Surface data output (forces, Cp, Cf, etc.)
class SurfaceDataWriter {
public:
    SurfaceDataWriter() = default;

    // Compute and write surface data
    bool write(const std::string& filename,
               const Mesh& mesh,
               const std::vector<ElementSolution>& solution,
               const SimParams& params,
               int bc_tag);  // Which boundary to output

    // Compute aerodynamic coefficients
    void computeAeroCoefficients(const Mesh& mesh,
                                 const std::vector<ElementSolution>& solution,
                                 const SimParams& params,
                                 int bc_tag,
                                 Real& Cl, Real& Cd, Real& Cm);

    // Get error message
    const std::string& errorMessage() const { return error_msg_; }

private:
    std::string error_msg_;
};

// Residual history writer
class ResidualWriter {
public:
    ResidualWriter() = default;

    // Open residual file
    bool open(const std::string& filename);

    // Write residual entry
    void write(int iteration, Real time, Real residual,
               const std::array<Real, N_VARS>& component_residuals);

    // Close file
    void close();

private:
    std::ofstream file_;
};

}  // namespace zhijian
