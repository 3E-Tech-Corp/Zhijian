#pragma once

#include "common/types.hpp"
#include "mesh/mesh.hpp"
#include "solver/fr_solver.hpp"
#include <string>
#include <vector>

namespace zhijian {

// HDF5 restart file writer/reader
class RestartIO {
public:
    RestartIO() = default;

    // Write restart file
    // Includes mesh connectivity, solution data, and simulation state
    bool write(const std::string& filename,
               const Mesh& mesh,
               const std::vector<ElementSolution>& solution,
               const SimParams& params,
               Real time,
               int iteration);

    // Read restart file
    bool read(const std::string& filename,
              Mesh& mesh,
              std::vector<ElementSolution>& solution,
              SimParams& params,
              Real& time,
              int& iteration);

    // Write parallel restart file (each rank writes its portion)
    bool writeParallel(const std::string& filename,
                       int rank, int num_procs,
                       const Mesh& local_mesh,
                       const std::vector<ElementSolution>& solution,
                       const SimParams& params,
                       Real time,
                       int iteration);

    // Read parallel restart file
    bool readParallel(const std::string& filename,
                      int rank, int num_procs,
                      Mesh& local_mesh,
                      std::vector<ElementSolution>& solution,
                      SimParams& params,
                      Real& time,
                      int& iteration);

    // Get error message
    const std::string& errorMessage() const { return error_msg_; }

private:
    std::string error_msg_;

    // Internal helper functions
    bool writeDataset(int file_id, const std::string& name,
                      const std::vector<Real>& data,
                      const std::vector<size_t>& dims);
    bool readDataset(int file_id, const std::string& name,
                     std::vector<Real>& data,
                     std::vector<size_t>& dims);
    bool writeAttribute(int file_id, const std::string& name, Real value);
    bool writeAttribute(int file_id, const std::string& name, int value);
    bool writeAttribute(int file_id, const std::string& name, const std::string& value);
    bool readAttribute(int file_id, const std::string& name, Real& value);
    bool readAttribute(int file_id, const std::string& name, int& value);
    bool readAttribute(int file_id, const std::string& name, std::string& value);
};

}  // namespace zhijian
