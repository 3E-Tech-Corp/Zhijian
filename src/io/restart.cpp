#include "io/restart.hpp"
#include <hdf5.h>
#include <iostream>

namespace zhijian {

bool RestartIO::write(const std::string& filename,
                       const Mesh& mesh,
                       const std::vector<ElementSolution>& solution,
                       const SimParams& params,
                       Real time,
                       int iteration) {
    // Create HDF5 file
    hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) {
        error_msg_ = "Failed to create HDF5 file: " + filename;
        return false;
    }

    // Write attributes
    hid_t attr_space = H5Screate(H5S_SCALAR);

    // Time
    hid_t attr_id = H5Acreate2(file_id, "time", H5T_NATIVE_DOUBLE, attr_space,
                                H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &time);
    H5Aclose(attr_id);

    // Iteration
    attr_id = H5Acreate2(file_id, "iteration", H5T_NATIVE_INT, attr_space,
                          H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_INT, &iteration);
    H5Aclose(attr_id);

    // Polynomial order
    int poly_order = params.poly_order;
    attr_id = H5Acreate2(file_id, "poly_order", H5T_NATIVE_INT, attr_space,
                          H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_INT, &poly_order);
    H5Aclose(attr_id);

    // Gamma
    attr_id = H5Acreate2(file_id, "gamma", H5T_NATIVE_DOUBLE, attr_space,
                          H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &params.gamma);
    H5Aclose(attr_id);

    H5Sclose(attr_space);

    // Write solution data
    size_t n_elem = solution.size();

    // Flatten solution data
    std::vector<Real> sol_data;
    std::vector<int> elem_n_sp;

    for (size_t i = 0; i < n_elem; ++i) {
        int n_sp = static_cast<int>(solution[i].sol_pts.size());
        elem_n_sp.push_back(n_sp);

        for (int sp = 0; sp < n_sp; ++sp) {
            for (int v = 0; v < N_VARS; ++v) {
                sol_data.push_back(solution[i].sol_pts[sp][v]);
            }
        }
    }

    // Write solution array
    hsize_t dims[1] = {sol_data.size()};
    hid_t dataspace = H5Screate_simple(1, dims, nullptr);
    hid_t dataset = H5Dcreate2(file_id, "solution", H5T_NATIVE_DOUBLE, dataspace,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, sol_data.data());
    H5Dclose(dataset);
    H5Sclose(dataspace);

    // Write element info
    dims[0] = elem_n_sp.size();
    dataspace = H5Screate_simple(1, dims, nullptr);
    dataset = H5Dcreate2(file_id, "elem_n_sp", H5T_NATIVE_INT, dataspace,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, elem_n_sp.data());
    H5Dclose(dataset);
    H5Sclose(dataspace);

    H5Fclose(file_id);

    std::cout << "Restart file written: " << filename << std::endl;
    return true;
}

bool RestartIO::read(const std::string& filename,
                      Mesh& mesh,
                      std::vector<ElementSolution>& solution,
                      SimParams& params,
                      Real& time,
                      int& iteration) {
    // Open HDF5 file
    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        error_msg_ = "Failed to open HDF5 file: " + filename;
        return false;
    }

    // Read attributes
    hid_t attr_id = H5Aopen(file_id, "time", H5P_DEFAULT);
    H5Aread(attr_id, H5T_NATIVE_DOUBLE, &time);
    H5Aclose(attr_id);

    attr_id = H5Aopen(file_id, "iteration", H5P_DEFAULT);
    H5Aread(attr_id, H5T_NATIVE_INT, &iteration);
    H5Aclose(attr_id);

    attr_id = H5Aopen(file_id, "poly_order", H5P_DEFAULT);
    H5Aread(attr_id, H5T_NATIVE_INT, &params.poly_order);
    H5Aclose(attr_id);

    attr_id = H5Aopen(file_id, "gamma", H5P_DEFAULT);
    H5Aread(attr_id, H5T_NATIVE_DOUBLE, &params.gamma);
    H5Aclose(attr_id);

    // Read element info
    hid_t dataset = H5Dopen2(file_id, "elem_n_sp", H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);
    hsize_t dims[1];
    H5Sget_simple_extent_dims(dataspace, dims, nullptr);

    std::vector<int> elem_n_sp(dims[0]);
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, elem_n_sp.data());
    H5Dclose(dataset);
    H5Sclose(dataspace);

    // Read solution data
    dataset = H5Dopen2(file_id, "solution", H5P_DEFAULT);
    dataspace = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace, dims, nullptr);

    std::vector<Real> sol_data(dims[0]);
    H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, sol_data.data());
    H5Dclose(dataset);
    H5Sclose(dataspace);

    H5Fclose(file_id);

    // Reconstruct solution
    size_t n_elem = elem_n_sp.size();
    solution.resize(n_elem);

    size_t idx = 0;
    for (size_t i = 0; i < n_elem; ++i) {
        int n_sp = elem_n_sp[i];
        solution[i].sol_pts.resize(n_sp);

        for (int sp = 0; sp < n_sp; ++sp) {
            for (int v = 0; v < N_VARS; ++v) {
                solution[i].sol_pts[sp][v] = sol_data[idx++];
            }
        }
    }

    std::cout << "Restart file read: " << filename
              << " (time=" << time << ", iter=" << iteration << ")" << std::endl;
    return true;
}

bool RestartIO::writeParallel(const std::string& filename,
                               int rank, int num_procs,
                               const Mesh& local_mesh,
                               const std::vector<ElementSolution>& solution,
                               const SimParams& params,
                               Real time,
                               int iteration) {
    // Use parallel HDF5
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

    hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

    if (file_id < 0) {
        error_msg_ = "Failed to create parallel HDF5 file: " + filename;
        return false;
    }

    // Write global attributes (from rank 0)
    if (rank == 0) {
        hid_t attr_space = H5Screate(H5S_SCALAR);

        hid_t attr_id = H5Acreate2(file_id, "time", H5T_NATIVE_DOUBLE, attr_space,
                                    H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &time);
        H5Aclose(attr_id);

        attr_id = H5Acreate2(file_id, "iteration", H5T_NATIVE_INT, attr_space,
                              H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr_id, H5T_NATIVE_INT, &iteration);
        H5Aclose(attr_id);

        attr_id = H5Acreate2(file_id, "num_procs", H5T_NATIVE_INT, attr_space,
                              H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr_id, H5T_NATIVE_INT, &num_procs);
        H5Aclose(attr_id);

        H5Sclose(attr_space);
    }

    // Each rank writes its data to a separate group
    std::string group_name = "/rank_" + std::to_string(rank);
    hid_t group_id = H5Gcreate2(file_id, group_name.c_str(), H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);

    // Flatten and write local solution
    std::vector<Real> sol_data;
    for (const auto& elem_sol : solution) {
        for (const auto& state : elem_sol.sol_pts) {
            for (int v = 0; v < N_VARS; ++v) {
                sol_data.push_back(state[v]);
            }
        }
    }

    hsize_t dims[1] = {sol_data.size()};
    hid_t dataspace = H5Screate_simple(1, dims, nullptr);
    hid_t dataset = H5Dcreate2(group_id, "solution", H5T_NATIVE_DOUBLE, dataspace,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, sol_data.data());
    H5Dclose(dataset);
    H5Sclose(dataspace);

    H5Gclose(group_id);
    H5Fclose(file_id);

    return true;
}

bool RestartIO::readParallel(const std::string& filename,
                              int rank, int num_procs,
                              Mesh& local_mesh,
                              std::vector<ElementSolution>& solution,
                              SimParams& params,
                              Real& time,
                              int& iteration) {
    // Open with parallel HDF5
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, plist_id);
    H5Pclose(plist_id);

    if (file_id < 0) {
        error_msg_ = "Failed to open parallel HDF5 file: " + filename;
        return false;
    }

    // Read global attributes
    hid_t attr_id = H5Aopen(file_id, "time", H5P_DEFAULT);
    H5Aread(attr_id, H5T_NATIVE_DOUBLE, &time);
    H5Aclose(attr_id);

    attr_id = H5Aopen(file_id, "iteration", H5P_DEFAULT);
    H5Aread(attr_id, H5T_NATIVE_INT, &iteration);
    H5Aclose(attr_id);

    // Read local data from this rank's group
    std::string group_name = "/rank_" + std::to_string(rank);
    hid_t group_id = H5Gopen2(file_id, group_name.c_str(), H5P_DEFAULT);

    hid_t dataset = H5Dopen2(group_id, "solution", H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);
    hsize_t dims[1];
    H5Sget_simple_extent_dims(dataspace, dims, nullptr);

    std::vector<Real> sol_data(dims[0]);
    H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, sol_data.data());

    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Gclose(group_id);
    H5Fclose(file_id);

    // Reconstruct solution (needs element info)
    // ...

    return true;
}

}  // namespace zhijian
