#pragma once

#include "common/types.hpp"
#include "mesh/mesh.hpp"
#include <mpi.h>
#include <vector>
#include <map>

namespace zhijian {

// MPI communication patterns for the FR solver
class MPICommunicator {
public:
    MPICommunicator();
    ~MPICommunicator();

    // Initialize MPI (call before using any MPI functions)
    static void init(int* argc, char*** argv);

    // Finalize MPI
    static void finalize();

    // Get rank and size
    int rank() const { return rank_; }
    int size() const { return size_; }

    // Is this the root process?
    bool isRoot() const { return rank_ == 0; }

    // Setup communication pattern based on mesh partitioning
    void setupCommunication(const Mesh& mesh,
                            const std::vector<Index>& local_to_global_elem);

    // Exchange halo (ghost) element data
    // data: solution data at flux points [n_local_elem * n_fp * N_VARS]
    // After exchange, ghost element data is appended
    void exchangeHalo(std::vector<Real>& data, int vars_per_elem);

    // Non-blocking exchange start
    void exchangeHaloBegin(const std::vector<Real>& send_data, int vars_per_elem);

    // Non-blocking exchange finish
    void exchangeHaloEnd(std::vector<Real>& recv_data);

    // GPU-aware MPI exchange (if supported)
    void exchangeHaloGPU(Real* d_data, int vars_per_elem);

    // Barrier synchronization
    void barrier() const { MPI_Barrier(comm_); }

    // Reduction operations
    Real reduceSum(Real local_value) const;
    Real reduceMax(Real local_value) const;
    Real reduceMin(Real local_value) const;

    // All-reduce operations
    void allReduceSum(Real* data, int count) const;
    void allReduceSum(std::vector<Real>& data) const;

    // Broadcast from root
    void broadcast(Real* data, int count, int root = 0) const;
    void broadcast(int* data, int count, int root = 0) const;
    void broadcast(std::vector<Real>& data, int root = 0) const;

    // Gather to root
    void gather(const std::vector<Real>& send_data,
                std::vector<Real>& recv_data, int root = 0) const;

    // Scatter from root
    void scatter(const std::vector<Real>& send_data,
                 std::vector<Real>& recv_data, int root = 0) const;

    // Get neighbor ranks
    const std::vector<int>& neighborRanks() const { return neighbor_ranks_; }

    // Get number of elements to send to each neighbor
    const std::map<int, std::vector<Index>>& sendMap() const { return send_map_; }

    // Get number of elements to receive from each neighbor
    const std::map<int, std::vector<Index>>& recvMap() const { return recv_map_; }

private:
    MPI_Comm comm_;
    int rank_;
    int size_;

    // Communication pattern
    std::vector<int> neighbor_ranks_;
    std::map<int, std::vector<Index>> send_map_;  // Elements to send to each neighbor
    std::map<int, std::vector<Index>> recv_map_;  // Elements to receive from each neighbor

    // Send/receive buffers
    std::vector<Real> send_buffer_;
    std::vector<Real> recv_buffer_;

    // Request handles for non-blocking communication
    std::vector<MPI_Request> send_requests_;
    std::vector<MPI_Request> recv_requests_;
};

// Parallel mesh distribution
class ParallelMesh {
public:
    ParallelMesh() = default;

    // Distribute mesh across MPI processes
    // On root: reads global mesh and partitions
    // On all: receives local partition
    void distribute(const std::string& mesh_file,
                    MPICommunicator& comm);

    // Get local mesh
    const Mesh& localMesh() const { return local_mesh_; }
    Mesh& localMesh() { return local_mesh_; }

    // Mapping from local to global indices
    const std::vector<Index>& localToGlobalElem() const { return local_to_global_elem_; }
    const std::vector<Index>& localToGlobalNode() const { return local_to_global_node_; }

    // Number of ghost elements
    int numGhostElements() const { return num_ghost_; }

private:
    Mesh local_mesh_;
    std::vector<Index> local_to_global_elem_;
    std::vector<Index> local_to_global_node_;
    int num_ghost_ = 0;
};

}  // namespace zhijian
