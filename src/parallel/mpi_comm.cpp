#include "parallel/mpi_comm.hpp"
#include <iostream>
#include <algorithm>
#include <numeric>

namespace zhijian {

MPICommunicator::MPICommunicator() : comm_(MPI_COMM_WORLD), rank_(0), size_(1) {
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &size_);
}

MPICommunicator::~MPICommunicator() = default;

void MPICommunicator::init(int* argc, char*** argv) {
    int initialized;
    MPI_Initialized(&initialized);
    if (!initialized) {
        MPI_Init(argc, argv);
    }
}

void MPICommunicator::finalize() {
    int finalized;
    MPI_Finalized(&finalized);
    if (!finalized) {
        MPI_Finalize();
    }
}

void MPICommunicator::setupCommunication(const Mesh& mesh,
                                          const std::vector<Index>& local_to_global_elem) {
    // Find neighbor ranks and elements to exchange
    neighbor_ranks_.clear();
    send_map_.clear();
    recv_map_.clear();

    std::set<int> neighbors;

    for (const auto& face : mesh.faces()) {
        if (!face.is_boundary && face.mpi_rank >= 0 && face.mpi_rank != rank_) {
            neighbors.insert(face.mpi_rank);

            // Element to send to this neighbor
            send_map_[face.mpi_rank].push_back(face.left_elem);
        }
    }

    neighbor_ranks_.assign(neighbors.begin(), neighbors.end());

    // Exchange send/receive sizes with neighbors
    for (int neighbor : neighbor_ranks_) {
        // Send number of elements we'll send to this neighbor
        int send_count = static_cast<int>(send_map_[neighbor].size());
        int recv_count;

        MPI_Sendrecv(&send_count, 1, MPI_INT, neighbor, 0,
                     &recv_count, 1, MPI_INT, neighbor, 0,
                     comm_, MPI_STATUS_IGNORE);

        recv_map_[neighbor].resize(recv_count);
    }

    std::cout << "Rank " << rank_ << ": " << neighbor_ranks_.size()
              << " neighbors" << std::endl;
}

void MPICommunicator::exchangeHalo(std::vector<Real>& data, int vars_per_elem) {
    if (neighbor_ranks_.empty()) return;

    // Prepare send buffers
    std::vector<std::vector<Real>> send_bufs(neighbor_ranks_.size());
    std::vector<std::vector<Real>> recv_bufs(neighbor_ranks_.size());

    for (size_t i = 0; i < neighbor_ranks_.size(); ++i) {
        int neighbor = neighbor_ranks_[i];
        const auto& send_elems = send_map_.at(neighbor);
        const auto& recv_elems = recv_map_.at(neighbor);

        send_bufs[i].resize(send_elems.size() * vars_per_elem);
        recv_bufs[i].resize(recv_elems.size() * vars_per_elem);

        // Pack send data
        for (size_t j = 0; j < send_elems.size(); ++j) {
            Index elem = send_elems[j];
            for (int v = 0; v < vars_per_elem; ++v) {
                send_bufs[i][j * vars_per_elem + v] = data[elem * vars_per_elem + v];
            }
        }
    }

    // Non-blocking sends and receives
    std::vector<MPI_Request> requests(2 * neighbor_ranks_.size());

    for (size_t i = 0; i < neighbor_ranks_.size(); ++i) {
        int neighbor = neighbor_ranks_[i];

        MPI_Isend(send_bufs[i].data(), static_cast<int>(send_bufs[i].size()),
                  MPI_DOUBLE, neighbor, 0, comm_, &requests[2*i]);
        MPI_Irecv(recv_bufs[i].data(), static_cast<int>(recv_bufs[i].size()),
                  MPI_DOUBLE, neighbor, 0, comm_, &requests[2*i + 1]);
    }

    // Wait for all communication to complete
    MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);

    // Unpack receive data (append to ghost region)
    // This assumes ghost elements are at the end of the data array
}

void MPICommunicator::exchangeHaloBegin(const std::vector<Real>& send_data,
                                         int vars_per_elem) {
    if (neighbor_ranks_.empty()) return;

    send_requests_.resize(neighbor_ranks_.size());
    recv_requests_.resize(neighbor_ranks_.size());

    // Prepare and send
    size_t offset = 0;
    for (size_t i = 0; i < neighbor_ranks_.size(); ++i) {
        int neighbor = neighbor_ranks_[i];
        const auto& send_elems = send_map_.at(neighbor);
        const auto& recv_elems = recv_map_.at(neighbor);

        int send_size = static_cast<int>(send_elems.size() * vars_per_elem);
        int recv_size = static_cast<int>(recv_elems.size() * vars_per_elem);

        // Would need proper buffer management here
        MPI_Isend(send_data.data() + offset, send_size, MPI_DOUBLE,
                  neighbor, 0, comm_, &send_requests_[i]);

        offset += send_size;
    }
}

void MPICommunicator::exchangeHaloEnd(std::vector<Real>& recv_data) {
    if (neighbor_ranks_.empty()) return;

    MPI_Waitall(static_cast<int>(send_requests_.size()), send_requests_.data(),
                MPI_STATUSES_IGNORE);
    MPI_Waitall(static_cast<int>(recv_requests_.size()), recv_requests_.data(),
                MPI_STATUSES_IGNORE);
}

void MPICommunicator::exchangeHaloGPU(Real* d_data, int vars_per_elem) {
    // GPU-aware MPI exchange
    // Requires CUDA-aware MPI implementation
    // For now, fall back to CPU transfer
    std::cerr << "GPU-aware MPI not implemented, using CPU fallback" << std::endl;
}

Real MPICommunicator::reduceSum(Real local_value) const {
    Real global_value;
    MPI_Allreduce(&local_value, &global_value, 1, MPI_DOUBLE, MPI_SUM, comm_);
    return global_value;
}

Real MPICommunicator::reduceMax(Real local_value) const {
    Real global_value;
    MPI_Allreduce(&local_value, &global_value, 1, MPI_DOUBLE, MPI_MAX, comm_);
    return global_value;
}

Real MPICommunicator::reduceMin(Real local_value) const {
    Real global_value;
    MPI_Allreduce(&local_value, &global_value, 1, MPI_DOUBLE, MPI_MIN, comm_);
    return global_value;
}

void MPICommunicator::allReduceSum(Real* data, int count) const {
    MPI_Allreduce(MPI_IN_PLACE, data, count, MPI_DOUBLE, MPI_SUM, comm_);
}

void MPICommunicator::allReduceSum(std::vector<Real>& data) const {
    allReduceSum(data.data(), static_cast<int>(data.size()));
}

void MPICommunicator::broadcast(Real* data, int count, int root) const {
    MPI_Bcast(data, count, MPI_DOUBLE, root, comm_);
}

void MPICommunicator::broadcast(int* data, int count, int root) const {
    MPI_Bcast(data, count, MPI_INT, root, comm_);
}

void MPICommunicator::broadcast(std::vector<Real>& data, int root) const {
    int size = static_cast<int>(data.size());
    broadcast(&size, 1, root);
    data.resize(size);
    broadcast(data.data(), size, root);
}

void MPICommunicator::gather(const std::vector<Real>& send_data,
                              std::vector<Real>& recv_data, int root) const {
    int send_count = static_cast<int>(send_data.size());
    std::vector<int> recv_counts(size_);
    std::vector<int> displs(size_);

    MPI_Gather(&send_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, root, comm_);

    if (rank_ == root) {
        displs[0] = 0;
        for (int i = 1; i < size_; ++i) {
            displs[i] = displs[i-1] + recv_counts[i-1];
        }
        recv_data.resize(displs[size_-1] + recv_counts[size_-1]);
    }

    MPI_Gatherv(send_data.data(), send_count, MPI_DOUBLE,
                recv_data.data(), recv_counts.data(), displs.data(),
                MPI_DOUBLE, root, comm_);
}

void MPICommunicator::scatter(const std::vector<Real>& send_data,
                               std::vector<Real>& recv_data, int root) const {
    int recv_count;
    std::vector<int> send_counts(size_);
    std::vector<int> displs(size_);

    if (rank_ == root) {
        int per_proc = static_cast<int>(send_data.size()) / size_;
        for (int i = 0; i < size_; ++i) {
            send_counts[i] = per_proc;
            displs[i] = i * per_proc;
        }
    }

    MPI_Scatter(send_counts.data(), 1, MPI_INT, &recv_count, 1, MPI_INT, root, comm_);
    recv_data.resize(recv_count);

    MPI_Scatterv(send_data.data(), send_counts.data(), displs.data(), MPI_DOUBLE,
                 recv_data.data(), recv_count, MPI_DOUBLE, root, comm_);
}

// ============================================================================
// ParallelMesh Implementation
// ============================================================================

void ParallelMesh::distribute(const std::string& mesh_file, MPICommunicator& comm) {
    Mesh global_mesh;

    // Root reads and partitions the mesh
    if (comm.isRoot()) {
        CGNSReader reader;
        if (!reader.read(mesh_file, global_mesh)) {
            std::cerr << "Failed to read mesh: " << reader.errorMessage() << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        std::cout << "Global mesh loaded:" << std::endl;
        global_mesh.printStats();

        // Partition the mesh
        MeshPartitioner partitioner;
        if (!partitioner.partition(global_mesh, comm.size())) {
            std::cerr << "Failed to partition mesh: " << partitioner.errorMessage() << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // Broadcast mesh size info
    int n_nodes = 0, n_elems = 0;
    if (comm.isRoot()) {
        n_nodes = static_cast<int>(global_mesh.numNodes());
        n_elems = static_cast<int>(global_mesh.numElements());
    }
    comm.broadcast(&n_nodes, 1);
    comm.broadcast(&n_elems, 1);

    // Each rank extracts its local partition
    if (comm.isRoot()) {
        for (int p = 0; p < comm.size(); ++p) {
            Mesh part_mesh;
            std::vector<Index> l2g_elem, l2g_node;

            MeshPartitioner partitioner;
            partitioner.extractLocalMesh(global_mesh, p, part_mesh, l2g_elem, l2g_node);

            if (p == 0) {
                local_mesh_ = std::move(part_mesh);
                local_to_global_elem_ = std::move(l2g_elem);
                local_to_global_node_ = std::move(l2g_node);
            } else {
                // Send to rank p
                // (Simplified - would need full serialization)
            }
        }
    } else {
        // Receive from root
        // (Simplified - would need full serialization)
    }

    // Setup MPI communication pattern
    comm.setupCommunication(local_mesh_, local_to_global_elem_);

    // Count ghost elements
    num_ghost_ = 0;
    for (const auto& face : local_mesh_.faces()) {
        if (face.mpi_rank >= 0 && face.mpi_rank != comm.rank()) {
            num_ghost_++;
        }
    }

    std::cout << "Rank " << comm.rank() << ": "
              << local_mesh_.numElements() << " elements, "
              << num_ghost_ << " ghost elements" << std::endl;
}

}  // namespace zhijian
