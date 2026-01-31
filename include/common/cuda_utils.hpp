#pragma once

#include <cuda_runtime.h>
#include <iostream>
#include <cstdlib>

namespace zhijian {

// atomicAdd for double — only compiled by nvcc (CUDA files)
#ifdef __CUDACC__
#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ < 600
static __inline__ __device__ double atomicAdd_double(double* address, double val) {
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val + __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}
#else
static __inline__ __device__ double atomicAdd_double(double* address, double val) {
    return atomicAdd(address, val);
}
#endif
#endif // __CUDACC__

// CUDA error checking macro
#define CUDA_CHECK(call)                                                       \
    do {                                                                        \
        cudaError_t err = call;                                                 \
        if (err != cudaSuccess) {                                               \
            std::cerr << "CUDA error in " << __FILE__ << ":" << __LINE__        \
                      << " - " << cudaGetErrorString(err) << std::endl;         \
            std::exit(EXIT_FAILURE);                                            \
        }                                                                       \
    } while (0)

// Kernel launch macro with error checking
#define CUDA_LAUNCH(kernel, grid, block, ...)                                  \
    do {                                                                        \
        kernel<<<grid, block>>>(__VA_ARGS__);                                   \
        CUDA_CHECK(cudaGetLastError());                                         \
    } while (0)

// Device memory management
template<typename T>
class DeviceArray {
public:
    DeviceArray() : data_(nullptr), size_(0) {}

    explicit DeviceArray(size_t n) : size_(n) {
        if (n > 0) {
            CUDA_CHECK(cudaMalloc(&data_, n * sizeof(T)));
        }
    }

    ~DeviceArray() {
        if (data_) {
            cudaFree(data_);
        }
    }

    // Disable copy
    DeviceArray(const DeviceArray&) = delete;
    DeviceArray& operator=(const DeviceArray&) = delete;

    // Enable move
    DeviceArray(DeviceArray&& other) noexcept
        : data_(other.data_), size_(other.size_) {
        other.data_ = nullptr;
        other.size_ = 0;
    }

    DeviceArray& operator=(DeviceArray&& other) noexcept {
        if (this != &other) {
            if (data_) cudaFree(data_);
            data_ = other.data_;
            size_ = other.size_;
            other.data_ = nullptr;
            other.size_ = 0;
        }
        return *this;
    }

    void resize(size_t n) {
        if (data_) cudaFree(data_);
        size_ = n;
        if (n > 0) {
            CUDA_CHECK(cudaMalloc(&data_, n * sizeof(T)));
        } else {
            data_ = nullptr;
        }
    }

    void copyToDevice(const T* host_data, size_t n) {
        CUDA_CHECK(cudaMemcpy(data_, host_data, n * sizeof(T), cudaMemcpyHostToDevice));
    }

    void copyToDevice(const std::vector<T>& host_data) {
        copyToDevice(host_data.data(), host_data.size());
    }

    void copyToHost(T* host_data, size_t n) const {
        CUDA_CHECK(cudaMemcpy(host_data, data_, n * sizeof(T), cudaMemcpyDeviceToHost));
    }

    void copyToHost(std::vector<T>& host_data) const {
        host_data.resize(size_);
        copyToHost(host_data.data(), size_);
    }

    void fill(T value) {
        CUDA_CHECK(cudaMemset(data_, 0, size_ * sizeof(T)));
    }

    T* data() { return data_; }
    const T* data() const { return data_; }
    size_t size() const { return size_; }

private:
    T* data_;
    size_t size_;
};

// Pinned (page-locked) host memory for faster transfers
template<typename T>
class PinnedArray {
public:
    PinnedArray() : data_(nullptr), size_(0) {}

    explicit PinnedArray(size_t n) : size_(n) {
        if (n > 0) {
            CUDA_CHECK(cudaMallocHost(&data_, n * sizeof(T)));
        }
    }

    ~PinnedArray() {
        if (data_) {
            cudaFreeHost(data_);
        }
    }

    PinnedArray(const PinnedArray&) = delete;
    PinnedArray& operator=(const PinnedArray&) = delete;

    PinnedArray(PinnedArray&& other) noexcept
        : data_(other.data_), size_(other.size_) {
        other.data_ = nullptr;
        other.size_ = 0;
    }

    PinnedArray& operator=(PinnedArray&& other) noexcept {
        if (this != &other) {
            if (data_) cudaFreeHost(data_);
            data_ = other.data_;
            size_ = other.size_;
            other.data_ = nullptr;
            other.size_ = 0;
        }
        return *this;
    }

    void resize(size_t n) {
        if (data_) cudaFreeHost(data_);
        size_ = n;
        if (n > 0) {
            CUDA_CHECK(cudaMallocHost(&data_, n * sizeof(T)));
        } else {
            data_ = nullptr;
        }
    }

    T& operator[](size_t i) { return data_[i]; }
    const T& operator[](size_t i) const { return data_[i]; }

    T* data() { return data_; }
    const T* data() const { return data_; }
    size_t size() const { return size_; }

private:
    T* data_;
    size_t size_;
};

// CUDA stream wrapper
class CudaStream {
public:
    CudaStream() {
        CUDA_CHECK(cudaStreamCreate(&stream_));
    }

    ~CudaStream() {
        cudaStreamDestroy(stream_);
    }

    CudaStream(const CudaStream&) = delete;
    CudaStream& operator=(const CudaStream&) = delete;

    void synchronize() {
        CUDA_CHECK(cudaStreamSynchronize(stream_));
    }

    cudaStream_t get() const { return stream_; }

private:
    cudaStream_t stream_;
};

// Get optimal block size for kernel
inline int getOptimalBlockSize(int n, int max_threads = 256) {
    return std::min(max_threads, (n + 31) / 32 * 32);
}

// Get grid size for n elements
inline int getGridSize(int n, int block_size) {
    return (n + block_size - 1) / block_size;
}

}  // namespace zhijian
