#pragma once

// CUDA device-only utilities.
// This header must ONLY be included from .cu files (compiled by nvcc).

#include "cuda_utils.hpp"

namespace zhijian {

// atomicAdd for double — required for compute capability < 6.0,
// and for compilers that don't provide the built-in double overload.
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

}  // namespace zhijian
