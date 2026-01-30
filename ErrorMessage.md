make -j
[  6%] Building C object CMakeFiles/hdf5_compat.dir/src/hdf5_compat.c.o
[ 20%] Building CXX object CMakeFiles/zhijian_lib.dir/src/mesh/mesh.cpp.o
[ 26%] Building CXX object CMakeFiles/zhijian_lib.dir/src/mesh/partitioner.cpp.o
[ 40%] Building CXX object CMakeFiles/zhijian_lib.dir/src/io/restart.cpp.o
[ 40%] Building CXX object CMakeFiles/zhijian_lib.dir/src/mesh/cgns_reader.cpp.o
[ 40%] Building CXX object CMakeFiles/zhijian_lib.dir/src/parallel/mpi_comm.cpp.o
[ 53%] Building CUDA object CMakeFiles/zhijian_lib.dir/src/gpu/kernels.cu.o
[ 60%] Building CXX object CMakeFiles/zhijian_lib.dir/src/basis/basis.cpp.o
[ 66%] Building CXX object CMakeFiles/zhijian_lib.dir/src/solver/fr_solver.cpp.o
[ 66%] Building CXX object CMakeFiles/zhijian_lib.dir/src/io/visualization.cpp.o
In file included from /kuhpc/sw/hdf5/1.14.5/gcc/11.4/include/H5public.h:31,
                 from /kuhpc/sw/hdf5/1.14.5/gcc/11.4/include/hdf5.h:21,
                 from /home/z651w035/codes/Zhijian/src/hdf5_compat.c:11:
/kuhpc/sw/hdf5/1.14.5/gcc/11.4/include/H5version.h:869:22: error: conflicting types for ‘H5Literate2’; have ‘herr_t(hid_t,  H5_index_t,  H5_iter_order_t,  hsize_t *, herr_t (*)(hid_t,  const char *, const H5L_info1_t *, void *), void *)’ {aka ‘int(long int,  H5_index_t,  H5_iter_order_t,  long unsigned int *, int (*)(long int,  const char *, const H5L_info1_t *, void *), void *)’}
  869 |   #define H5Literate H5Literate2
      |                      ^~~~~~~~~~~
/home/z651w035/codes/Zhijian/src/hdf5_compat.c:16:8: note: in expansion of macro ‘H5Literate’
   16 | herr_t H5Literate(hid_t grp_id, H5_index_t idx_type, H5_iter_order_t order,
      |        ^~~~~~~~~~
In file included from /kuhpc/sw/hdf5/1.14.5/gcc/11.4/include/H5Gpublic.h:26,
                 from /kuhpc/sw/hdf5/1.14.5/gcc/11.4/include/hdf5.h:29,
                 from /home/z651w035/codes/Zhijian/src/hdf5_compat.c:11:
/kuhpc/sw/hdf5/1.14.5/gcc/11.4/include/H5Lpublic.h:923:15: note: previous declaration of ‘H5Literate2’ with type ‘herr_t(hid_t,  H5_index_t,  H5_iter_order_t,  hsize_t *, herr_t (*)(hid_t,  const char *, const H5L_info2_t *, void *), void *)’ {aka ‘int(long int,  H5_index_t,  H5_iter_order_t,  long unsigned int *, int (*)(long int,  const char *, const H5L_info2_t *, void *), void *)’}
  923 | H5_DLL herr_t H5Literate2(hid_t grp_id, H5_index_t idx_type, H5_iter_order_t order, hsize_t *idx,
      |               ^~~~~~~~~~~
make[2]: *** [CMakeFiles/hdf5_compat.dir/build.make:76: CMakeFiles/hdf5_compat.dir/src/hdf5_compat.c.o] Error 1
make[1]: *** [CMakeFiles/Makefile2:113: CMakeFiles/hdf5_compat.dir/all] Error 2
make[1]: *** Waiting for unfinished jobs....
/home/z651w035/codes/Zhijian/src/gpu/kernels.cu(125): warning #177-D: variable "var" was declared but never referenced
      int var = idx % N_VARS;
          ^

Remark: The warnings can be suppressed with "-diag-suppress <warning-number>"

/home/z651w035/codes/Zhijian/src/gpu/kernels.cu(805): error: no instance of overloaded function "atomicAdd" matches the argument list
            argument types are: (zhijian::Real *__restrict__, zhijian::Real)
          atomicAdd(norm_sq, shared[0]);
          ^
/kuhpc/sw/nvhpc/Linux_x86_64/25.3/cuda/12.8/targets/x86_64-linux/include/sm_20_atomic_functions.hpp(82): note #3326-D: function "atomicAdd(float *, float)" does not match because argument #1 does not match parameter
  static __inline__ __attribute__((device)) float atomicAdd(float *address, float val)
                                                  ^
/kuhpc/sw/nvhpc/Linux_x86_64/25.3/cuda/12.8/targets/x86_64-linux/include/device_atomic_functions.hpp(224): note #3326-D: function "atomicAdd(unsigned long long *, unsigned long long)" does not match because argument #1 does not match parameter
  static __inline__ __attribute__((device)) unsigned long long int atomicAdd(unsigned long long int *address, unsigned long long int val)
                                                                   ^
/kuhpc/sw/nvhpc/Linux_x86_64/25.3/cuda/12.8/targets/x86_64-linux/include/device_atomic_functions.hpp(110): note #3326-D: function "atomicAdd(unsigned int *, unsigned int)" does not match because argument #1 does not match parameter
  static __inline__ __attribute__((device)) unsigned int atomicAdd(unsigned int *address, unsigned int val)
                                                         ^
/kuhpc/sw/nvhpc/Linux_x86_64/25.3/cuda/12.8/targets/x86_64-linux/include/device_atomic_functions.hpp(105): note #3326-D: function "atomicAdd(int *, int)" does not match because argument #1 does not match parameter
  static __inline__ __attribute__((device)) int atomicAdd(int *address, int val)
                                                ^

1 error detected in the compilation of "/home/z651w035/codes/Zhijian/src/gpu/kernels.cu".
make[2]: *** [CMakeFiles/zhijian_lib.dir/build.make:189: CMakeFiles/zhijian_lib.dir/src/gpu/kernels.cu.o] Error 1
make[2]: *** Waiting for unfinished jobs....
make[1]: *** [CMakeFiles/Makefile2:87: CMakeFiles/zhijian_lib.dir/all] Error 2
make: *** [Makefile:136: all] Error 2

