make -j
[  6%] Building C object CMakeFiles/hdf5_compat.dir/src/hdf5_compat.c.o
[ 20%] Building CXX object CMakeFiles/zhijian_lib.dir/src/mesh/cgns_reader.cpp.o
[ 20%] Building CXX object CMakeFiles/zhijian_lib.dir/src/mesh/mesh.cpp.o
[ 53%] Building CXX object CMakeFiles/zhijian_lib.dir/src/mesh/partitioner.cpp.o
[ 53%] Building CUDA object CMakeFiles/zhijian_lib.dir/src/gpu/kernels.cu.o
[ 66%] Building CXX object CMakeFiles/zhijian_lib.dir/src/parallel/mpi_comm.cpp.o
[ 66%] Building CXX object CMakeFiles/zhijian_lib.dir/src/io/restart.cpp.o
[ 66%] Building CXX object CMakeFiles/zhijian_lib.dir/src/io/visualization.cpp.o
[ 66%] Building CXX object CMakeFiles/zhijian_lib.dir/src/basis/basis.cpp.o
[ 66%] Building CXX object CMakeFiles/zhijian_lib.dir/src/solver/fr_solver.cpp.o
[ 73%] Linking C static library libhdf5_compat.a
[ 73%] Built target hdf5_compat
In file included from /home/z651w035/codes/Zhijian/include/solver/fr_solver.hpp:4,
                 from /home/z651w035/codes/Zhijian/include/io/visualization.hpp:5,
                 from /home/z651w035/codes/Zhijian/src/io/visualization.cpp:1:
/home/z651w035/codes/Zhijian/include/common/cuda_utils.hpp: In function ‘double zhijian::atomicAdd_double(double*, double)’:
/home/z651w035/codes/Zhijian/include/common/cuda_utils.hpp:18:52: error: ‘__longlong_as_double’ was not declared in this scope
   18 |                         __double_as_longlong(val + __longlong_as_double(assumed)));
      |                                                    ^~~~~~~~~~~~~~~~~~~~
/home/z651w035/codes/Zhijian/include/common/cuda_utils.hpp:18:25: error: ‘__double_as_longlong’ was not declared in this scope
   18 |                         __double_as_longlong(val + __longlong_as_double(assumed)));
      |                         ^~~~~~~~~~~~~~~~~~~~
/home/z651w035/codes/Zhijian/include/common/cuda_utils.hpp:17:15: error: ‘atomicCAS’ was not declared in this scope
   17 |         old = atomicCAS(address_as_ull, assumed,
      |               ^~~~~~~~~
/home/z651w035/codes/Zhijian/include/common/cuda_utils.hpp:20:12: error: ‘__longlong_as_double’ was not declared in this scope
   20 |     return __longlong_as_double(old);
      |            ^~~~~~~~~~~~~~~~~~~~
In file included from /home/z651w035/codes/Zhijian/include/solver/fr_solver.hpp:4,
                 from /home/z651w035/codes/Zhijian/src/solver/fr_solver.cpp:1:
/home/z651w035/codes/Zhijian/include/common/cuda_utils.hpp: In function ‘double zhijian::atomicAdd_double(double*, double)’:
/home/z651w035/codes/Zhijian/include/common/cuda_utils.hpp:18:52: error: ‘__longlong_as_double’ was not declared in this scope
   18 |                         __double_as_longlong(val + __longlong_as_double(assumed)));
      |                                                    ^~~~~~~~~~~~~~~~~~~~
/home/z651w035/codes/Zhijian/include/common/cuda_utils.hpp:18:25: error: ‘__double_as_longlong’ was not declared in this scope
   18 |                         __double_as_longlong(val + __longlong_as_double(assumed)));
      |                         ^~~~~~~~~~~~~~~~~~~~
/home/z651w035/codes/Zhijian/include/common/cuda_utils.hpp:17:15: error: ‘atomicCAS’ was not declared in this scope
   17 |         old = atomicCAS(address_as_ull, assumed,
      |               ^~~~~~~~~
/home/z651w035/codes/Zhijian/include/common/cuda_utils.hpp:20:12: error: ‘__longlong_as_double’ was not declared in this scope
   20 |     return __longlong_as_double(old);
      |            ^~~~~~~~~~~~~~~~~~~~
make[2]: *** [CMakeFiles/zhijian_lib.dir/build.make:160: CMakeFiles/zhijian_lib.dir/src/io/visualization.cpp.o] Error 1
make[2]: *** Waiting for unfinished jobs....
In file included from /home/z651w035/codes/Zhijian/include/solver/fr_solver.hpp:4,
                 from /home/z651w035/codes/Zhijian/include/io/restart.hpp:5,
                 from /home/z651w035/codes/Zhijian/src/io/restart.cpp:1:
/home/z651w035/codes/Zhijian/include/common/cuda_utils.hpp: In function ‘double zhijian::atomicAdd_double(double*, double)’:
/home/z651w035/codes/Zhijian/include/common/cuda_utils.hpp:18:52: error: ‘__longlong_as_double’ was not declared in this scope
   18 |                         __double_as_longlong(val + __longlong_as_double(assumed)));
      |                                                    ^~~~~~~~~~~~~~~~~~~~
/home/z651w035/codes/Zhijian/include/common/cuda_utils.hpp:18:25: error: ‘__double_as_longlong’ was not declared in this scope
   18 |                         __double_as_longlong(val + __longlong_as_double(assumed)));
      |                         ^~~~~~~~~~~~~~~~~~~~
/home/z651w035/codes/Zhijian/include/common/cuda_utils.hpp:17:15: error: ‘atomicCAS’ was not declared in this scope
   17 |         old = atomicCAS(address_as_ull, assumed,
      |               ^~~~~~~~~
/home/z651w035/codes/Zhijian/include/common/cuda_utils.hpp:20:12: error: ‘__longlong_as_double’ was not declared in this scope
   20 |     return __longlong_as_double(old);
      |            ^~~~~~~~~~~~~~~~~~~~
make[2]: *** [CMakeFiles/zhijian_lib.dir/build.make:132: CMakeFiles/zhijian_lib.dir/src/solver/fr_solver.cpp.o] Error 1
make[2]: *** [CMakeFiles/zhijian_lib.dir/build.make:146: CMakeFiles/zhijian_lib.dir/src/io/restart.cpp.o] Error 1
/home/z651w035/codes/Zhijian/src/gpu/kernels.cu(125): warning #177-D: variable "var" was declared but never referenced
      int var = idx % N_VARS;
          ^

Remark: The warnings can be suppressed with "-diag-suppress <warning-number>"

make[1]: *** [CMakeFiles/Makefile2:87: CMakeFiles/zhijian_lib.dir/all] Error 2
make: *** [Makefile:136: all] Error 2

