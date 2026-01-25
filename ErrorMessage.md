[z651w035@login1 release]$ make -j
[  7%] Building CXX object CMakeFiles/zhijian_lib.dir/src/solver/fr_solver.cpp.o
[ 23%] Building CXX object CMakeFiles/zhijian_lib.dir/src/io/visualization.cpp.o
[ 23%] Building CUDA object CMakeFiles/zhijian_lib.dir/src/gpu/kernels.cu.o
nvcc warning : Support for offline compilation for architectures prior to '<compute/sm/lto>_75' will be removed in a future release (Use -Wno-deprecated-gpu-targets to suppress warning).
[ 30%] Building CXX object CMakeFiles/zhijian_lib.dir/src/io/restart.cpp.o
/home/z651w035/codes/Zhijian/src/io/restart.cpp: In member function ‘bool zhijian::RestartIO::writeParallel(const string&, int, int, const zhijian::Mesh&, const std::vector<zhijian::ElementSolution>&, const zhijian::SimParams&, zhijian::Real, int)’:
/home/z651w035/codes/Zhijian/src/io/restart.cpp:176:5: error: ‘H5Pset_fapl_mpio’ was not declared in this scope; did you mean ‘H5Pset_fapl_stdio’?
  176 |     H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
      |     ^~~~~~~~~~~~~~~~
      |     H5Pset_fapl_stdio
/home/z651w035/codes/Zhijian/src/io/restart.cpp: In member function ‘bool zhijian::RestartIO::readParallel(const string&, int, int, zhijian::Mesh&, std::vector<zhijian::ElementSolution>&, zhijian::SimParams&, zhijian::Real&, int&)’:
/home/z651w035/codes/Zhijian/src/io/restart.cpp:246:5: error: ‘H5Pset_fapl_mpio’ was not declared in this scope; did you mean ‘H5Pset_fapl_stdio’?
  246 |     H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
      |     ^~~~~~~~~~~~~~~~
      |     H5Pset_fapl_stdio
make[2]: *** [CMakeFiles/zhijian_lib.dir/build.make:146: CMakeFiles/zhijian_lib.dir/src/io/restart.cpp.o] Error 1
make[2]: *** Waiting for unfinished jobs....
/home/z651w035/codes/Zhijian/src/gpu/kernels.cu(740): error: identifier "DeviceArray" is undefined
      DeviceArray<Real> dt_local(n_elem);
      ^

/home/z651w035/codes/Zhijian/src/gpu/kernels.cu(740): error: type name is not allowed
      DeviceArray<Real> dt_local(n_elem);
                  ^

/home/z651w035/codes/Zhijian/src/gpu/kernels.cu(740): error: identifier "dt_local" is undefined
      DeviceArray<Real> dt_local(n_elem);
                        ^

/home/z651w035/codes/Zhijian/src/gpu/kernels.cu(754): error: type name is not allowed
      DeviceArray<Real> partial_min(reduce_blocks);
                  ^

/home/z651w035/codes/Zhijian/src/gpu/kernels.cu(754): error: identifier "partial_min" is undefined
      DeviceArray<Real> partial_min(reduce_blocks);
                        ^

/home/z651w035/codes/Zhijian/src/gpu/kernels.cu(761): error: type name is not allowed
          DeviceArray<Real> temp_min(new_blocks);
                      ^

/home/z651w035/codes/Zhijian/src/gpu/kernels.cu(761): error: identifier "temp_min" is undefined
          DeviceArray<Real> temp_min(new_blocks);

