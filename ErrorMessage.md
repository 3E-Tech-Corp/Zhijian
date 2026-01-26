--
-- === Zhijian Configuration Summary ===
-- Build type: Release
-- C++ compiler: /usr/bin/g++
-- CUDA compiler: /kuhpc/sw/nvhpc/Linux_x86_64/25.3/cuda/12.8/bin/nvcc
-- CUDA architectures: 70;75;80;86;89;90
-- MPI: /kuhpc/sw/openmpi/5.0.5/gcc/11.4/lib/libmpi.so
-- HDF5: /kuhpc/sw/hdf5/1.14.5/gcc/11.4/lib/libhdf5.so;/usr/lib64/libz.so;/usr/lib64/libdl.a;/usr/lib64/libm.so;/kuhpc/sw/hdf5/1.14.5/gcc/11.4/lib/libhdf5_cpp.so;/kuhpc/sw/hdf5/1.14.5/gcc/11.4/lib/libhdf5.so;/usr/lib64/libz.so;/usr/lib64/libdl.a;/usr/lib64/libm.so
-- CGNS: /home/z651w035/hpmusic/contrib/cgnslib/lib/libcgns.a
-- METIS: /home/z651w035/local/metis/lib
-- =====================================
--
-- Configuring done (4.3s)
CMake Warning at CMakeLists.txt:123 (target_link_libraries):
  Target "zhijian" requests linking to directory
  "/home/z651w035/local/metis/lib".  Targets may link only to libraries.
  CMake is dropping the item.


CMake Warning at CMakeLists.txt:123 (target_link_libraries):
  Target "zhijian" requests linking to directory
  "/home/z651w035/local/metis/lib".  Targets may link only to libraries.
  CMake is dropping the item.


-- Generating done (3.4s)
-- Build files have been written to: /home/z651w035/codes/Zhijian/release
[ 76%] Built target zhijian_lib
[ 84%] Linking CUDA device code CMakeFiles/zhijian.dir/cmake_device_link.o
nvlink warning : Skipping incompatible '/usr/lib64/libdl.a' when searching for -ldl (target: sm_70)
nvlink warning : Skipping incompatible '/usr/lib64/librt.a' when searching for -lrt (target: sm_70)
nvlink warning : Skipping incompatible '/usr/lib64/libpthread.a' when searching for -lpthread (target: sm_70)
nvlink warning : Skipping incompatible '/usr/lib64/libdl.a' when searching for -ldl (target: sm_75)
nvlink warning : Skipping incompatible '/usr/lib64/librt.a' when searching for -lrt (target: sm_75)
nvlink warning : Skipping incompatible '/usr/lib64/libpthread.a' when searching for -lpthread (target: sm_75)
nvlink warning : Skipping incompatible '/usr/lib64/libdl.a' when searching for -ldl (target: sm_80)
nvlink warning : Skipping incompatible '/usr/lib64/librt.a' when searching for -lrt (target: sm_80)
nvlink warning : Skipping incompatible '/usr/lib64/libpthread.a' when searching for -lpthread (target: sm_80)
nvlink warning : Skipping incompatible '/usr/lib64/libdl.a' when searching for -ldl (target: sm_86)
nvlink warning : Skipping incompatible '/usr/lib64/librt.a' when searching for -lrt (target: sm_86)
nvlink warning : Skipping incompatible '/usr/lib64/libpthread.a' when searching for -lpthread (target: sm_86)
nvlink warning : Skipping incompatible '/usr/lib64/libdl.a' when searching for -ldl (target: sm_89)
nvlink warning : Skipping incompatible '/usr/lib64/librt.a' when searching for -lrt (target: sm_89)
nvlink warning : Skipping incompatible '/usr/lib64/libpthread.a' when searching for -lpthread (target: sm_89)
nvlink warning : Skipping incompatible '/usr/lib64/libdl.a' when searching for -ldl (target: sm_90)
nvlink warning : Skipping incompatible '/usr/lib64/librt.a' when searching for -lrt (target: sm_90)
nvlink warning : Skipping incompatible '/usr/lib64/libpthread.a' when searching for -lpthread (target: sm_90)
[ 92%] Linking CXX executable zhijian
/usr/bin/ld: libzhijian_lib.a(partitioner.cpp.o): in function `zhijian::MeshPartitioner::partition(zhijian::Mesh&, int)':
partitioner.cpp:(.text+0x19c7): undefined reference to `METIS_SetDefaultOptions'
/usr/bin/ld: partitioner.cpp:(.text+0x1a20): undefined reference to `METIS_PartGraphKway'
/usr/bin/ld: /home/z651w035/hpmusic/contrib/cgnslib/lib/libcgns.a(ADFH.c.o): in function `ADFH_Children_IDs':
ADFH.c:(.text+0x441f): undefined reference to `H5Literate'
/usr/bin/ld: ADFH.c:(.text+0x453b): undefined reference to `H5Literate'
/usr/bin/ld: /home/z651w035/hpmusic/contrib/cgnslib/lib/libcgns.a(ADFH.c.o): in function `ADFH_Children_Names':
ADFH.c:(.text+0x581f): undefined reference to `H5Literate'
/usr/bin/ld: ADFH.c:(.text+0x592b): undefined reference to `H5Literate'
collect2: error: ld returned 1 exit status
make[2]: *** [CMakeFiles/zhijian.dir/build.make:153: zhijian] Error 1
make[1]: *** [CMakeFiles/Makefile2:111: CMakeFiles/zhijian.dir/all] Error 2
make: *** [Makefile:136: all] Error 2

