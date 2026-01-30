cmake_minimum_required(VERSION 3.20)
project(Zhijian LANGUAGES C CXX CUDA)

# C++ and CUDA standards
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# CUDA architecture (adjust based on target GPU)
# Default to modern architectures; override with -DCMAKE_CUDA_ARCHITECTURES="80;86;89;90"
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 75 80 86 89 90)
endif()

# Suppress deprecation warnings for older architectures
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Wno-deprecated-gpu-targets")

# Find required packages
find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(CUDAToolkit REQUIRED)

# Explicitly find HDF5 C library to ensure it's linked (needed for H5Literate)
# Some systems don't include the C library in HDF5_LIBRARIES properly
get_filename_component(_hdf5_lib_dir "${HDF5_C_LIBRARIES}" DIRECTORY)
if(NOT _hdf5_lib_dir OR _hdf5_lib_dir STREQUAL "")
    # Fallback: extract directory from HDF5_LIBRARIES
    list(GET HDF5_LIBRARIES 0 _hdf5_first_lib)
    get_filename_component(_hdf5_lib_dir "${_hdf5_first_lib}" DIRECTORY)
endif()

find_library(HDF5_C_LIB
    NAMES hdf5 hdf5-shared
    HINTS ${_hdf5_lib_dir} ${HDF5_ROOT}/lib ${HDF5_DIR}/lib
    PATH_SUFFIXES lib lib64)

if(HDF5_C_LIB)
    message(STATUS "Found HDF5 C library: ${HDF5_C_LIB}")
else()
    message(STATUS "Using HDF5_LIBRARIES: ${HDF5_LIBRARIES}")
endif()

# Find CGNS
find_path(CGNS_INCLUDE_DIR cgnslib.h
    HINTS ${CGNS_DIR}/include /usr/local/include /usr/include)
find_library(CGNS_LIBRARY cgns
    HINTS ${CGNS_DIR}/lib /usr/local/lib /usr/lib)

if(NOT CGNS_INCLUDE_DIR OR NOT CGNS_LIBRARY)
    message(FATAL_ERROR "CGNS library not found. Set CGNS_DIR.")
endif()

# Find METIS for mesh partitioning
find_path(METIS_INCLUDE_DIR metis.h
    HINTS ${METIS_DIR}/include /usr/local/include /usr/include)
find_library(METIS_LIBRARY metis
    HINTS ${METIS_DIR}/lib /usr/local/lib /usr/lib)

# Validate METIS_LIBRARY - if user passed a directory instead of library file, search inside it
if(METIS_LIBRARY AND IS_DIRECTORY "${METIS_LIBRARY}")
    message(STATUS "METIS_LIBRARY points to directory, searching for library inside...")
    set(_metis_search_dir "${METIS_LIBRARY}")
    unset(METIS_LIBRARY CACHE)
    find_library(METIS_LIBRARY NAMES metis
        PATHS "${_metis_search_dir}"
        NO_DEFAULT_PATH)
    if(NOT METIS_LIBRARY)
        # Try again with default paths
        find_library(METIS_LIBRARY metis
            HINTS ${METIS_DIR}/lib /usr/local/lib /usr/lib)
    endif()
endif()

# Find GKlib (METIS dependency) - try multiple possible library names
# GKlib is usually installed alongside METIS
get_filename_component(_metis_lib_dir "${METIS_LIBRARY}" DIRECTORY)
find_library(GKLIB_LIBRARY
    NAMES GKlib gklib gk
    HINTS ${_metis_lib_dir} ${METIS_DIR}/lib ${GKLIB_DIR}/lib /usr/local/lib /usr/lib)

if(NOT METIS_INCLUDE_DIR OR NOT METIS_LIBRARY)
    message(FATAL_ERROR "METIS library not found. Set METIS_DIR to the METIS installation directory.")
endif()

# Verify METIS_LIBRARY is a file
if(NOT EXISTS "${METIS_LIBRARY}" OR IS_DIRECTORY "${METIS_LIBRARY}")
    message(FATAL_ERROR "METIS_LIBRARY='${METIS_LIBRARY}' is not a valid library file. "
        "Please set METIS_DIR to the installation prefix (e.g., /usr/local) or "
        "METIS_LIBRARY to the full path of libmetis.a or libmetis.so")
endif()

if(GKLIB_LIBRARY)
    message(STATUS "Found GKlib: ${GKLIB_LIBRARY}")
else()
    message(FATAL_ERROR "GKlib library not found. METIS requires GKlib for linking.\n"
        "Please set -DGKLIB_LIBRARY=/path/to/libGKlib.a when running cmake.\n"
        "GKlib is typically built alongside METIS. Check:\n"
        "  - ${_metis_lib_dir}/libGKlib.a\n"
        "  - The METIS build directory (e.g., build/libGKlib.a)\n"
        "  - Or rebuild METIS with: make config shared=1 prefix=/your/prefix && make install")
endif()

# Include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${MPI_CXX_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
    ${CGNS_INCLUDE_DIR}
    ${METIS_INCLUDE_DIR}
    ${CUDAToolkit_INCLUDE_DIRS}
)

# Source files
file(GLOB_RECURSE COMMON_SOURCES "src/common/*.cpp")
file(GLOB_RECURSE MESH_SOURCES "src/mesh/*.cpp")
file(GLOB_RECURSE BASIS_SOURCES "src/basis/*.cpp")
file(GLOB_RECURSE PHYSICS_SOURCES "src/physics/*.cpp")
file(GLOB_RECURSE SOLVER_SOURCES "src/solver/*.cpp")
file(GLOB_RECURSE TIME_SOURCES "src/time/*.cpp")
file(GLOB_RECURSE BC_SOURCES "src/bc/*.cpp")
file(GLOB_RECURSE IO_SOURCES "src/io/*.cpp")
file(GLOB_RECURSE PARALLEL_SOURCES "src/parallel/*.cpp")

# CUDA source files
file(GLOB_RECURSE CUDA_SOURCES "src/gpu/*.cu")

# Create library
add_library(zhijian_lib STATIC
    ${COMMON_SOURCES}
    ${MESH_SOURCES}
    ${BASIS_SOURCES}
    ${PHYSICS_SOURCES}
    ${SOLVER_SOURCES}
    ${TIME_SOURCES}
    ${BC_SOURCES}
    ${IO_SOURCES}
    ${PARALLEL_SOURCES}
    ${CUDA_SOURCES}
)

# Set CUDA properties
set_target_properties(zhijian_lib PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    POSITION_INDEPENDENT_CODE ON
)

# Suppress CUB/Thrust version conflict when ROCm headers are in the include path
target_compile_definitions(zhijian_lib PRIVATE THRUST_IGNORE_CUB_VERSION_CHECK)

# Link libraries to zhijian_lib — use MPI::MPI_CXX target for proper include propagation
target_link_libraries(zhijian_lib
    CUDA::cudart
    CUDA::cublas
    MPI::MPI_CXX
)

# HDF5 compatibility shim for older CGNS libraries
# Provides H5Literate symbol for CGNS built against HDF5 < 1.12
add_library(hdf5_compat STATIC src/hdf5_compat.c)
target_include_directories(hdf5_compat PRIVATE ${HDF5_INCLUDE_DIRS})

# Main executable
add_executable(zhijian src/main.cpp)

# For static library linking with circular dependencies, use linker groups
# This tells the linker to iterate through libraries until all symbols are resolved
# Format: -Wl,--start-group <libs> -Wl,--end-group

# Build HDF5 link list - ensure C library is included
set(HDF5_LINK_LIBS ${HDF5_LIBRARIES})
if(HDF5_C_LIB)
    # Prepend explicit C library to ensure H5Literate is found
    list(INSERT HDF5_LINK_LIBS 0 ${HDF5_C_LIB})
endif()

target_link_libraries(zhijian
    zhijian_lib
    ${MPI_CXX_LIBRARIES}
    # Use linker group for static libraries with interdependencies
    -Wl,--start-group
    ${METIS_LIBRARY}
    ${GKLIB_LIBRARY}
    ${CGNS_LIBRARY}
    hdf5_compat
    ${HDF5_LINK_LIBS}
    -Wl,--end-group
    m
)

# Install
install(TARGETS zhijian DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)

# Print configuration summary
message(STATUS "")
message(STATUS "=== Zhijian Configuration Summary ===")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "CUDA compiler: ${CMAKE_CUDA_COMPILER}")
message(STATUS "CUDA architectures: ${CMAKE_CUDA_ARCHITECTURES}")
message(STATUS "MPI: ${MPI_CXX_LIBRARIES}")
message(STATUS "HDF5: ${HDF5_LINK_LIBS}")
message(STATUS "CGNS: ${CGNS_LIBRARY}")
message(STATUS "METIS: ${METIS_LIBRARY}")
message(STATUS "GKlib: ${GKLIB_LIBRARY}")
message(STATUS "=====================================")
message(STATUS "")
