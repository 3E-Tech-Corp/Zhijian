cmake_minimum_required(VERSION 3.20)
project(Zhijian LANGUAGES CXX CUDA)

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
set(CMAKE_CUDA_ARCHITECTURES 70 75 80 86 89 90)

# Find required packages
find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(CUDAToolkit REQUIRED)

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

if(NOT METIS_INCLUDE_DIR OR NOT METIS_LIBRARY)
    message(FATAL_ERROR "METIS library not found. Set METIS_DIR.")
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

# Link libraries to zhijian_lib
target_link_libraries(zhijian_lib
    ${MPI_CXX_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${CGNS_LIBRARY}
    ${METIS_LIBRARY}
    CUDA::cudart
    CUDA::cublas
)

# Main executable
add_executable(zhijian src/main.cpp)
target_link_libraries(zhijian zhijian_lib)

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
message(STATUS "HDF5: ${HDF5_LIBRARIES}")
message(STATUS "CGNS: ${CGNS_LIBRARY}")
message(STATUS "METIS: ${METIS_LIBRARY}")
message(STATUS "=====================================")
message(STATUS "")
