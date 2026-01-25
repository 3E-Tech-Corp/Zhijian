# Zhijian - 2D High-Order Flow Solver

A GPU-accelerated 2D compressible Navier-Stokes solver using the Flux Reconstruction (FR) method.

## Features

- **High-Order Accuracy**: Flux Reconstruction method with polynomial orders P1 to P5
- **GPU Acceleration**: CUDA kernels for all compute-intensive operations
- **MPI Parallelization**: Distributed computing across multiple nodes and GPUs
- **Mesh Support**: CGNS format with Q1/Q2 triangular and quadrilateral elements
- **Time Integration**: 3rd-order Strong Stability Preserving Runge-Kutta (SSP-RK3)
- **Boundary Conditions**: Wall, slip wall, symmetry, far-field, inflow/outflow
- **Riemann Solvers**: Rusanov (LLF), Roe, and HLLC
- **I/O**: HDF5 restart files, VTK visualization output

## Requirements

- CMake >= 3.20
- CUDA Toolkit >= 11.0
- MPI implementation (OpenMPI, MPICH, etc.)
- HDF5 with C++ bindings
- CGNS library
- METIS for mesh partitioning
- CUB library (included with CUDA Toolkit)

## Building

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release \
         -DCGNS_DIR=/path/to/cgns \
         -DMETIS_DIR=/path/to/metis
make -j
```

## Usage

```bash
# Single GPU
./zhijian -m mesh.cgns -c config.txt -o output

# Multiple GPUs with MPI
mpirun -np 4 ./zhijian -m mesh.cgns -c config.txt -o output

# Restart from checkpoint
./zhijian -m mesh.cgns -c config.txt -r restart.h5 -o output
```

## Configuration File

See `examples/config.txt` for a complete example. Key parameters:

```
# Physical parameters
gamma       1.4         # Ratio of specific heats
Re          5000        # Reynolds number
Mach        0.5         # Free-stream Mach number
AoA         2.0         # Angle of attack (degrees)

# Numerical parameters
poly_order  3           # Polynomial order (1-5)
flux_type   HU          # FR flux: DG, SD, HU, GA
riemann     Roe         # Riemann solver: Rusanov, Roe, HLLC
CFL         0.5         # CFL number
```

## Flux Reconstruction Schemes

The solver supports several FR correction functions:

- **DG**: Standard Discontinuous Galerkin (c=0)
- **SD**: Spectral Difference
- **HU**: Huynh's g2 scheme (energy stable)
- **GA**: Gauss scheme

## Governing Equations

The solver solves the 2D compressible Navier-Stokes equations:

$$\frac{\partial \mathbf{U}}{\partial t} + \nabla \cdot (\mathbf{F}_c - \mathbf{F}_v) = 0$$

where:
- **U** = [ρ, ρu, ρv, ρE]ᵀ (conservative variables)
- **F_c** = convective (Euler) fluxes
- **F_v** = viscous fluxes

With the ideal gas equation of state: p = (γ-1)ρe

## Output Formats

- **VTK/VTU**: Visualization files for ParaView, VisIt
- **HDF5**: Restart files for checkpoint/restart capability
- **CSV**: Residual history

## Directory Structure

```
Zhijian/
├── include/           # Header files
│   ├── common/        # Common types and utilities
│   ├── mesh/          # Mesh data structures
│   ├── basis/         # FR basis functions
│   ├── physics/       # Euler and NS equations
│   ├── solver/        # FR solver
│   ├── time/          # Time integration
│   ├── bc/            # Boundary conditions
│   ├── io/            # I/O operations
│   ├── parallel/      # MPI communication
│   └── gpu/           # CUDA kernels
├── src/               # Source files
├── examples/          # Example configurations
└── tests/             # Unit tests
```

## References

1. Huynh, H.T. "A flux reconstruction approach to high-order schemes including discontinuous Galerkin methods." AIAA 2007.
2. Vincent, P.E., et al. "A new class of high-order energy stable flux reconstruction schemes." J. Sci. Comput. 2011.
3. Castonguay, P., et al. "Energy stable flux reconstruction schemes for advection-diffusion problems." Comput. Methods Appl. Mech. Eng. 2013.

## License

MIT License



nvcc warning : Support for offline compilation for architectures prior to '<compute/sm/lto>_75' will be removed in a future release (Use -Wno-deprecated-gpu-targets to suppress warning).
In file included from /home/z651w035/codes/Zhijian/include/io/visualization.hpp:5,
                 from /home/z651w035/codes/Zhijian/src/io/visualization.cpp:1:
/home/z651w035/codes/Zhijian/include/solver/fr_solver.hpp:73:41: error: 'function' in namespace 'std' does not name a template type
   73 |     void setInitialCondition(const std::function<State(Vec2)>& ic_func);
      |                                         ^~~~~~~~
/home/z651w035/codes/Zhijian/include/solver/fr_solver.hpp:9:1: note: 'std::function' is defined in header '<functional>'; did you forget to '#include <functional>'?
    8 | #include "physics/navierstokes.hpp"
  +++ |+#include <functional>
    9 | #include <memory>
/home/z651w035/codes/Zhijian/include/solver/fr_solver.hpp:73:49: error: expected ',' or '...' before '<' token
   73 |     void setInitialCondition(const std::function<State(Vec2)>& ic_func);
      |                                                 ^
In file included from /home/z651w035/codes/Zhijian/src/solver/fr_solver.cpp:1:
/home/z651w035/codes/Zhijian/include/solver/fr_solver.hpp:73:41: error: 'function' in namespace 'std' does not name a template type
   73 |     void setInitialCondition(const std::function<State(Vec2)>& ic_func);
      |                                         ^~~~~~~~
/home/z651w035/codes/Zhijian/include/solver/fr_solver.hpp:9:1: note: 'std::function' is defined in header '<functional>'; did you forget to '#include <functional>'?
    8 | #include "physics/navierstokes.hpp"
  +++ |+#include <functional>
    9 | #include <memory>
/home/z651w035/codes/Zhijian/include/solver/fr_solver.hpp:73:49: error: expected ',' or '...' before '<' token
   73 |     void setInitialCondition(const std::function<State(Vec2)>& ic_func);
      |                                                 ^
/home/z651w035/codes/Zhijian/src/solver/fr_solver.cpp:158:47: error: 'function' in namespace 'std' does not name a template type
  158 | void FRSolver::setInitialCondition(const std::function<State(Vec2)>& ic_func) {
      |                                               ^~~~~~~~
/home/z651w035/codes/Zhijian/src/solver/fr_solver.cpp:7:1: note: 'std::function' is defined in header '<functional>'; did you forget to '#include <functional>'?
    6 | #include <cmath>
  +++ |+#include <functional>
    7 | 
/home/z651w035/codes/Zhijian/src/solver/fr_solver.cpp:158:55: error: expected ',' or '...' before '<' token
  158 | void FRSolver::setInitialCondition(const std::function<State(Vec2)>& ic_func) {
      |                                                       ^
/home/z651w035/codes/Zhijian/src/solver/fr_solver.cpp: In member function 'void zhijian::FRSolver::setInitialCondition(int)':
/home/z651w035/codes/Zhijian/src/solver/fr_solver.cpp:180:40: error: 'ic_func' was not declared in this scope
  180 |             solution_[i].sol_pts[sp] = ic_func(pos);
      |                                        ^~~~~~~
In file included from /home/z651w035/codes/Zhijian/include/io/restart.hpp:5,
                 from /home/z651w035/codes/Zhijian/src/io/restart.cpp:1:
/home/z651w035/codes/Zhijian/include/solver/fr_solver.hpp:73:41: error: 'function' in namespace 'std' does not name a template type
   73 |     void setInitialCondition(const std::function<State(Vec2)>& ic_func);
      |                                         ^~~~~~~~
/home/z651w035/codes/Zhijian/include/solver/fr_solver.hpp:9:1: note: 'std::function' is defined in header '<functional>'; did you forget to '#include <functional>'?
    8 | #include "physics/navierstokes.hpp"
  +++ |+#include <functional>
    9 | #include <memory>
/home/z651w035/codes/Zhijian/include/solver/fr_solver.hpp:73:49: error: expected ',' or '...' before '<' token
   73 |     void setInitialCondition(const std::function<State(Vec2)>& ic_func);
      |                                                 ^
make[2]: *** [CMakeFiles/zhijian_lib.dir/build.make:160: CMakeFiles/zhijian_lib.dir/src/io/visualization.cpp.o] Error 1
make[2]: *** Waiting for unfinished jobs....
/home/z651w035/codes/Zhijian/src/io/restart.cpp: In member function 'bool zhijian::RestartIO::writeParallel(const string&, int, int, const zhijian::Mesh&, const std::vector<zhijian::ElementSolution>&, const zhijian::SimParams&, zhijian::Real, int)':
/home/z651w035/codes/Zhijian/src/io/restart.cpp:175:32: error: 'MPI_COMM_WORLD' was not declared in this scope
  175 |     H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
      |                                ^~~~~~~~~~~~~~
/home/z651w035/codes/Zhijian/src/io/restart.cpp:175:48: error: 'MPI_INFO_NULL' was not declared in this scope
  175 |     H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
      |                                                ^~~~~~~~~~~~~
/home/z651w035/codes/Zhijian/src/io/restart.cpp:175:5: error: 'H5Pset_fapl_mpio' was not declared in this scope; did you mean 'H5Pset_fapl_stdio'?
  175 |     H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
      |     ^~~~~~~~~~~~~~~~
      |     H5Pset_fapl_stdio
/home/z651w035/codes/Zhijian/src/io/restart.cpp: In member function 'bool zhijian::RestartIO::readParallel(const string&, int, int, zhijian::Mesh&, std::vector<zhijian::ElementSolution>&, zhijian::SimParams&, zhijian::Real&, int&)':
/home/z651w035/codes/Zhijian/src/io/restart.cpp:245:32: error: 'MPI_COMM_WORLD' was not declared in this scope
  245 |     H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
      |                                ^~~~~~~~~~~~~~
make[2]: *** [CMakeFiles/zhijian_lib.dir/build.make:132: CMakeFiles/zhijian_lib.dir/src/solver/fr_solver.cpp.o] Error 1
/home/z651w035/codes/Zhijian/src/io/restart.cpp:245:48: error: 'MPI_INFO_NULL' was not declared in this scope
  245 |     H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
      |                                                ^~~~~~~~~~~~~
/home/z651w035/codes/Zhijian/src/io/restart.cpp:245:5: error: 'H5Pset_fapl_mpio' was not declared in this scope; did you mean 'H5Pset_fapl_stdio'?
  245 |     H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
      |     ^~~~~~~~~~~~~~~~
      |     H5Pset_fapl_stdio
make[2]: *** [CMakeFiles/zhijian_lib.dir/build.make:146: CMakeFiles/zhijian_lib.dir/src/io/restart.cpp.o] Error 1
In file included from /kuhpc/sw/rocm/6.3.1/include/thrust/system/cuda/detail/execution_policy.h:35,
                 from /kuhpc/sw/rocm/6.3.1/include/thrust/iterator/detail/device_system_tag.h:23,
                 from /kuhpc/sw/rocm/6.3.1/include/thrust/iterator/detail/iterator_facade_category.h:22,
                 from /kuhpc/sw/rocm/6.3.1/include/thrust/iterator/iterator_facade.h:37,
                 from /kuhpc/sw/nvhpc/Linux_x86_64/25.3/cuda/12.8/targets/x86_64-linux/include/cub/iterator/cache_modified_input_iterator.cuh:58,
                 from /kuhpc/sw/nvhpc/Linux_x86_64/25.3/cuda/12.8/targets/x86_64-linux/include/cub/block/block_load.cuh:45,
                 from /kuhpc/sw/nvhpc/Linux_x86_64/25.3/cuda/12.8/targets/x86_64-linux/include/cub/cub.cuh:52,
                 from /home/z651w035/codes/Zhijian/src/gpu/kernels.cu:4:
/kuhpc/sw/rocm/6.3.1/include/thrust/system/cuda/config.h:116:2: error: #error The version of CUB in your include path is not compatible with this release of Thrust. CUB is now included in the CUDA Toolkit, so you no longer need to use your own checkout of CUB. Define THRUST_IGNORE_CUB_VERSION_CHECK to ignore this.
  116 | #error The version of CUB in your include path is not compatible with this release of Thrust. CUB is now included in the CUDA Toolkit, so you no longer need to use your own checkout of CUB. Define THRUST_IGNORE_CUB_VERSION_CHECK to ignore this.
      |  ^~~~~
make[2]: *** [CMakeFiles/zhijian_lib.dir/build.make:189: CMakeFiles/zhijian_lib.dir/src/gpu/kernels.cu.o] Error 1
make[1]: *** [CMakeFiles/Makefile2:85: CMakeFiles/zhijian_lib.dir/all] Error 2
make: *** [Makefile:136: all] Error 2

