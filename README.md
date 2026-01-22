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
