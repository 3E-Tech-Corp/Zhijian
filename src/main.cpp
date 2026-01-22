#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <memory>

#include "common/types.hpp"
#include "mesh/mesh.hpp"
#include "solver/fr_solver.hpp"
#include "bc/boundary_condition.hpp"
#include "parallel/mpi_comm.hpp"
#include "io/restart.hpp"
#include "io/visualization.hpp"

using namespace zhijian;

// Parse command line arguments
struct ProgramOptions {
    std::string mesh_file;
    std::string config_file;
    std::string restart_file;
    std::string output_dir = "./output";
    int num_threads = 1;
    bool help = false;
};

void printUsage(const char* prog_name) {
    std::cout << "Usage: " << prog_name << " [options]\n"
              << "\n"
              << "Options:\n"
              << "  -m, --mesh <file>       Input mesh file (CGNS format)\n"
              << "  -c, --config <file>     Configuration file\n"
              << "  -r, --restart <file>    Restart from file\n"
              << "  -o, --output <dir>      Output directory (default: ./output)\n"
              << "  -h, --help              Show this help message\n"
              << "\n"
              << "Example:\n"
              << "  mpirun -np 4 " << prog_name << " -m naca0012.cgns -c config.txt\n";
}

ProgramOptions parseArgs(int argc, char* argv[]) {
    ProgramOptions opts;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            opts.help = true;
        } else if ((arg == "-m" || arg == "--mesh") && i + 1 < argc) {
            opts.mesh_file = argv[++i];
        } else if ((arg == "-c" || arg == "--config") && i + 1 < argc) {
            opts.config_file = argv[++i];
        } else if ((arg == "-r" || arg == "--restart") && i + 1 < argc) {
            opts.restart_file = argv[++i];
        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            opts.output_dir = argv[++i];
        }
    }

    return opts;
}

// Parse configuration file
SimParams parseConfig(const std::string& filename) {
    SimParams params;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Warning: Could not open config file: " << filename << std::endl;
        return params;
    }

    std::string line;
    while (std::getline(file, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string key;
        iss >> key;

        if (key == "gamma") iss >> params.gamma;
        else if (key == "Pr") iss >> params.Pr;
        else if (key == "Re") iss >> params.Re;
        else if (key == "Mach") iss >> params.Mach_inf;
        else if (key == "AoA") iss >> params.AoA;
        else if (key == "T_inf") iss >> params.T_inf;
        else if (key == "p_inf") iss >> params.p_inf;
        else if (key == "rho_inf") iss >> params.rho_inf;
        else if (key == "poly_order") iss >> params.poly_order;
        else if (key == "CFL") iss >> params.CFL;
        else if (key == "dt") iss >> params.dt;
        else if (key == "t_final") iss >> params.t_final;
        else if (key == "max_iter") iss >> params.max_iter;
        else if (key == "viscous") {
            std::string val;
            iss >> val;
            params.viscous = (val == "true" || val == "1");
        }
        else if (key == "output_freq") iss >> params.output_freq;
        else if (key == "restart_freq") iss >> params.restart_freq;
        else if (key == "case_name") iss >> params.case_name;
        else if (key == "flux_type") {
            std::string type;
            iss >> type;
            if (type == "DG") params.flux_type = FluxType::DG;
            else if (type == "SD") params.flux_type = FluxType::SD;
            else if (type == "HU") params.flux_type = FluxType::HU;
            else if (type == "GA") params.flux_type = FluxType::GA;
        }
        else if (key == "riemann") {
            std::string solver;
            iss >> solver;
            if (solver == "Rusanov") params.riemann = RiemannSolver::Rusanov;
            else if (solver == "Roe") params.riemann = RiemannSolver::Roe;
            else if (solver == "HLLC") params.riemann = RiemannSolver::HLLC;
        }
    }

    return params;
}

// Uniform flow initial condition
State uniformFlowIC(Vec2 pos, const SimParams& params) {
    IdealGas gas(params.gamma);

    Real rho = params.rho_inf;
    Real c = std::sqrt(params.gamma * params.p_inf / rho);
    Real vel = params.Mach_inf * c;
    Real u = vel * std::cos(params.AoA * M_PI / 180.0);
    Real v = vel * std::sin(params.AoA * M_PI / 180.0);
    Real p = params.p_inf;

    return gas.primToConserv(rho, u, v, p);
}

int main(int argc, char* argv[]) {
    // Initialize MPI
    MPICommunicator::init(&argc, &argv);
    MPICommunicator mpi;

    // Parse command line arguments
    ProgramOptions opts = parseArgs(argc, argv);

    if (opts.help) {
        if (mpi.isRoot()) {
            printUsage(argv[0]);
        }
        MPICommunicator::finalize();
        return 0;
    }

    // Print banner
    if (mpi.isRoot()) {
        std::cout << "======================================\n";
        std::cout << "  Zhijian - 2D FR Flow Solver\n";
        std::cout << "  GPU-Accelerated with MPI\n";
        std::cout << "======================================\n\n";
        std::cout << "Running on " << mpi.size() << " MPI processes\n\n";
    }

    // Load or create default parameters
    SimParams params;
    if (!opts.config_file.empty()) {
        params = parseConfig(opts.config_file);
        if (mpi.isRoot()) {
            std::cout << "Configuration loaded from: " << opts.config_file << "\n";
        }
    }
    params.output_dir = opts.output_dir;

    // Print parameters
    if (mpi.isRoot()) {
        std::cout << "\nSimulation Parameters:\n";
        std::cout << "  Polynomial order: P" << params.poly_order << "\n";
        std::cout << "  Gamma: " << params.gamma << "\n";
        std::cout << "  Mach: " << params.Mach_inf << "\n";
        std::cout << "  AoA: " << params.AoA << " deg\n";
        std::cout << "  Reynolds: " << params.Re << "\n";
        std::cout << "  CFL: " << params.CFL << "\n";
        std::cout << "  Viscous: " << (params.viscous ? "yes" : "no") << "\n";
        std::cout << "  Final time: " << params.t_final << "\n";
        std::cout << "  Max iterations: " << params.max_iter << "\n\n";
    }

    // Load and partition mesh
    Mesh mesh;
    if (!opts.mesh_file.empty()) {
        if (mpi.isRoot()) {
            std::cout << "Loading mesh: " << opts.mesh_file << "\n";
        }

        if (mpi.size() > 1) {
            // Parallel mesh distribution
            ParallelMesh par_mesh;
            par_mesh.distribute(opts.mesh_file, mpi);
            mesh = par_mesh.localMesh();
        } else {
            // Single process: read mesh directly
            CGNSReader reader;
            if (!reader.read(opts.mesh_file, mesh)) {
                std::cerr << "Error: " << reader.errorMessage() << std::endl;
                MPICommunicator::finalize();
                return 1;
            }
        }

        if (mpi.isRoot()) {
            mesh.printStats();
        }
    } else {
        if (mpi.isRoot()) {
            std::cerr << "Error: No mesh file specified. Use -m <file>\n";
            printUsage(argv[0]);
        }
        MPICommunicator::finalize();
        return 1;
    }

    // Initialize solver
    FRSolver solver;
    solver.initialize(mesh, params);

    // Set boundary conditions
    auto farfield_bc = FarFieldBC::fromMach(params.Mach_inf, params.AoA,
                                             params.rho_inf, params.p_inf,
                                             params.gamma);
    solver.setBoundaryCondition(1, farfield_bc);

    auto wall_bc = std::make_shared<SlipWallBC>(params.gamma);
    solver.setBoundaryCondition(2, wall_bc);

    // Set initial condition or load restart
    Real start_time = 0.0;
    int start_iter = 0;

    if (!opts.restart_file.empty()) {
        RestartIO restart_io;
        std::vector<ElementSolution> sol;
        SimParams restart_params;
        if (restart_io.read(opts.restart_file, mesh, sol, restart_params,
                            start_time, start_iter)) {
            // Would need to copy sol to solver
            if (mpi.isRoot()) {
                std::cout << "Restarting from: " << opts.restart_file << "\n";
                std::cout << "  Time: " << start_time << ", Iteration: " << start_iter << "\n";
            }
        }
    } else {
        // Uniform flow initial condition
        solver.setInitialCondition([&params](Vec2 pos) {
            return uniformFlowIC(pos, params);
        });
    }

    // Initialize output
    VTKWriter vtk_writer;
    RestartIO restart_io;
    ResidualWriter res_writer;

    if (mpi.isRoot()) {
        res_writer.open(params.output_dir + "/" + params.case_name + "_residual.csv");
    }

    // Write initial condition
    if (mpi.isRoot()) {
        std::string init_file = params.output_dir + "/" + params.case_name + "_0000.vtu";
        vtk_writer.writeVTU(init_file, mesh, solver.solution(), params);
        std::cout << "\nInitial condition written to: " << init_file << "\n\n";
    }

    // Main time stepping loop
    if (mpi.isRoot()) {
        std::cout << "Starting time integration...\n";
        std::cout << "-------------------------------------------\n";
        std::cout << "  Iter        Time          Residual\n";
        std::cout << "-------------------------------------------\n";
    }

    Real time = start_time;
    int output_count = 1;

    for (int iter = start_iter; iter < params.max_iter; ++iter) {
        // Compute time step
        Real dt_local = solver.computeTimeStep();
        Real dt = mpi.reduceMin(dt_local);

        if (params.dt > 0) {
            dt = std::min(dt, params.dt);
        }

        // Check if we've reached final time
        if (time + dt > params.t_final) {
            dt = params.t_final - time;
        }

        // Advance solution
        solver.advance(dt);
        time += dt;

        // Compute and report residual
        Real res_local = solver.computeResidual();
        Real residual = mpi.reduceSum(res_local * res_local);
        residual = std::sqrt(residual);

        if (mpi.isRoot() && iter % 10 == 0) {
            std::cout << std::setw(6) << iter << "  "
                      << std::scientific << std::setprecision(4)
                      << std::setw(12) << time << "  "
                      << std::setw(12) << residual << "\n";

            std::array<Real, N_VARS> comp_res = {0, 0, 0, 0};  // Would compute properly
            res_writer.write(iter, time, residual, comp_res);
        }

        // Output solution
        if (iter % params.output_freq == 0 && iter > start_iter) {
            solver.copyToHost();

            std::ostringstream fname;
            fname << params.output_dir << "/" << params.case_name
                  << "_" << std::setfill('0') << std::setw(4) << output_count << ".vtu";

            if (mpi.size() > 1) {
                vtk_writer.writeParallelVTU(params.output_dir + "/" + params.case_name +
                                             "_" + std::to_string(output_count),
                                             mpi.rank(), mpi.size(),
                                             mesh, solver.solution(), params);
            } else {
                vtk_writer.writeVTU(fname.str(), mesh, solver.solution(), params);
            }

            if (mpi.isRoot()) {
                std::cout << "  --> Output written: " << fname.str() << "\n";
            }
            output_count++;
        }

        // Write restart file
        if (iter % params.restart_freq == 0 && iter > start_iter) {
            solver.copyToHost();

            std::string restart_fname = params.output_dir + "/" + params.case_name + "_restart.h5";
            if (mpi.size() > 1) {
                restart_io.writeParallel(restart_fname, mpi.rank(), mpi.size(),
                                          mesh, solver.solution(), params, time, iter);
            } else {
                restart_io.write(restart_fname, mesh, solver.solution(), params, time, iter);
            }

            if (mpi.isRoot()) {
                std::cout << "  --> Restart written: " << restart_fname << "\n";
            }
        }

        // Check convergence or final time
        if (time >= params.t_final) {
            if (mpi.isRoot()) {
                std::cout << "\nFinal time reached.\n";
            }
            break;
        }
    }

    // Final output
    solver.copyToHost();

    std::string final_file = params.output_dir + "/" + params.case_name + "_final.vtu";
    if (mpi.size() > 1) {
        vtk_writer.writeParallelVTU(params.output_dir + "/" + params.case_name + "_final",
                                     mpi.rank(), mpi.size(),
                                     mesh, solver.solution(), params);
    } else {
        vtk_writer.writeVTU(final_file, mesh, solver.solution(), params);
    }

    // Final restart
    std::string final_restart = params.output_dir + "/" + params.case_name + "_final.h5";
    if (mpi.size() > 1) {
        restart_io.writeParallel(final_restart, mpi.rank(), mpi.size(),
                                  mesh, solver.solution(), params, time, solver.iteration());
    } else {
        restart_io.write(final_restart, mesh, solver.solution(), params, time, solver.iteration());
    }

    if (mpi.isRoot()) {
        std::cout << "\n===========================================\n";
        std::cout << "Simulation completed successfully!\n";
        std::cout << "  Final time: " << time << "\n";
        std::cout << "  Iterations: " << solver.iteration() << "\n";
        std::cout << "  Output: " << final_file << "\n";
        std::cout << "  Restart: " << final_restart << "\n";
        std::cout << "===========================================\n";
    }

    res_writer.close();
    MPICommunicator::finalize();

    return 0;
}
