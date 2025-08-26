using Elmer
using OrderedCollections

# Parameter definition
frequency = 50

# Set up simulation
simulation = Simulation(7, Elmer.CoordinateCartesian(), Elmer.SimulationSteady(), 1, OrderedDict("Frequency" => frequency))
sif = SolverInformationFile("case", simulation, data_path="examples/TEAM7/simdata/", mesh_db="team7")

solver_coil = load_solver!(sif, "CoilSolver", Elmer.ExecBeforeAll(), "examples/TEAM7/solvers.yml")
solver_mgharm = load_solver!(sif, "MGDynHarm", Elmer.ExecAlways(), "examples/TEAM7/solvers.yml")
solver_mgpost = load_solver!(sif, "MGDynPost", Elmer.ExecAlways(), "examples/TEAM7/solvers.yml")
solver_result = load_solver!(sif, "ResultOutput", Elmer.ExecAfterStep(), "examples/TEAM7/solvers.yml")

eq_coil = add_equation!(sif, "coil", [solver_coil, solver_mgharm, solver_mgpost, solver_result])
eq_domain = add_equation!(sif, "domain", [solver_mgharm, solver_mgpost, solver_result])

# Define physical bodies
mat_air = Elmer.load_material!(sif, "Air", "examples/TEAM7/materials.yml")
mat_aluminium = Elmer.load_material!(sif, "Aluminium", "examples/TEAM7/materials.yml")

body_coil = add_body!(sif, "Coil", [2]; equation=eq_coil, material=mat_air)
body_plate = add_body!(sif, "Plate", [1]; equation=eq_domain, material=mat_aluminium)
body_domain = add_body!(sif, "Domain", [3]; equation=eq_domain, material=mat_air)

bc_inf = add_boundary_condition!(sif, "Inf", data=OrderedDict("AV {e}" => 0), [1])

component = add_component!(sif, "Coil", data=OrderedDict("Coil Type" => "\"Stranded\"", "Number of Turns" => 1, "Desired Coil Current" => -2742); master_bodies=[body_coil])

# Write SIF & run
Elmer.write(sif)
Elmer.write_startinfo(sif)

# Elmer.elmergrid_gmsh(sif.data_path, "team7.msh")
Elmer.run_elmer_solver(sif)