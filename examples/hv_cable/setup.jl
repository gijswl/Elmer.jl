using Elmer
using OrderedCollections

data_path = joinpath(@__DIR__, "simdata\\")
solvers_db = joinpath(@__DIR__, "solvers.yml")
materials_db = joinpath(@__DIR__, "materials.yml")

frequency = 50
Icond1 = 1000 * exp(0π / 3 * 1im)
Icond2 = 1000 * exp(2π / 3 * 1im)
Icond3 = 1000 * exp(-2π / 3 * 1im)

# Electromagnetic simulation
sim_em = Simulation(7, Elmer.CoordinateCartesian(), Elmer.SimulationSteady(), 1e-3, OrderedDict(
    "Timestep Intervals" => length(frequency),
    "Output Intervals" => 1,
    "Steady State Max Iterations" => 1,
    "Frequency" => frequency
))

sif_em = SolverInformationFile("electromagnetic", sim_em, data_path=data_path, mesh_db="electromagnetic")

solver_circuit = load_solver!(sif_em, "CircuitSolver", Elmer.ExecAlways(), solvers_db)
solver_mgharm = load_solver!(sif_em, "MGDynHarm2D", Elmer.ExecAlways(), solvers_db)
post_mgharm = load_solver!(sif_em, "MGDynPost2D", Elmer.ExecAlways(), solvers_db)
post_circuit = load_solver!(sif_em, "CircuitPost", Elmer.ExecAlways(), solvers_db)
save_result = load_solver!(sif_em, "ResultOutput", Elmer.ExecAfterStep(), solvers_db)
save_scalars = load_solver!(sif_em, "SaveScalars", Elmer.ExecAlways(), solvers_db)

update_solver_data!(sif_em, solver_mgharm, "Quadratic Approximation", true)

update_solver_data!(sif_em, save_scalars, "Filename", "electromagnetic.dat")
update_solver_data!(sif_em, save_scalars, OrderedDict(
    "Variable 1" => "joule heating e",
    "Mask Name 1" => "Conductor1",
    "Operator 1" => "body int",
    "Variable 2" => "joule heating e",
    "Mask Name 2" => "Conductor2",
    "Operator 2" => "body int",
    "Variable 3" => "joule heating e",
    "Mask Name 3" => "Conductor3",
    "Operator 3" => "body int",
    "Variable 4" => "joule heating e",
    "Mask Name 4" => "Sheath",
    "Operator 4" => "body int"
))

update_solver_data!(sif_em, save_result, "Output File Name", "electromagnetic")

eq = add_equation!(sif_em, "Domain", [solver_circuit, solver_mgharm, post_mgharm, post_circuit, save_result, save_scalars])

mat_air = load_material!(sif_em, "Air", materials_db)
mat_Al = load_material!(sif_em, "Aluminium", materials_db)
mat_Cu = load_material!(sif_em, "Copper", materials_db)
mat_XLPE = load_material!(sif_em, "XLPE", materials_db)
mat_HDPE = load_material!(sif_em, "HDPE", materials_db)

body_conductor1 = add_body!(sif_em, "Conductor1", [1]; equation=eq, material=mat_Al, mask = ["Conductor1"])
body_conductor2 = add_body!(sif_em, "Conductor2", [2]; equation=eq, material=mat_Al, mask = ["Conductor2"])
body_conductor3 = add_body!(sif_em, "Conductor3", [3]; equation=eq, material=mat_Al, mask = ["Conductor3"])
body_insulation = add_body!(sif_em, "Insulation", [4]; equation=eq, material=mat_XLPE)
body_sheath = add_body!(sif_em, "Sheath", [5]; equation=eq, material=mat_Cu, mask = ["Sheath"])
body_jacket = add_body!(sif_em, "Jacket", [6]; equation=eq, material=mat_HDPE)
body_air = add_body!(sif_em, "Air", [7]; equation=eq, material=mat_air)

bc_inf = add_boundary_condition!(sif_em, "Inf", [1]; data=OrderedDict("Infinity BC" => true))

Elmer.write(sif_em)

# Elmer.elmergrid_gmsh(sif_em.data_path, "electromagnetic.msh")
# Elmer.run_elmer_solver(sif_em)

# TODO read loss from .dat