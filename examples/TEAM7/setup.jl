using Elmer
using OrderedCollections

# Parameter definition
frequency = [50, 100, 200]

# Set up simulation
simulation = Simulation(7, Elmer.CoordinateCartesian(), Elmer.SimulationScanning(); data = OrderedDict{String, Any}(
	"Output Intervals" => 1,
	"Steady State Max Iterations" => 1,
))
update_simulation_data!(simulation, Elmer.format_frequency(frequency))

data_path = joinpath(@__DIR__, "simdata\\")
solvers_db = joinpath(@__DIR__, "solvers.yml")
materials_db = joinpath(@__DIR__, "materials.yml")

sif = SolverInformationFile("case", simulation, data_path = data_path, mesh_db = "team7")

solver_coil = load_solver!(sif, "CoilSolver", Elmer.ExecBeforeAll(), solvers_db)
solver_mgharm = load_solver!(sif, "MGDynHarm", Elmer.ExecAlways(), solvers_db)
update_solver_data!(sif, solver_mgharm, "Linear System Nullify Guess", true)

post_mgharm = load_solver!(sif, "MGDynPost", Elmer.ExecAlways(), solvers_db)
save_scalars = load_solver!(sif, "SaveScalars", Elmer.ExecAlways(), solvers_db)
save_result = load_solver!(sif, "ResultOutput", Elmer.ExecAfterStep(), solvers_db)
update_solver_data!(sif, save_scalars, "Filename", "scalars.dat")
update_solver_data!(sif, save_result, "Output File Name", "case")

eq_coil = add_equation!(sif, "coil", [solver_coil, solver_mgharm, post_mgharm, save_scalars])
eq_domain = add_equation!(sif, "domain", [solver_mgharm, post_mgharm])

# Define physical bodies
mat_coil = add_material!(sif, "Coil", data = OrderedDict("Relative Permittivity" => 1.0, "Relative Permeability" => 1.0))
mat_aluminium = Elmer.load_material!(sif, "Aluminium", materials_db)
mat_air = Elmer.load_material!(sif, "Air", materials_db)

body_coil = add_body!(sif, "Coil", [2]; equation = eq_coil, material = mat_coil)
body_plate = add_body!(sif, "Plate", [1]; equation = eq_domain, material = mat_aluminium)
body_domain = add_body!(sif, "Domain", [3]; equation = eq_domain, material = mat_air)

bc_inf = add_boundary_condition!(sif, "Inf", [1]; data = OrderedDict("AV {e}" => 0.0))

component = add_component!(sif, "Coil", master_bodies = [body_coil], data = OrderedDict(
	"Coil Type" => "Stranded",
	"Number of Turns" => 1.0,
	"Desired Coil Current" => 2742.0,
	"Coil Use W Vector" => true,
	"W Vector Variable Name" => "String \"CoilCurrent e\""))

# Write SIF & run
Elmer.write(sif)

# Elmer.elmergrid_gmsh(sif.data_path, "team7.msh")
Elmer.elmer_solver(sif)

dat_em = load_dat(joinpath(data_path, "results", "scalars.dat"))
