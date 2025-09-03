using Elmer
using OrderedCollections

# Parameter definition
# frequency = [50, 100, 200]
frequency = [50]

# Set up simulation
simulation = Simulation(7, Elmer.CoordinateCartesian(), Elmer.SimulationScanning(), 1, OrderedDict(
    "Timestep Intervals" => length(frequency),
    "Output Intervals" => 1,
    "Steady State Max Iterations" => 1,
    "\$ f" => join(frequency, " "),
    "Frequency" => "Variable timestep; Real MATC \"f(tx-1)\""
))

data_path = joinpath(@__DIR__, "simdata\\")
solvers_db = joinpath(@__DIR__, "solvers.yml")
materials_db = joinpath(@__DIR__, "materials.yml")

sif = SolverInformationFile("case", simulation, data_path=data_path, mesh_db="team7")

solver_coil = load_solver!(sif, "CoilSolver", Elmer.ExecBeforeAll(), solvers_db)
solver_mgharm = load_solver!(sif, "MGDynHarm", Elmer.ExecAlways(), solvers_db)
post_mgharm = load_solver!(sif, "MGDynPost", Elmer.ExecAlways(), solvers_db)
save_scalars = load_solver!(sif, "SaveScalars", Elmer.ExecAlways(), solvers_db)
save_result = load_solver!(sif, "ResultOutput", Elmer.ExecAfterStep(), solvers_db)
update_solver_data!(sif, save_scalars, "Filename", "scalars.dat")
update_solver_data!(sif, save_result, "Output File Name", "case")

eq_coil = add_equation!(sif, "coil", [solver_coil, solver_mgharm, post_mgharm, save_scalars])
eq_domain = add_equation!(sif, "domain", [solver_mgharm, post_mgharm])

# Define physical bodies
mat_coil = add_material!(sif, "Coil", data=OrderedDict("Relative Permittivity" => 1, "Relative Permeability" => 1))
mat_aluminium = Elmer.load_material!(sif, "Aluminium", materials_db)
mat_air = Elmer.load_material!(sif, "Air", materials_db)

body_coil = add_body!(sif, "Coil", [2]; equation=eq_coil, material=mat_coil)
body_plate = add_body!(sif, "Plate", [1]; equation=eq_domain, material=mat_aluminium)
body_domain = add_body!(sif, "Domain", [3]; equation=eq_domain, material=mat_air)

bc_inf = add_boundary_condition!(sif, "Inf", [1]; data=OrderedDict("AV {e}" => 0))

component = add_component!(sif, "Coil", data=OrderedDict(
        "Coil Type" => "Stranded",
        "Number of Turns" => 1,
        "Desired Coil Current" => 2742,
        "Coil Use W Vector" => true,
        "W Vector Variable Name" => "String \"CoilCurrent e\""); master_bodies=[body_coil])

# Write SIF & run
Elmer.write(sif)
Elmer.write_startinfo(sif)

# Elmer.elmergrid_gmsh(sif.data_path, "team7.msh")
Elmer.run_elmer_solver(sif)