using Elmer
using OrderedCollections

using Unitful
using PhysicalConstants.CODATA2022

# Parameter definition
f = 2.5e9
a = 10e-2
b = 5e-2

ω = 2π * f
c0 = ustrip(SpeedOfLightInVacuum)
k0 = ω / c0
kc = π / a
β0 = √(k0^2 - kc^2)

# Set up simulation
data_path = joinpath(@__DIR__, "simdata\\")
solvers_db = joinpath(@__DIR__, "solvers.yml")

simulation = Simulation(7, Elmer.CoordinateCartesian(), Elmer.SimulationSteady())
update_simulation_data!(simulation, Elmer.format_frequency(f))

sif = SolverInformationFile("case", simulation, data_path=data_path, mesh_db="bent_waveguide")

solver_helmholtz = load_solver!(sif, "VectorHelmholtzSolver", Elmer.ExecAlways(), solvers_db)
update_solver_data!(sif, solver_helmholtz, "Linear System Preconditioning", "ILUT")
update_solver_data!(sif, solver_helmholtz, "Quadratic Approximation", true)
post_helmholtz = load_solver!(sif, "VectorHelmholtzPost", Elmer.ExecAlways(), solvers_db)
solver_result = load_solver!(sif, "ResultOutputSolver", Elmer.ExecAfterStep(), solvers_db)

eqn = add_equation!(sif, "main", [solver_helmholtz, post_helmholtz, solver_result])

# Define physical bodies
mat_air = add_material!(sif, "Air"; data=OrderedDict("Relative Permittivity" => 1.0, "Relative Permeability" => 1.0))

bdy_wg = add_body!(sif, "Waveguide", [1]; equation=eqn, material=mat_air)

bnd_pec = add_boundary_condition!(sif, "PEC", [3]; data=OrderedDict("E re {e}" => 0.0, "E im {e}" => 0.0))
bnd_in = add_boundary_condition!(sif, "Port in", [1]; data=OrderedDict(
    "Magnetic Boundary Load 2" => "Variable Coordinate 1; Real MATC \"-2 * $β0 * $k0 / $kc * sin($kc * (tx + $a/2))\"",
    "Electric Robin Coefficient im" => β0))
bnd_out = add_boundary_condition!(sif, "Port out", [2]; data=OrderedDict("Electric Robin Coefficient im" => β0))

# Write SIF & run
Elmer.write(sif)

Elmer.elmergrid_gmsh(sif.data_path, "bent_waveguide.msh")
Elmer.elmer_solver(sif)
