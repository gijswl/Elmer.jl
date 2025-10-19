using Elmer
using OrderedCollections
using PyCall

cb = pyimport("elmer_circuitbuilder")

data_path = joinpath(@__DIR__, "simdata\\")
solvers_db = joinpath(@__DIR__, "solvers.yml")
materials_db = joinpath(@__DIR__, "materials.yml")
circuit_file = "circuit.definition"

frequency = 50
Icond = 1000

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

body_conductor1 = add_body!(sif_em, "Conductor1", [1]; equation=eq, material=mat_Al, mask=["Conductor1"])
body_conductor2 = add_body!(sif_em, "Conductor2", [2]; equation=eq, material=mat_Al, mask=["Conductor2"])
body_conductor3 = add_body!(sif_em, "Conductor3", [3]; equation=eq, material=mat_Al, mask=["Conductor3"])
body_insulation = add_body!(sif_em, "Insulation", [4]; equation=eq, material=mat_XLPE)
body_sheath = add_body!(sif_em, "Sheath", [5]; equation=eq, material=mat_Cu, mask=["Sheath"])
body_jacket = add_body!(sif_em, "Jacket", [6]; equation=eq, material=mat_HDPE)
body_air = add_body!(sif_em, "Air", [7]; equation=eq, material=mat_air)

bc_inf = add_boundary_condition!(sif_em, "Inf", [1]; data=OrderedDict("Infinity BC" => true))

# Circuit definition
c = cb.number_of_circuits(1)

c[1].ref_node = 1
I1 = cb.I("I1", 1, 2, Icond * exp(+0π / 3 * 1im))
I2 = cb.I("I2", 1, 3, Icond * exp(+2π / 3 * 1im))
I3 = cb.I("I3", 1, 4, Icond * exp(-2π / 3 * 1im))
cond1 = cb.ElmerComponent("T1", 1, 2, 1, [PyObject(body_conductor1)])
cond2 = cb.ElmerComponent("T2", 1, 3, 2, [PyObject(body_conductor2)])
cond3 = cb.ElmerComponent("T3", 1, 4, 3, [PyObject(body_conductor3)])
c[1].components = [[I1, I2, I3, cond1, cond2, cond3]]

cb.generate_elmer_circuits(c, joinpath(data_path, circuit_file))
add_include!(sif_em, circuit_file)

Elmer.write(sif_em)

# Elmer.elmergrid_gmsh(sif_em.data_path, "electromagnetic.msh")
Elmer.run_elmer_solver(sif_em)

dat_em = load_dat(joinpath(data_path, "results", "electromagnetic.dat"))
Pj_cond1 = dat_em[1, "body int: joule heating e mask conductor1"]
Pj_cond2 = dat_em[1, "body int: joule heating e mask conductor2"]
Pj_cond3 = dat_em[1, "body int: joule heating e mask conductor3"]
Pj_sheath = dat_em[1, "body int: joule heating e mask sheath"]
Peddy = dat_em[1, "res: eddy current power"]

@assert Pj_cond1 + Pj_cond2 + Pj_cond3 + Pj_sheath ≈ Peddy "Total eddy current losses should be equal to the sum of the defined joule heating body integrals"