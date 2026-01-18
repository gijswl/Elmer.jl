using Elmer
using Elmer.CircuitBuilder
using OrderedCollections

data_path = joinpath(@__DIR__, "simdata\\")
solvers_db = joinpath(@__DIR__, "solvers.yml")
materials_db = joinpath(@__DIR__, "materials.yml")

Elmer.elmergrid_gmsh(data_path, "electromagnetic.msh")
Elmer.elmergrid_gmsh(data_path, "thermal.msh")

"""
    setup_em(frequency::F, Icond::Real; sim_name::String="electromagnetic", circuit_file::String = "circuit.definition") where {F<:Union{Real,Vector{<:Real}}}

Set up the electromagnetic part of the model
"""
function setup_em(frequency::F, Icond::Real; sim_name::String="electromagnetic", circuit_file::String="circuit.definition") where {F <: Union{Real, Vector{<:Real}}}
    # Simulation parameters
    coord_scaling = 1e-3
    sim_em = Simulation(7, Elmer.CoordinateCartesian(), Elmer.SimulationScanning(); coordinate_scale=coord_scaling, data=OrderedDict{String, Any}(
        "Output Intervals" => 1,
        "Steady State Max Iterations" => 1,
    ))
    update_simulation_data!(sim_em, Elmer.format_frequency(frequency))

    sif_em = SolverInformationFile(sim_name, sim_em, data_path=data_path, mesh_db=sim_name)

    # Solvers
    solver_circuit = load_solver!(sif_em, "CircuitSolver", Elmer.ExecAlways(), solvers_db)
    solver_mgharm = load_solver!(sif_em, "MGDynHarm2D", Elmer.ExecAlways(), solvers_db)
    post_mgharm = load_solver!(sif_em, "MGDynPost2D", Elmer.ExecAlways(), solvers_db)
    post_circuit = load_solver!(sif_em, "CircuitPost", Elmer.ExecAlways(), solvers_db)
    save_result = load_solver!(sif_em, "ResultOutput", Elmer.ExecAfterStep(), solvers_db)
    save_scalars = load_solver!(sif_em, "SaveScalars", Elmer.ExecAlways(), solvers_db)

    update_solver_data!(sif_em, solver_mgharm, "Quadratic Approximation", true)

    update_solver_data!(sif_em, save_scalars, "Filename", "$sim_name.dat")
    update_solver_data!(
        sif_em,
        save_scalars,
        OrderedDict(
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
            "Operator 4" => "body int",
        ),
    )

    update_solver_data!(sif_em, save_result, "Output File Name", sim_name)

    # Equations
    eq = add_equation!(sif_em, "Domain", [solver_circuit, solver_mgharm, post_mgharm, post_circuit, save_result, save_scalars])

    # Materials
    mat_air = load_material!(sif_em, "Air", materials_db)
    mat_Al = load_material!(sif_em, "Aluminium", materials_db)
    mat_Cu = load_material!(sif_em, "Copper", materials_db)
    mat_XLPE = load_material!(sif_em, "XLPE", materials_db)
    mat_HDPE = load_material!(sif_em, "HDPE", materials_db)

    # Bodies
    body_conductor1 = add_body!(sif_em, "Conductor1", [1]; equation=eq, material=mat_Al, mask=["Conductor1"])
    body_conductor2 = add_body!(sif_em, "Conductor2", [2]; equation=eq, material=mat_Al, mask=["Conductor2"])
    body_conductor3 = add_body!(sif_em, "Conductor3", [3]; equation=eq, material=mat_Al, mask=["Conductor3"])
    body_insulation = add_body!(sif_em, "Insulation", [4]; equation=eq, material=mat_XLPE)
    body_sheath = add_body!(sif_em, "Sheath", [5]; equation=eq, material=mat_Cu, mask=["Sheath"])
    body_jacket = add_body!(sif_em, "Jacket", [6]; equation=eq, material=mat_HDPE)
    body_air = add_body!(sif_em, "Air", [7]; equation=eq, material=mat_air)

    # Boundry conditions
    bc_inf = add_boundary_condition!(sif_em, "Inf", [1]; data=OrderedDict("Infinity BC" => true))

    # Circuit definition
    c = create_circuits(1)

    insert_component!(c[1], CurrentSource("I1", (1, 2), Icond * exp(+0im * 2π / 3)))
    insert_component!(c[1], CurrentSource("I2", (1, 3), Icond * exp(+1im * 2π / 3)))
    insert_component!(c[1], CurrentSource("I3", (1, 4), Icond * exp(-1im * 2π / 3)))
    insert_component!(c[1], ElmerComponent("T1", (1, 2), 1, [body_conductor1]))
    insert_component!(c[1], ElmerComponent("T2", (1, 3), 2, [body_conductor2]))
    insert_component!(c[1], ElmerComponent("T3", (1, 4), 3, [body_conductor3]))

    write_circuits(c, circuit_file; path=data_path)

    add_include!(sif_em, circuit_file)
    components, body_forces = add_circuits!(sif_em, c)

    add_coil_data!(sif_em, 1, CoilMassive())
    add_coil_data!(sif_em, 2, CoilMassive())
    add_coil_data!(sif_em, 3, CoilMassive())

    # Generate SIF
    Elmer.write(sif_em)

    return sif_em
end

function setup_thermal(loss::Dict; sim_name::String="thermal")
    coord_scaling = 1e-3
    sim_thermal = Simulation(7, Elmer.CoordinateCartesian(), Elmer.SimulationSteady(); coordinate_scale=coord_scaling, data=OrderedDict{String, Any}(
        "Output Intervals" => 1,
        "Steady State Max Iterations" => 1,
    ))

    sif_thermal = SolverInformationFile(sim_name, sim_thermal, data_path=data_path, mesh_db=sim_name)

    # Solvers
    solver_thermal = load_solver!(sif_thermal, "HeatSolver", Elmer.ExecAlways(), solvers_db)
    save_result = load_solver!(sif_thermal, "ResultOutput", Elmer.ExecAfterStep(), solvers_db)

    update_solver_data!(sif_thermal, save_result, "Output File Name", sim_name)

    # Equations
    eq = add_equation!(sif_thermal, "Domain", [solver_thermal, save_result])

    # Materials
    mat_soil = load_material!(sif_thermal, "Soil", materials_db)
    mat_Al = load_material!(sif_thermal, "Aluminium", materials_db)
    mat_XLPE = load_material!(sif_thermal, "XLPE", materials_db)
    mat_Cu = load_material!(sif_thermal, "Copper", materials_db)
    mat_HDPE = load_material!(sif_thermal, "HDPE", materials_db)

    # Body forces
    bf_conductor1 = add_body_force!(sif_thermal, "Conductor1"; data=OrderedDict("Heat Source" => 1.0, "Integral Heat Source" => loss["Conductor1"]))
    bf_conductor2 = add_body_force!(sif_thermal, "Conductor2"; data=OrderedDict("Heat Source" => 1.0, "Integral Heat Source" => loss["Conductor2"]))
    bf_conductor3 = add_body_force!(sif_thermal, "Conductor3"; data=OrderedDict("Heat Source" => 1.0, "Integral Heat Source" => loss["Conductor3"]))
    bf_sheath = add_body_force!(sif_thermal, "Sheath"; data=OrderedDict("Heat Source" => 1.0, "Integral Heat Source" => loss["Sheath"]))
    bf_dielectric = add_body_force!(sif_thermal, "Dielectric"; data=OrderedDict("Heat Source" => 1.0, "Integral Heat Source" => loss["Dielectric"]))

    # Bodies
    body_conductor1 = add_body!(sif_thermal, "Conductor1", [1]; equation=eq, material=mat_Al, body_force=bf_conductor1)
    body_conductor2 = add_body!(sif_thermal, "Conductor2", [2]; equation=eq, material=mat_Al, body_force=bf_conductor2)
    body_conductor3 = add_body!(sif_thermal, "Conductor3", [3]; equation=eq, material=mat_Al, body_force=bf_conductor3)
    body_insulation = add_body!(sif_thermal, "Insulation", [4]; equation=eq, material=mat_XLPE, body_force=bf_dielectric)
    body_sheath = add_body!(sif_thermal, "Sheath", [5]; equation=eq, material=mat_Cu, body_force=bf_sheath)
    body_jacket = add_body!(sif_thermal, "Jacket", [6]; equation=eq, material=mat_HDPE)
    body_air = add_body!(sif_thermal, "Soil", [7]; equation=eq, material=mat_soil)

    # Boundary conditions
    bc_surf = add_boundary_condition!(sif_thermal, "Surface", [1]; data=OrderedDict("Temperature" => 40.0))
    bc_inf = add_boundary_condition!(sif_thermal, "Inf", [2]; data=OrderedDict("Infinity BC Temp" => true))

    # Generate sif
    Elmer.write(sif_thermal)

    return sif_thermal
end

frequency = 50.0
Icond = 515.0

sif_em = setup_em(frequency, Icond)
Elmer.elmer_solver(sif_em, log_file="solver_em.log")

dat_em = load_dat(joinpath(data_path, "results", "electromagnetic.dat"))
Pj_cond1 = dat_em[1, "body int: joule heating e mask conductor1"]
Pj_cond2 = dat_em[1, "body int: joule heating e mask conductor2"]
Pj_cond3 = dat_em[1, "body int: joule heating e mask conductor3"]
Pj_sheath = dat_em[1, "body int: joule heating e mask sheath"]
Peddy = dat_em[1, "res: eddy current power"]

@assert Pj_cond1 + Pj_cond2 + Pj_cond3 + Pj_sheath ≈ Peddy "Total eddy current losses should be equal to the sum of the defined joule heating body integrals"

Pj_diel = 3 * 0.128; # TODO calculate based on voltage and tan δ

sif_thermal = setup_thermal(Dict(
    "Conductor1" => Pj_cond1,
    "Conductor2" => Pj_cond2,
    "Conductor3" => Pj_cond3,
    "Sheath" => Pj_sheath,
    "Dielectric" => Pj_diel,
))
Elmer.elmer_solver(sif_thermal, log_file="solver_thermal.log")
