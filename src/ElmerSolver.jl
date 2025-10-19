import Base: write

abstract type CoordinateSystem end;
struct CoordinateCartesian <: CoordinateSystem end
struct CoordinateAxi <: CoordinateSystem end
struct CoordinateCylindric <: CoordinateSystem end

abstract type SimulationType end;
struct SimulationSteady <: SimulationType end
struct SimulationTransient <: SimulationType end
struct SimulationScanning <: SimulationType end

mutable struct Simulation
    output_level::Int
    coordinate_system::CoordinateSystem
    simulation_type::SimulationType
    coordinate_scaling::Real
    data::OrderedDict
end

abstract type ExecSolver end
struct ExecAlways <: ExecSolver end
struct ExecBeforeAll <: ExecSolver end
struct ExecAfterAll <: ExecSolver end
struct ExecBeforeSim <: ExecSolver end
struct ExecAfterSim <: ExecSolver end
struct ExecBeforeStep <: ExecSolver end
struct ExecAfterStep <: ExecSolver end
struct ExecBeforeSave <: ExecSolver end
struct ExecAfterSave <: ExecSolver end
struct ExecNever <: ExecSolver end

mutable struct Solver
    id::Int
    name::String
    exec::ExecSolver
    data::OrderedDict
end

mutable struct Equation
    id::Int
    name::String
    solvers::Vector{Int}
end

mutable struct Material
    id::Int
    name::String
    data::OrderedDict
end

mutable struct Body
    id::Int
    name::String
    target_bodies::Vector{Int}
    equation::Union{Int,Missing}
    material::Union{Int,Missing}
    body_force::Union{Int,Missing}
    initial_condition::Union{Int,Missing}
    data::OrderedDict
end

mutable struct BodyForce
    id::Int
    name::String
    data::OrderedDict
end

mutable struct InitialCondition
    id::Int
    name::String
    data::OrderedDict
end

mutable struct BoundaryCondition
    id::Int
    name::String
    target_boundaries::Vector{Int}
    data::OrderedDict
end

mutable struct Component
    id::Int
    name::String
    master_bodies::Vector{Int}
    master_boundaries::Vector{Int}
    data::OrderedDict
end

mutable struct SolverInformationFile
    filename::String
    data_path::String
    result_dir::String
    mesh_db::String

    simulation::Simulation
    constants::Dict
    includes::Vector{String}

    equations::Vector{Equation}
    solvers::Vector{Solver}
    materials::Vector{Material}
    bodies::Vector{Body}
    body_forces::Vector{BodyForce}
    initial_conditions::Vector{InitialCondition}
    boundary_conditions::Vector{BoundaryCondition}
    components::Vector{Component}

    function SolverInformationFile(filename::String, simulation::Simulation; data_path::String="simdata/", result_dir::String="results", mesh_db::String=".")
        return new(filename * ".sif", data_path, result_dir, mesh_db, simulation, get_default_constants(), String[], Equation[], Solver[], Material[], Body[], BodyForce[], InitialCondition[], BoundaryCondition[], Component[])
    end
end

include("Database.jl")

format_property(name::String, value::String) = name * " = " * value
format_property(name::String, len::Int, value::String) = name * "(" * string(len) * ") = " * value
format_property(name::String, value::Any) = format_property(name, format_value(value))
format_property(name::String, value::Vector) = format_property(name, length(value), format_value(value))

format_value(value::String) = value
format_value(value::Real) = string(value)
format_value(value::Bool) = value ? "Logical True" : "Logical False"
format_value(value::Vector) = join(format_value.(value), " ")

format_value(value::ExecAlways) = "Always"
format_value(value::ExecBeforeAll) = "Before All"
format_value(value::ExecAfterAll) = "After All"
format_value(value::ExecBeforeSim) = "Before Simulation"
format_value(value::ExecAfterSim) = "After Simulation"
format_value(value::ExecBeforeStep) = "Before Timestep"
format_value(value::ExecAfterStep) = "After Timestep"
format_value(value::ExecBeforeSave) = "Before Saving"
format_value(value::ExecAfterSave) = "After Saving"
format_value(value::ExecNever) = "Never"

format_value(value::CoordinateCartesian) = "Cartesian"
format_value(value::CoordinateAxi) = "Axi Symmetric"
format_value(value::CoordinateCylindric) = "Cylindric Symmetric"

format_value(value::SimulationSteady) = "Steady"
format_value(value::SimulationTransient) = "Transient"
format_value(value::SimulationScanning) = "Scanning"

function add_constant!(sif::SolverInformationFile, name::String, value)
    sif.constants[name] = value
end

function add_include!(sif::SolverInformationFile, file::String)
    push!(sif.includes, file)
end

function add_solver!(sif::SolverInformationFile, name::String, exec::ExecSolver; data::OrderedDict=OrderedDict())
    id = length(sif.solvers) + 1
    solver = Solver(id, name, exec, data)
    push!(sif.solvers, solver)
    return id
end

function add_equation!(sif::SolverInformationFile, name::String, solvers::Vector{Int})
    id = length(sif.equations) + 1
    eq = Equation(id, name, solvers)
    push!(sif.equations, eq)
    return id
end

function add_material!(sif::SolverInformationFile, name::String; data::OrderedDict=OrderedDict())
    id = length(sif.materials) + 1
    material = Material(id, name, data)
    push!(sif.materials, material)
    return id
end

function add_body!(sif::SolverInformationFile, name::String, target_bodies::Vector{Int}=Int[]; equation=missing, material=missing, body_force=missing, initial_condition=missing, mask=missing, data::OrderedDict=OrderedDict())
    id = length(sif.bodies) + 1

    if (~ismissing(mask))
        for m ∈ mask
            data[m] = true
        end
    end

    body = Body(id, name, target_bodies, equation, material, body_force, initial_condition, data)
    push!(sif.bodies, body)
    return id
end

function add_body_force!(sif::SolverInformationFile, name::String; data::OrderedDict=OrderedDict())
    id = length(sif.body_forces) + 1
    body_force = BodyForce(id, name, data)
    push!(sif.body_forces, body_force)
    return id
end

function add_initial_condition!(sif::SolverInformationFile, name::String; data::OrderedDict=OrderedDict())
    id = length(sif.initial_conditions) + 1
    initial_condition = InitialCondition(id, name, data)
    push!(sif.initial_conditions, initial_condition)
    return id
end

function add_boundary_condition!(sif::SolverInformationFile, name::String, target_boundaries::Vector{Int}=Int[]; data::OrderedDict=OrderedDict())
    id = length(sif.boundary_conditions) + 1
    boundary_condition = BoundaryCondition(id, name, target_boundaries, data)
    push!(sif.boundary_conditions, boundary_condition)
    return id
end

function add_component!(sif::SolverInformationFile, name::String; master_bodies::Vector{Int}=Int[], master_boundaries::Vector{Int}=Int[], data::OrderedDict=OrderedDict())
    id = length(sif.components) + 1
    component = Component(id, name, master_bodies, master_boundaries, data)
    push!(sif.components, component)
    return id
end

function update_solver_data!(sif::SolverInformationFile, solver_id::Int, data::OrderedDict)
    solver = sif.solvers[solver_id]
    for (key, val) ∈ data
        solver.data[key] = val
    end
end

function update_solver_data!(sif::SolverInformationFile, solver_id::Int, key::String, value)
    solver = sif.solvers[solver_id]
    solver.data[key] = value
end

function update_material_data!(sif::SolverInformationFile, material_id::Int, data::OrderedDict)
    material = sif.materials[material_id]
    for (key, val) ∈ data
        material.data[key] = val
    end
end

function update_material_data!(sif::SolverInformationFile, material_id::Int, key::String, value)
    material = sif.materials[material_id]
    material.data[key] = value
end

function update_body_data!(sif::SolverInformationFile, body_id::Int, data::OrderedDict)
    body = sif.bodies[body_id]
    for (key, val) ∈ data
        body.data[key] = val
    end
end

function update_body_data!(sif::SolverInformationFile, body_id::Int, key::String, value)
    body = sif.bodies[body_id]
    body.data[key] = value
end

function update_bf_data!(sif::SolverInformationFile, bf_id::Int, key::String, value)
    body_force = sif.body_forces[bf_id]
    body_force.data[key] = value
end

function update_ic_data!(sif::SolverInformationFile, ic_id::Int, key::String, value)
    initial_condition = sif.initial_conditions[ic_id]
    initial_condition.data[key] = value
end

function update_bc_data!(sif::SolverInformationFile, bc_id::Int, key::String, value)
    boundary_condition = sif.boundary_conditions[bc_id]
    boundary_condition.data[key] = value
end

function update_component_data!(sif::SolverInformationFile, component_id::Int, key::String, value)
    component = sif.components[component_id]
    component.data[key] = value
end

function write_header(io::IOStream, sif::SolverInformationFile)
    Base.write(io, "Header\n")
    Base.write(io, "  CHECK KEYWORDS \"Warn\"\n")
    Base.write(io, "  Mesh DB \".\" \"", sif.mesh_db, "\"\n")
    Base.write(io, "  Results Directory \"", sif.result_dir, "\"\n")
    Base.write(io, "End\n\n")
end

function write_simulation(io::IOStream, sif::SolverInformationFile)
    sim = sif.simulation

    Base.write(io, "Simulation\n")
    Base.write(io, "  ", format_property("Max Output Level", sim.output_level), "\n")
    Base.write(io, "  ", format_property("Coordinate System", sim.coordinate_system), "\n")
    Base.write(io, "  ", format_property("Simulation Type", sim.simulation_type), "\n")
    if (~isapprox(sim.coordinate_scaling, 1.0))
        Base.write(io, "  ", format_property("Coordinate Scaling", sim.coordinate_scaling), "\n")
    end
    for (key, val) ∈ sim.data
        Base.write(io, "  ", format_property(key, val), "\n")
    end
    Base.write(io, "End\n\n")
end

function write_constants(io::IOStream, sif::SolverInformationFile)
    constants = sort(collect(sif.constants), by=first)

    Base.write(io, "Constants\n")
    for (name, value) ∈ constants
        Base.write(io, "  ", format_property(name, value), "\n")
    end
    Base.write(io, "End\n\n")
end

function write_includes(io::IOStream, sif::SolverInformationFile)
    for include ∈ sif.includes
        Base.write(io, "Include \"$include\"")
    end
    Base.write(io, "\n\n")
end

function write_equations(io::IOStream, sif::SolverInformationFile)
    for equation ∈ sif.equations
        solver_names = [sif.solvers[s].name for s ∈ equation.solvers]

        Base.write(io, "! ", equation.name, "\n")
        Base.write(io, "Equation ", string(equation.id), "\n")
        Base.write(io, "  ", format_property("Active Solvers", equation.solvers), "   ! ", join(solver_names, ", "), "\n")
        Base.write(io, "End\n\n")
    end
end

function write_solvers(io::IOStream, sif::SolverInformationFile)
    for solver ∈ sif.solvers
        Base.write(io, "! ", solver.name, "\n")
        Base.write(io, "Solver ", string(solver.id), "\n")
        Base.write(io, "  ", format_property("Equation", solver.name), "\n")
        Base.write(io, "  ", format_property("Exec Solver", solver.exec), "\n")
        for (key, val) ∈ solver.data
            Base.write(io, "  ", format_property(key, val), "\n")
        end
        Base.write(io, "End\n\n")
    end
end

function write_materials(io::IOStream, sif::SolverInformationFile)
    for material ∈ sif.materials
        # Base.write(io, "! ", material.name, "\n")
        Base.write(io, "Material ", string(material.id), "\n")
        Base.write(io, "  ", format_property("Name", material.name), "\n")
        for (key, val) ∈ material.data
            Base.write(io, "  ", format_property(key, val), "\n")
        end
        Base.write(io, "End\n\n")
    end
end

function write_bodies(io::IOStream, sif::SolverInformationFile)
    for body ∈ sif.bodies
        # Base.write(io, "! ", body.name, "\n")
        Base.write(io, "Body ", string(body.id), "\n")
        Base.write(io, "  ", format_property("Name", body.name), "\n")
        if (length(body.target_bodies) > 0)
            Base.write(io, "  ", format_property("Target Bodies", body.target_bodies), "\n")
        end
        if (~ismissing(body.equation))
            name = sif.equations[body.equation].name
            Base.write(io, "  ", format_property("Equation", body.equation), "  ! ", name, "\n")
        end
        if (~ismissing(body.material))
            name = sif.materials[body.material].name
            Base.write(io, "  ", format_property("Material", body.material), "  ! ", name, "\n")
        end
        if (~ismissing(body.body_force))
            name = sif.body_forces[body.body_force].name
            Base.write(io, "  ", format_property("Body Force", body.body_force), "  ! ", name, "\n")
        end
        if (~ismissing(body.initial_condition))
            name = sif.initial_conditions[body.initial_condition].name
            Base.write(io, "  ", format_property("Initial Condition", body.initial_condition), "  ! ", name, "\n")
        end
        for (key, val) ∈ body.data
            Base.write(io, "  ", format_property(key, val), "\n")
        end
        Base.write(io, "End\n\n")
    end
end

function write_body_forces(io::IOStream, sif::SolverInformationFile)
    for body_force ∈ sif.body_forces
        # Base.write(io, "! ", body_force.name, "\n")
        Base.write(io, "Body Force ", string(body_force.id), "\n")
        Base.write(io, "  ", format_property("Name", body_force.name), "\n")
        for (key, val) ∈ body_force.data
            Base.write(io, "  ", format_property(key, val), "\n")
        end
        Base.write(io, "End\n\n")
    end
end

function write_initial_conditions(io::IOStream, sif::SolverInformationFile)
    for initial_condition ∈ sif.initial_conditions
        # Base.write(io, "! ", initial_condition.name, "\n")
        Base.write(io, "Initial Condition ", string(initial_condition.id), "\n")
        Base.write(io, "  ", format_property("Name", initial_condition.name), "\n")
        for (key, val) ∈ initial_condition.data
            Base.write(io, "  ", format_property(key, val), "\n")
        end
        Base.write(io, "End\n\n")
    end
end

function write_boundary_conditions(io::IOStream, sif::SolverInformationFile)
    for boundary_condition ∈ sif.boundary_conditions
        # Base.write(io, "! ", boundary_condition.name, "\n")
        Base.write(io, "Boundary Condition ", string(boundary_condition.id), "\n")
        Base.write(io, "  ", format_property("Name", boundary_condition.name), "\n")
        if (length(boundary_condition.target_boundaries) > 0)
            Base.write(io, "  ", format_property("Target Boundaries", boundary_condition.target_boundaries), "\n")
        end
        for (key, val) ∈ boundary_condition.data
            Base.write(io, "  ", format_property(key, val), "\n")
        end
        Base.write(io, "End\n\n")
    end
end

function write_components(io::IOStream, sif::SolverInformationFile)
    for component ∈ sif.components
        # Base.write(io, "! ", component.name, "\n")
        Base.write(io, "Component ", string(component.id), "\n")
        Base.write(io, "  ", format_property("Name", component.name), "\n")
        if (length(component.master_bodies) > 0)
            body_names = [sif.bodies[b].name for b ∈ component.master_bodies]
            Base.write(io, "  ", format_property("Master Bodies", component.master_bodies), "  ! ", join(body_names, ", "), "\n")
        end
        if (length(component.master_boundaries) > 0)
            boundary_names = [sif.boundary_conditions[b].name for b ∈ component.master_boundaries]
            Base.write(io, "  ", format_property("Master Boundaries", component.master_boundaries), "  ! ", join(boundary_names, ", "), "\n")
        end
        for (key, val) ∈ component.data
            Base.write(io, "  ", format_property(key, val), "\n")
        end
        Base.write(io, "End\n\n")
    end
end

function write(sif::SolverInformationFile)
    ~isdir(sif.data_path) && mkdir(sif.data_path)

    open(joinpath(sif.data_path, sif.filename), "w") do file
        write_header(file, sif)
        write_includes(file, sif)

        Base.write(file, "!---------------------------------------------------------\n")
        Base.write(file, "! Simulation setup\n")
        Base.write(file, "!---------------------------------------------------------\n")
        write_simulation(file, sif)
        write_constants(file, sif)
        Base.write(file, "\n")

        Base.write(file, "!---------------------------------------------------------\n")
        Base.write(file, "! Solver configuration\n")
        Base.write(file, "!---------------------------------------------------------\n")
        write_equations(file, sif)
        write_solvers(file, sif)
        Base.write(file, "\n")

        Base.write(file, "!---------------------------------------------------------\n")
        Base.write(file, "! Materials\n")
        Base.write(file, "!---------------------------------------------------------\n")
        write_materials(file, sif)
        Base.write(file, "\n")

        Base.write(file, "!---------------------------------------------------------\n")
        Base.write(file, "! Bodies & boundary conditions\n")
        Base.write(file, "!---------------------------------------------------------\n")
        write_bodies(file, sif)
        write_body_forces(file, sif)
        write_initial_conditions(file, sif)
        write_boundary_conditions(file, sif)
        Base.write(file, "\n")

        Base.write(file, "!---------------------------------------------------------\n")
        Base.write(file, "! Components\n")
        Base.write(file, "!---------------------------------------------------------\n")
        write_components(file, sif)
    end
end

function write_startinfo(sif::SolverInformationFile)
    @assert isdir(sif.data_path) "STARTINFO can only be written after generating the SIF"

    open(joinpath(sif.data_path, "ELMERSOLVER_STARTINFO"), "w") do file
        Base.write(file, "$(sif.filename)\n")
    end
end

function run_elmer_solver(sim_dir::String; solver_path="ElmerSolver", log_file="elmersolver.log")
    logfile = joinpath(sim_dir, log_file)
    cmd = Cmd(`$solver_path`, dir=sim_dir)
    pipe = pipeline(cmd, stdout=logfile)
    run(pipe)
end

function run_elmer_solver(sif::SolverInformationFile; solver_path="ElmerSolver", log_file="elmersolver.log")
    logfile = joinpath(sif.data_path, log_file)
    cmd = Cmd(`$solver_path $(sif.filename)`, dir=sif.data_path)
    pipe = pipeline(cmd, stdout=logfile)
    run(pipe)
end

function get_default_constants()
    return Dict(
        "Permittivity of Vacuum" => ustrip(VacuumElectricPermittivity),
        "Permeability of Vacuum" => ustrip(VacuumMagneticPermeability),
        "Stefan Boltzmann" => ustrip(StefanBoltzmannConstant),
        "Boltzmann Constant" => ustrip(BoltzmannConstant),
        "Unit Charge" => ustrip(ElementaryCharge)
    )
end