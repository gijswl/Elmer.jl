abstract type AbstractCoilType end
struct CoilMassive <: AbstractCoilType end

struct CoilStranded <: AbstractCoilType
    number_turns::Real
    resistance::Union{Nothing, Real}
end
CoilStranded(number_turns::Real) = CoilStranded(number_turns, nothing)

struct CoilFoil <: AbstractCoilType
    number_turns::Real
    thickness::Real
end

abstract type AbstractTerminalType end
struct TerminalClosed <: AbstractTerminalType
    area::Real
end
struct TerminalOpen <: AbstractTerminalType
    boundaries::NTuple{2, Int}
    symmetry::Real
end
TerminalOpen(bnd1::Int, bnd2::Int, symmetry::Real) = TerminalOpen((bnd1, bnd2), symmetry)

function add_circuits!(sif::SolverInformationFile, circuits::Vector{Circuit})
    components = Dict()
    body_forces = Dict()

    for c ∈ circuits
        bf_data = OrderedDict()

        for component ∈ get_components(c)
            if (typeof(component) <: ElmerComponent)
                id = add_component!(sif, component.name; master_bodies=component.master_bodies)
                @assert id == component.component_id "Circuit component numbering should be monotonically increasing."
                components[id] = component
            elseif (typeof(component) <: SourceComponent)
                entry = format_source(component.name)

                if (typeof(component.value) <: Complex)
                    bf_data["$entry re"] = real(component.value)
                    bf_data["$entry im"] = imag(component.value)
                elseif (typeof(component.value) <: Real)
                    bf_data[entry] = component.value
                end
            end
        end

        body_forces[c.index] = add_body_force!(sif, "Circuit $(c.index)"; data=bf_data)
    end

    return components, body_forces
end

function add_coil_data!(sif, id, coil::AbstractCoilType; symmetry::Real=1.0)
    data = OrderedDict()

    coil_data!(data, coil)

    data["Symmetry Coefficient"] = symmetry

    update_component_data!(sif, id, data)
end

function add_coil_data!(sif, id, coil::AbstractCoilType, terminal::AbstractTerminalType)
    data = OrderedDict()

    coil_data!(data, coil)
    terminal_data!(data, terminal)

    update_component_data!(sif, id, data)
end

function coil_data!(data, coil::CoilMassive)
    data["Coil Type"] = "Massive"
end

function coil_data!(data, coil::CoilStranded)
    data["Coil Type"] = "Stranded"
    data["Number of Turns"] = Float64(coil.number_turns)
    if (typeof(coil.resistance) <: Real)
        data["Resistance"] = coil.resistance
    end
end

function coil_data!(data, coil::CoilFoil)
    data["Coil Type"] = "Foil Winding"
    data["Coil Thickness"] = coil.thickness
    data["Number of Turns"] = Float64(coil.number_turns)
end

function terminal_data!(data, terminal::TerminalClosed)
    data["Coil Use W Vector"] = true
    data["W Vector Variable Name"] = "String \"CoilCurrent e\""
    data["Electrode Area"] = terminal.area
end

function terminal_data!(data, terminal::TerminalOpen)
    data["Coil Use W Vector"] = true
    data["W Vector Variable Name"] = "String \"CoilCurrent e\""
    data["Electrode Boundaries"] = terminal.boundaries
    data["Circuit Equation Voltage Factor"] = terminal.symmetry
end