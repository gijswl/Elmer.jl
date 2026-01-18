"""
    abstract type AbstractComponent

Represents an electrical circuit component.
"""
abstract type AbstractComponent end

"""
    mutable struct Resistor
    Resistor(name, nodes, value)

Represents a resistor named `name` with a `value` in [Ω]. It is connected between two circuit `nodes`.

It is used to build the Elmer circuit's B-matrix.
"""
mutable struct Resistor <: AbstractComponent
    name::String
    nodes::NTuple{2, Int}
    value::Real
end

"""
    mutable struct Inductor
    Inductor(name, nodes, value)

Represents a inductor named `name` with a `value` in [H]. It is connected between two circuit `nodes`.

It is used to build the Elmer circuit's A-matrix.
"""
mutable struct Inductor <: AbstractComponent
    name::String
    nodes::NTuple{2, Int}
    value::Real
end

"""
    mutable struct Capacitor
    Capacitor(name, nodes, value)

Represents a capacitor named `name` with a `value` in [F]. It is connected between two circuit `nodes`.

It is used to build the Elmer circuit's A-matrix.
"""
mutable struct Capacitor <: AbstractComponent
    name::String
    nodes::NTuple{2, Int}
    value::Real
end

"""
    mutable struct VoltageSource
    VoltageSource(name, nodes, value)

Represents a voltage source named `name` with a `value` in [V]. It is connected between two circuit `nodes`.

It is used to build the Elmer circuit's B-matrix.
"""
mutable struct VoltageSource <: AbstractComponent
    name::String
    nodes::NTuple{2, Int}
    value::Union{Nothing, Real, Complex}

    function VoltageSource(name, nodes, value=nothing)
        new(name, nodes, value)
    end
end

"""
    mutable struct CurrentSource
    CurrentSource(name, nodes, value)

Represents a current source named `name` with a `value` in [A]. It is connected between two circuit `nodes`.
Current is defined as flowing _out_ of the positive terminal `nodes[1]`.

It is used to build the Elmer circuit's B-matrix.
"""
mutable struct CurrentSource <: AbstractComponent
    name::String
    nodes::NTuple{2, Int}
    value::Union{Nothing, Real, Complex}

    function CurrentSource(name, nodes, value=nothing)
        new(name, nodes, value)
    end
end

"""
    mutable struct ElmerComponent
    ElmerComponent(name, nodes, component_id, master_bodies)

Represents a massive 2D or 3D coil named `name` connected between two circuit `nodes`.
The component is associated with the Component (`component_id`) section in the SIF and the physical bodies in `master_bodies`.

It is used to build the Elmer circuit's B-matrix.
"""
mutable struct ElmerComponent <: AbstractComponent
    name::String
    nodes::NTuple{2, Int}
    component_id::Int
    master_bodies::Vector{Int}

    function ElmerComponent(name, nodes, component_id, master_bodies)
        new(name, nodes, component_id, master_bodies)
    end
end

const SourceComponent = Union{VoltageSource, CurrentSource}