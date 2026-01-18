module CircuitBuilder

import Dates: format, now

using Elmer
using LinearAlgebra: Diagonal
using OrderedCollections: OrderedDict

export Resistor, Inductor, Capacitor, VoltageSource, CurrentSource, ElmerComponent
include("Components.jl")

export Circuit
export create_circuits, insert_component!
include("Circuits.jl")

export write_circuits
include("Write.jl")

export CoilMassive, CoilStranded, CoilFoil
export TerminalClosed, TerminalOpen
export add_circuits!, add_coil_data!
include("InterfaceSIF.jl")

end # module CircuitBuilder