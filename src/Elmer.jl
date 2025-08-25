module Elmer

include("Utils.jl")

export SolverInformationFile
export Simulation, Solver, Equation, Material, Body, BodyForce, InitialCondition, BoundaryCondition, Component
export add_constant!, add_solver!, add_equation!, add_material!, add_body!, add_body_force!, add_initial_condition!, add_boundary_condition!, add_component!
export load_solver!
include("ElmerSolver.jl")
include("ElmerGrid.jl")

end # module Elmer
