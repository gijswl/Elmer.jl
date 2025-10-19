module Elmer

using OrderedCollections
using PhysicalConstants.CODATA2022
using Printf
using Unitful

export SolverInformationFile
export Simulation, Solver, Equation, Material, Body, BodyForce, InitialCondition, BoundaryCondition, Component
export add_constant!, add_include!, add_solver!, add_equation!, add_material!, add_body!, add_body_force!, add_initial_condition!, add_boundary_condition!, add_component!
export update_simulation_data!, update_solver_data!, update_material_data!, update_body_data!, update_bf_data!, update_ic_data!, update_bc_data!, update_component_data!
export load_solver!, load_material!
include("ElmerSolver.jl")
include("ElmerGrid.jl")

export load_dat
include("PostProcessing.jl")

end # module Elmer
