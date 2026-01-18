"""
    mutable struct Circuit

A circuit holds a list of `components` and reference node `ref_node`.
Elmer supports the definition of multiple circuits, numbered through `index`.
"""
mutable struct Circuit
    index::Int
    components::Vector{AbstractComponent}
    ref_node::Int

    function Circuit(index; components=AbstractComponent[], ref_node=1)
        new(index, components, ref_node)
    end
end

"""
    create_circuits(n)

Creates an array of empty `Circuit`s with indices from 1 to `n`.
"""
create_circuits(n::Int) = [Circuit(i) for i ∈ 1:n]

function get_nodes(c::Circuit)
    nodes = Int[]
    for component ∈ c.components
        node1 = component.nodes[1]
        node2 = component.nodes[2]

        node1 ∈ nodes || push!(nodes, node1)
        node2 ∈ nodes || push!(nodes, node2)
    end

    return nodes
end

"""
    get_nnodes(c)

Returns the number of circuit nodes in circuit `c`.
"""
get_nnodes(c::Circuit) = length(get_nodes(c))

"""
    get_nedges(c)

Returns the number of edges (= number of components) in circuit `c`.
"""
get_nedges(c::Circuit) = length(c.components)

"""
    get_components(c)

Returns the components of circuit `c`.
"""
get_components(c::Circuit) = c.components

"""
    insert_component!(c, component)
    insert_component!(c, components)

Adds a single `component` or array of `components` to circuit `c`.
"""
@inline function insert_component!(c::Circuit, component::AbstractComponent)
    push!(c.components, component)
end

function insert_component!(c::Circuit, components::Vector{AbstractComponent})
    for component ∈ components
        insert_component!(c, component)
    end
end

"""
    get_incidence_matrix(c)

Returns the incidence matrix of circuit `c`.

The matrix is constructed using nodes to represent rows, and edges columns.
    See: https://en.wikipedia.org/wiki/Incidence_matrix
    The directed graph implies that the direction of each edge is given by the nodes, 
    for which you have a positive and a negative terminal for each component.
"""
function get_incidence_matrix(c::Circuit)
    A = zeros(get_nnodes(c), get_nedges(c))
    for (edge, component) ∈ enumerate(get_components(c))
        A[component.nodes[1], edge] = 1
        A[component.nodes[2], edge] = -1
    end

    # Remove the reference node row
    node_idx = 1:get_nnodes(c) .!= c.ref_node
    return A[node_idx, :]
end

"""
    add_resistance_entry!(R, i, component)

Defines the contribution of each `component` type to the resistance matrix `R`.
"""
add_resistance_entry!(::Vector{<:Real}, ::Int, ::AbstractComponent) = nothing
add_resistance_entry!(R::Vector{<:Real}, i::Int, component::Resistor) = (R[i] = component.value)
add_resistance_entry!(R::Vector{<:Real}, i::Int, ::CurrentSource) = (R[i] = 1)
add_resistance_entry!(R::Vector{<:Real}, i::Int, ::Capacitor) = (R[i] = 1)

"""
    get_resistance_matrix(c)

Returns the resistance matrix of circuit `c`. 
The matrix contains contributions from `Resistor`, `CurrentSource`, and `Capacitor`.
"""
function get_resistance_matrix(c::Circuit)
    R = zeros(get_nedges(c))
    for (edge, component) ∈ enumerate(get_components(c))
        add_resistance_entry!(R, edge, component)
    end

    return Diagonal(R)
end

"""
    add_conductance_entry!(G, i, component)

Defines the contribution of each `component` type to the conductance matrix `G`.
"""
add_conductance_entry!(::Vector{<:Real}, ::Int, ::AbstractComponent) = nothing
add_conductance_entry!(G::Vector{<:Real}, i::Int, ::Resistor) = (G[i] = -1)
add_conductance_entry!(G::Vector{<:Real}, i::Int, ::VoltageSource) = (G[i] = 1)
add_conductance_entry!(G::Vector{<:Real}, i::Int, ::Inductor) = (G[i] = 1)

"""
    get_conductance_matrix(c)

Returns the conductance matrix of circuit `c`. 
The matrix contains contributions from `Resistor`, `VoltageSource`, and `Inductor`.
"""
function get_conductance_matrix(c::Circuit)
    G = zeros(get_nedges(c))
    for (edge, component) ∈ enumerate(get_components(c))
        add_conductance_entry!(G, edge, component)
    end

    return Diagonal(G)
end

"""
    add_inductance_entry!(L, i, component)

Defines the contribution of each `component` type to the inductance matrix `L`.
"""
add_inductance_entry!(::Vector{<:Real}, ::Int, ::AbstractComponent) = nothing
add_inductance_entry!(L::Vector{<:Real}, i::Int, component::Inductor) = (L[i] = -component.value)

"""
    get_inductance_matrix(c)

Returns the inductance matrix of circuit `c`. 
The matrix contains contributions from `Inductor`.
"""
function get_inductance_matrix(c::Circuit)
    L = zeros(get_nedges(c))
    for (edge, component) ∈ enumerate(get_components(c))
        add_inductance_entry!(L, edge, component)
    end

    return Diagonal(L)
end

"""
    add_capacitance_entry!(C, i, component)

Defines the contribution of each `component` type to the capacitance matrix `L`.
"""
add_capacitance_entry!(::Vector{<:Real}, ::Int, ::AbstractComponent) = nothing
add_capacitance_entry!(C::Vector{<:Real}, i::Int, component::Capacitor) = (C[i] = -component.value)

"""
    get_capacitance_matrix(c)

Returns the capacitance matrix of circuit `c`. 
The matrix contains contributions from `Capacitor`.
"""
function get_capacitance_matrix(c::Circuit)
    C = zeros(get_nedges(c))
    for (edge, component) ∈ enumerate(get_components(c))
        add_capacitance_entry!(C, edge, component)
    end

    return Diagonal(C)
end

"""
    add_rhs_entry!(rhs, i, component)

Defines the contribution of each `component` type to the right-hand side vector `rhs`.
"""
add_rhs_entry!(::Vector, ::Int, ::AbstractComponent) = nothing
add_rhs_entry!(rhs::Vector, i::Int, component::CurrentSource) = (rhs[i] = component.name)
add_rhs_entry!(rhs::Vector, i::Int, component::VoltageSource) = (rhs[i] = component.name)

"""
    get_rhs(c)

Returns the right-hand side vector of circuit `c`. 
The vector contains contributions from `VoltageSource` and `CurrentSource`.
"""
function get_rhs(c::Circuit)
    rhs = ["" for _ ∈ 1:get_nedges(c)]
    for (edge, component) ∈ enumerate(get_components(c))
        add_rhs_entry!(rhs, edge, component)
    end

    return rhs
end

"""
    swap_rows!(v, i, j)
    swap_rows!(m, i, j)

Helper functions to swap rows `i` and `j` of a vector `v` or matrix `m`
"""
function swap_rows!(v::Vector, i::Int, j::Int)
    v[i], v[j] = v[j], v[i]
end
swap_rows!(m::Matrix, i::Int, j::Int) = Base.swaprows!(m, i, j)

"""
    swap_rows!(Bmat, Amat, svec, c)

Ensures that voltage components rows are zero rows by systematically swapping rows.
In order to couple lumped circuit networks to Elmer, the voltage rows for Elmer Components
    need to be empty. This way the matrix is completed by using the Component keyword in the .sif file.
"""
function swap_rows!(Bmat::Matrix, Amat::Matrix, svec::Vector, c::Circuit)
    Z = Bmat .+ Amat
    idx = axes(Z, 1)
    zero_rows = [all(Z[row, :] .== 0) for row ∈ idx]
    zero_rows = idx[zero_rows]

    components = get_components(c)
    idx = axes(components, 1)
    vcomp_rows = [typeof(component) <: ElmerComponent for component ∈ components]
    vcomp_rows = idx[vcomp_rows] .+ get_nedges(c)

    for (zrow, vrow) ∈ zip(zero_rows, vcomp_rows)
        swap_rows!(Bmat, zrow, vrow)
        swap_rows!(Amat, zrow, vrow)
        swap_rows!(svec, zrow, vrow)
    end

    return Bmat, Amat, svec
end

"""
    fix_kvl_sign!(Bmat, c)

Inverts the sign of KVL entries corresponding to a current or voltage source.
"""
function fix_kvl_sign!(Bmat::Matrix, c::Circuit)
    components = get_components(c)
    idx = axes(components, 1)
    source_rows = [typeof(component) <: SourceComponent for component ∈ components]
    source_rows = idx[source_rows]

    for srow ∈ source_rows
        idx1 = srow + get_nnodes(c) - 1
        idx2 = srow + get_nedges(c)
        Bmat[idx1, idx2] *= -1
    end
end

"""
    get_tableau_matrix(c)

Returns the Elmer circuit matrices `A`, `B`, and vector `source` according to the Sparse Tableau Method.
"""
function get_tableau_matrix(c::Circuit)
    Amat = get_incidence_matrix(c)
    Rmat = get_resistance_matrix(c)
    Gmat = get_conductance_matrix(c)
    Lmat = get_inductance_matrix(c)
    Cmat = get_capacitance_matrix(c)
    fvec = get_rhs(c)

    nnodes = get_nnodes(c)
    nedges = get_nedges(c)
    mnodes = nnodes - 1 # Number of nodes minus the reference node

    M_kcl = hcat(Amat, zeros(mnodes, nedges), zeros(mnodes, mnodes))
    M_kvl = hcat(zeros(nedges, nedges), -Diagonal(ones(nedges)), transpose(Amat))
    M_comp = hcat(Rmat, Gmat, zeros(nedges, mnodes))

    ElmerB = vcat(M_kcl, M_kvl, M_comp)
    ElmerA = vcat(zeros(mnodes + nedges, mnodes + 2 * nedges), hcat(Lmat, Cmat, zeros(nedges, mnodes)))

    src_zero = ["" for _ ∈ 1:(mnodes+nedges)]
    ElmerSrc = vcat(src_zero, fvec)

    fix_kvl_sign!(ElmerB, c)
    swap_rows!(ElmerB, ElmerA, ElmerSrc, c)

    return ElmerB, ElmerA, ElmerSrc
end

"""
    get_unknown_names(c)

Returns a list of names of the unknowns (degrees of freedom).
"""
function get_unknown_names(c::Circuit)
    unknowns = String[]
    nodes = setdiff(get_nodes(c), c.ref_node)

    for component ∈ get_components(c)
        if (typeof(component) <: ElmerComponent)
            push!(unknowns, "i_component($(component.component_id))")
        else
            push!(unknowns, "i_$(component.name)")
        end
    end

    for component ∈ get_components(c)
        if (typeof(component) <: ElmerComponent)
            push!(unknowns, "v_component($(component.component_id))")
        else
            push!(unknowns, "v_$(component.name)")
        end
    end

    for node ∈ nodes
        push!(unknowns, "u_$(node)_circuit_$(c.index)")
    end

    return unknowns
end