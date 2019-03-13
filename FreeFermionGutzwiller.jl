# Module containing functions and types relevant to the
# gutzwiller projected free fermion system.

module FreeFermionGutzwiller

include("mcbase.jl")

using StatsBase: sample!
using IterTools: product
using LightGraphs:
    SimpleGraph, edges, src, dst, rem_edge!, add_edge!, neighbors
using LinearAlgebra: dot

# abstract type definitions
abstract type FreeFermionGutzwillerPolicy <: AbstractPolicy end
abstract type FreeFermionGutzwillerMove <: AbstractMove end

# concrete type defs
mutable struct GutzwillerChain <: AbstractChain
    basechain::MarkovChain
    Gutzwiller_Chain() = new()
end

struct SpinConfiguration
    sc::Array
end

struct GutzwillerState <: AbstractState
     r_up::Tuple # position of spin up electrons
     r_down::Tuple # position of spin down electrons
     bonds::Tuple # which bonds have disaligned spins
end

function get_init_state(chain::GutzwillerChain)::GutzwillerState
    # get all the model information from the chain object
    model = get_model(chain)
    lattice = model["lattice"]
    filling = model["filling"]
    hamiltonian = model["hamiltonian"]
    dims = model["dims"]

    # Distribute electrons over the sites of the lattice
    num_sites = prod(dims)
    sites = collect(1:num_sites)

    # for now just assume fillign = 1 and paramagnetic state
    # later i can add functionality for magnetic states
    # or states with holes
    n_up = Int(num_sites*filling/2)
    n_down = Int(num_sites*filling/2)

    # vectors to store the location of the up and down electrons
    R_up = zeros(Int64, n_up)
    R_down = zeros(Int64, n_down)


    sample!(sites,R_up,replace=false)
    leftover_sites = collect(setdiff(Set(R_up),Set(sites)))
    sample!(leftover_sites,R_down,replace=false)

    spin_configuration = get_spin_config(R_up,R_down)

    # create list of neighbors which can be swapped in current configuration
    bonds = get_bonds(lattice.graph,spin_configuration)

    # make the state object and return
    state = GutzwillerState(R_up,R_down,bonds)
    return state
end

function get_spin_config(R_up,R_down)::SpinConfiguration
    nsites = length([R_up R_down])
    sc = zeros(Int64,nsites)
    sc[R_up] .= +1
    sc[R_down] .= -1

    return SpinConfiguration(sc)
end

function get_init_bonds(
        graph::SimpleGraph{Int64},
        spin_config::SpinConfiguration)::SimpleGraph{Int64}

    # when the graph is created, it is fully connected.
    # now we loop over the edges and use the spin configuration to remove bonds
    # where necessary. we do not need to add any bonds in the initial step

    # we need to make a deep copy of the graph because as we operate
    # on the graph its list of edges will change. instead we operate on a
    # copy of graph, and return that
    bonds = deepcopy(graph)

    for edge in deepcopy(edges(bonds))
        u, v = src(edge), dst(edge)
        if spin_config.sc[u] == spin_config.sc[v]
            rem_edge!(bonds,edge)
        end
    end

    return bonds
end



function get_filled_k_states(lattice,fermi_energy,dispersion)::Array{Tuple}
    # for now let's limit to lattices with a single site unit cell.
    # we need the dispersion relation, which is
    # sum of lattice vectors j e^{i k r_j} + h.c.
    # if ϵ(k) < ϵ_f, then we should include that k value in the list
    # of filled k states
    # dispersion is a function that calculates the energy of a particular
    # k state based on the lattice

    lattice_vectors = lattice.lattice_vectors

    all_k_states = get_all_k_states(lattice)
    all_k_energies = [dispersion(lattice_vectors,k) for k in all_k_states]

    filled_k_states = all_k_states[all_k_energies .< fermi_energy]

    return filled_k_states
end


function get_all_k_states(lattice)::Array{Tuple}
    lattice_dimensions = lattice.dims
    single_dimension_k_lists = [2*pi/d*(1:d) for d in lattice_dimensions]
    all_k_states_obj = product(single_dimension_k_lists...)
    list_of_k_states = reshape(collect(all_k_states_obj),(:,1))

    return list_of_k_states
end

function tight_binding_dispersion(lattice_vectors::Array{Int64}, k::Tuple)::Float64
    # using t = 1
    # dispersion relation for hopping -t \sum_{ij} c^\dagger_i c_j
    k_vector = collect(k)
    terms = [exp(1im*dot(k_vector,lattice_vectors[i,:,:])) for i in 1:length(k)]
    return -real(sum(terms+conj(terms)))
end


end
