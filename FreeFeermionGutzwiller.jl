module FreeFermionGutzwiller

include("mcbase.jl")

using StatsBase: sample!
using IterTools: product
using LightGraphs:
    SimpleGraph, edges, src, dst, rem_edge!, add_edge!, neighbors

# abstract type definitions
abstract type FreeFermionGutzwillerPolicy <: AbstractPolicy end
abstract type FreeFermionGutzwillerMove <: AbstractMove end

# concrete type defs
mutable struct GutzwillerChain <: AbstractChain
    basechain::MarkovChain
    Gutzwiller_Chain() = new()
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

    # create list of neighbors which can be swapped in current configuration
    bonds = get_bonds(lattice.graph,R_up,R_down)

    # make the state object and return
    state = GutzwillerState(R_up,R_down,bonds)
    return state
end


function get_filled_k_states(lattice,fermi_energy,hamiltonian)::Tuple
# for now let's limit to lattices with a single site unit cell.
# we need the dispersion relation, which is
# sum of lattice vectors j e^{i k r_j} + h.c.
# if ϵ(k) < ϵ_f, then we should include that k value in the list
# of filled k states

lattice_vectors = lattice.lattice_vectors

all_k_states = get_all_k_states(lattice)
all_k_energies = tight_binding_dispersion(lattice_vectors,all_k_states)

filled_k_states = all_k_states(all_k_energies < fermi_energy)

return filled_k_states
end

function get_all_k_states(lattice)::Tuple
    lattice_dimensions = lattice.dims
    single_dimension_k_lists = [2*pi/d * collect(1:d) for d in dims]

    all_k_states = product([k_list for k_list in single_dimension_k_lists])
    return all_k_states
end

function tight_binding_dispersion(lattice_vectors::Array{Int64,3}, k::Tuple)::Float64
    k_vector = collect(k)
    terms = [exp(1im*dot(k_vector,lattice_vectors[i,:,:])) for i in 1:length(k)]
end


end
