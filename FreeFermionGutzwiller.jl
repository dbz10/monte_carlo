# Module containing functions and types relevant to the
# gutzwiller projected free fermion system.

# needed functions:
# computing determinant, computing inverse, computing ratio
# doing a swap operation
# yellowbook which is "inverse" of R_up and R_down: it gives the spin
# and fermion label # on each site
# so for example, if we had 4 sites with R_up = [1,3] and R_down = [2,4]
# then the log book would look like [(1,1), (1,-1), (2,1), (2,-1)]

# another example: let R_up = [2 5 1] R_down = [4 3 6]
# book = [(3,1), (1,1), (2,-1), (1,-1), (2,1), (3,-1)]
# so the format is (electron #, spin)
# should it be the other way around? maybe thats better
# or maybe even better is to write it as
# book = [(3,"Up"), (1,"Up"), (2,"Down"), (1,"Down"), (2,"Up"), (3,"Down")]

module FreeFermionGutzwiller

include("mcbase.jl")

using StatsBase: sample, sample!
using IterTools: product
using LightGraphs:
    SimpleGraph, edges, src, dst, rem_edge!, add_edge!, neighbors
using LinearAlgebra: dot

# abstract type definitions
abstract type FreeFermionGutzwillerPolicy <: AbstractPolicy end
abstract type SwapMove <: AbstractMove end

# concrete type defs
mutable struct GutzwillerChain <: AbstractChain
    basechain::MarkovChain
    GutzwillerChain() = new()
end
get_MarkovChain(c::GutzwillerChain) = c.basechain


struct SpinConfiguration
    sc::Array
end

struct SwapNeighborMove <: SwapMove
    sites::Tuple
end

struct SwapNeighborsPolicy <: FreeFermionGutzwillerPolicy
end


function get_SwapNeighborMove(bonds::SimpleGraph{Int64})::SwapNeighborMove
    # draw a random bond and return the two sites to exchange
    list_of_edges = collect(1:ne(bonds))
    edge = sample(list_of_edges)
    move_sites = (src(edge),dst(edge))
    return SwapNeighborMove(move_sites)
end

struct GutzwillerState <: AbstractState
     r_up::Array{Int64} # position of spin up electrons
     r_down::Array{Int64} # position of spin down electrons
     spin_config::SpinConfiguration # a map of the spin configuration onto the graph
     bonds::SimpleGraph{Int64} # which bonds have disaligned spins
end

# high level function
function do_move(chain::GutzwillerChain)::GutzwillerChain
    # propose a move
    move = get_move(chain)
    # propose the move and decide whether to accept
    attempt_update_state!(chain,move)

    return chain
end

function get_move(chain::GutzwillerChain)
    policy = get_Policy(chain)
    move = get_move_from_policy(chain,policy)
    return move
end

function get_move_from_policy(chain::GutzwillerChain,policy::SwapNeighborsPolicy)::SwapNeighborMove
    bonds = get_State(chain).bonds
    move = get_SwapNeighborMove(bonds)
    return move
end

function attempt_update_state!(chain::GutzwillerChain,move)
    state = get_State(chain)


function get_init_state(chain::GutzwillerChain)::GutzwillerState
    # get all the model information from the chain object
    model = get_Model(chain)
    lattice = model["lattice"]
    filling = model["filling"]
    hamiltonian = model["hamiltonian"]
    dims = model["dims"]

    # Distribute electrons over the sites of the lattice
    num_sites = prod(dims)
    sites = collect(1:num_sites)

    # for now just assume filling = 1 and paramagnetic state
    # later i can add functionality for magnetic states
    # or states with holes
    n_up = Int(num_sites*filling/2)
    n_down = Int(num_sites*filling/2)

    # vectors to store the location of the up and down electrons
    R_up = zeros(Int64, n_up)
    R_down = zeros(Int64, n_down)


    sample!(sites,R_up,replace=false)
    leftover_sites = collect(setdiff(Set(sites),Set(R_up)))
    sample!(leftover_sites,R_down,replace=false)

    spin_configuration = get_spin_config(R_up,R_down)

    # create list of neighbors which can be swapped in current configuration
    bonds = get_init_bonds(lattice.graph,spin_configuration)

    # make the state object and return
    state = GutzwillerState(R_up,R_down,spin_configuration,bonds)
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
