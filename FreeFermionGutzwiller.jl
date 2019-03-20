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
include("FFGmath.jl")
include("FFGhelpers.jl")

using StatsBase: sample, sample!
using IterTools: nth
using LightGraphs:
    SimpleGraph, edges, src, dst, rem_edge!, add_edge!, neighbors,
    adjacency_matrix, ne
using LinearAlgebra: dot, det, inv, eigen

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
    """ draw a random bond and return the two sites attached to that bond """
    index = rand(1:ne(bonds))
    edge = nth(edges(bonds),index)
    move_sites = (src(edge),dst(edge))
    return SwapNeighborMove(move_sites)
end

struct GutzwillerState <: AbstractState
     r_up::Array{Int64} # position of spin up electrons
     r_down::Array{Int64} # position of spin down electrons
     spin_config::SpinConfiguration # a map of the spin configuration onto the graph
     bonds::SimpleGraph{Int64} # which bonds have disaligned spins
     wavefunctions::Array # wavefunctions ϕ_j(r_i) where ϵ_j < ϵ_fermi
     det_A_up::Float64
     det_A_down::Float64
     A_inv_up::Array
     A_inv_down::Array
end

# high level function


function get_Move(chain::GutzwillerChain)
    """ Get a move from the policy of the chain """
    policy = get_Policy(chain)
    move = get_move_from_policy(chain,policy)
    return move
end

function get_move_from_policy(chain::GutzwillerChain,policy::SwapNeighborsPolicy)::SwapNeighborMove
    """ Specific function for gutzwiller chain with swap neighbor policy"""
    bonds = get_State(chain).bonds
    move = get_SwapNeighborMove(bonds)
    return move
end



function update_state!(chain::GutzwillerChain,move)
end


function compute_ratio(chain::GutzwillerChain,move::SwapNeighborMove)

    """ Computes transition amplitude"""
    state = get_State(chain)
    wavefunctions = state.wavefunctions

    u, v = get_update_vectors(wavefunctions,move)
    # note that u to update the down determinant is -u
    # that was used to update the up determinant since
    # u_i(r1-> r2) = ϕ_i(r2) - ϕi(r1)
    # hence u_i(r2-> r1) = -u_i(r1-> r2)

    # up and down electrons have same wavefunction.
    # we do need to figure out whether move[1] is
    # an up or down spin though
    spin_A = get_spin_of_site(state,move[1])

    if spin_A == +1
        det_ratio_up = det_ratio_factor(state.A_inv_up,u,v)
        det_ratio_down = det_ratio_factor(state.A_inv_down,-u,v)
    else
        det_ratio_up = det_ratio_factor(state.A_inv_up,-u,v)
        det_ratio_down = det_ratio_factor(state.A_inv_down,u,v)
    end

    # we also need the ratio of proposal factors
    pf_ratio = get_proposal_factor_ratio(chain,move)

    return ratio_up*ratio_down*pf_ratio
end

function get_proposal_factor_ratio(chain::GutzwillerChain,move::SwapNeighborMove)
    state = get_State(chain)
    bonds = state.bonds
    configuration_factor_old = count_bonds(bonds)
    new_bonds = get_updated_bonds(state,move)
    configuration_factor_new = count_bonds(new_bonds)
    return configuration_factor_new/configuration_factor_old
end


count_bonds(graph::SimpleGraph{Int64}) = ne(graph)

function get_updated_bonds(chain::GutzwillerChain,move::SwapNeighborMove)::SimpleGraph{Int64}
    """ Returns a new graph object with the bonds of the updated state"""
    state = get_State(chain)
    model = get_Model(chain)
    lattice = model["lattice"]

    sc = deepcopy(state.spin_config.sc)
    bonds = deepcopy(state.bonds)
    sc[move.sites[1]], sc[move.sites[2]] = sc[move.sites[2]], sc[move.sites[1]]

    for site in move.sites
        for neighbor in neighbors(lattice.graph,site)
            if sc[site] == sc[neighbor]
                rem_edge!(bonds,site,neighbor)
            else
                add_edge!(bonds,site,neighbor)
            end
        end
    end
    return bonds
end



function get_init_state(chain::GutzwillerChain)::GutzwillerState
    # get all the model information from the chain object
    model = get_Model(chain)
    lattice = model["lattice"]
    filling = model["filling"]
    fermi_energy = model["fermi_energy"]
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

    # get occupied wavefunctions from lattice and fermi energy
    filled_states = get_Wavefunctions(model)

    # calculate determinants
    det_A_up = get_Determinant(R_up,filled_states)
    det_A_down = get_Determinant(R_down,filled_states)

    # calculate inverse matrices for A_up and A_down
    A_inv_up = get_Inverse_Matrix(R_up,filled_states)
    A_inv_down = get_Inverse_Matrix(R_down,filled_states)

    # make the state object and return
    state = GutzwillerState(R_up,R_down,spin_configuration,bonds, filled_states,
                            det_A_up, det_A_down, A_inv_up,A_inv_down)
    return state
end

function get_Wavefunctions(model::Dict)::Array
    lattice = model["lattice"]
    fermi_energy = model["fermi_energy"]
    hamiltonian = model["hamiltonian"]
    # calculate the determinants of the up and down spin matrices
    EO = hamiltonian(lattice) # eigen object
    wavefunctions = EO.vectors # wavefunctions[i,j] =   ϕ_j(r_i)
    eigenvalues = EO.values
    filled_mask = eigenvalues.<fermi_energy # find index of filled wavefunctions
    filled_states = wavefunctions[:,filled_mask] # ϕ_{ϵ(j) < ϵ_f} (r_i)
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

function get_Determinant(r,states)
    """ Form the matrix ϕ(r_i, k_j) and calculate its determinant
    this is not entirely general, could stand to make it
    more modular. different models could have different wavefunctions """
    mat = states[r,:]
    return det(mat)
end

function get_Inverse_Matrix(r,states)
    """ Compute the inverse of the matrix ϕ(r_i, k_j) """
    mat = states[r,:]
    return inv(mat)
end




end
