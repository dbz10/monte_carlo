# Module containing functions and types relevant to the
# gutzwiller projected free fermion system.

# needed functions:
#  computing ratio
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
using IterTools: nth
using LightGraphs:
    SimpleGraph, edges, src, dst, rem_edge!, add_edge!, neighbors,
    adjacency_matrix, ne, nv
using LinearAlgebra: dot, det, inv, eigen, cond

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
    """ This is borderline pointless, it's just a typed array"""
    sc::Array
end

struct SwapNeighborMove <: SwapMove
    sites::Tuple
end

struct SwapNeighborsPolicy <: FreeFermionGutzwillerPolicy
end

mutable struct GutzwillerState <: AbstractState
     r_up::Array{Int64} # position of spin up electrons
     r_down::Array{Int64} # position of spin down electrons
     business_directory::Array{Int64} # map from site to electron index.
     spin_config::SpinConfiguration # a map of the spin configuration onto the graph
     bonds::SimpleGraph{Int64} # which bonds have disaligned spins
     wavefunctions::Array # wavefunctions ϕ_j(r_i) where ϵ_j < ϵ_fermi
     det_A_up::Float64
     det_A_down::Float64
     A_inv_up::Array
     A_inv_down::Array
end

# Load in helper functions
include("FFGmath.jl")
include("FFGhelpers.jl")

""" High level functions """
function get_Move(chain::GutzwillerChain)
    """ Get a move from the policy of the chain """
    policy = get_Policy(chain)
    move = get_move_from_policy(chain,policy)
    return move
end

function update_state!(chain::GutzwillerChain,move)
    """ Annoyingly, the order in which things are updated is important.
    In particular, update_Rs! needs the current business directory and
    update_Bonds! needs the current spin config"""
    state = get_State(chain)
    update_Rs!(state,move)
    update_Bonds!(chain,move)
    update_Determinants!(chain,move)
    update_Inverses!(state,move)
    update_Spin_config!(state,move)
    update_Business_directory!(state,move) # important that BD is updated last
end

function update_Rs!(state::GutzwillerState,move::SwapNeighborMove)
    r_up = state.r_up
    r_down = state.r_down
    bd = state.business_directory
    # use bd to find what electrons are being exchanged
    l1, l2 = move.sites[1], move.sites[2]
    electron_A = bd[l1]
    electron_B = bd[l2]

    if electron_B > electron_A # then A is spin up and B is spin down
        electron_B = get_electron_index(electron_B,state)
        r_up[electron_A] = l2
        r_down[electron_B] = l1
    else # A is down and B is up
        electron_A = get_electron_index(electron_A,state)
        r_down[electron_A] = l2
        r_up[electron_B] = l1
    end
end

function update_Spin_config!(state::GutzwillerState,move::SwapNeighborMove)
    sc = state.spin_config.sc
    l1 = move.sites[1]
    l2 = move.sites[2]
    sc[l1], sc[l2] = sc[l2], sc[l1]
end

function update_Bonds!(chain::GutzwillerChain,move::SwapNeighborMove)
    state = get_State(chain)
    tmp = get_updated_bonds(chain,move)
    bonds = state.bonds
    bonds = tmp
end

function update_Determinants!(chain::GutzwillerChain,move::SwapNeighborMove)
    state = get_State(chain)
    trash, det_ratio_up, det_ratio_down = compute_ratio(chain,move)
    state.det_A_up *= det_ratio_up
    state.det_A_down *= det_ratio_down
end

function update_Inverses!(state::GutzwillerState,move::SwapNeighborMove)
    u1,v1,u2,v2 = get_update_vectors(state,move)
    # recall u1 v1 are for up, u2 v2 are for down
    t = sherman_morrison_inverse_update(state.A_inv_up,u1,v1)
    r = sherman_morrison_inverse_update(state.A_inv_down,u2,v2)

    state.A_inv_up = t
    state.A_inv_down = r
end

function update_Business_directory!(state::GutzwillerState,move::SwapNeighborMove)
    bd = state.business_directory
    l1 = move.sites[1]
    l2 = move.sites[2]
    bd[l1], bd[l2] = bd[l2], bd[l1]
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

    # find a well conditioned initial state
    conditioning_tol = 10^6

    R_up, R_down, filled_states = get_conditioned_state(
        n_up,n_down,
        model,conditioning_tol)

    # make the business directory
    business_directory = make_business_directory(R_up,R_down)

    # make the spin configuration
    spin_configuration = make_spin_config(R_up,R_down)

    # create list of neighbors which can be swapped in current configuration
    bonds = get_init_bonds(lattice.graph,spin_configuration)

    # calculate determinants
    det_A_up = get_Determinant(R_up,filled_states)
    det_A_down = get_Determinant(R_down,filled_states)

    # calculate inverse matrices for A_up and A_down
    A_inv_up = get_Inverse_Matrix(R_up,filled_states)
    A_inv_down = get_Inverse_Matrix(R_down,filled_states)

    # make the state object and return
    state = GutzwillerState(
        R_up,R_down, business_directory, spin_configuration,bonds,
        filled_states, det_A_up, det_A_down, A_inv_up,A_inv_down
        )
    return state
end

""" Lower level functions """



function get_SwapNeighborMove(bonds::SimpleGraph{Int64})::SwapNeighborMove
    """ draw a random bond and return the two sites attached to that bond """
    index = rand(1:ne(bonds))
    edge = nth(edges(bonds),index) #nth is from IterTools
    move_sites = (src(edge),dst(edge))
    return SwapNeighborMove(move_sites)
end

function get_move_from_policy(chain::GutzwillerChain,policy::SwapNeighborsPolicy)::SwapNeighborMove
    """ Specific function for gutzwiller chain with swap neighbor policy"""
    bonds = get_State(chain).bonds
    move = get_SwapNeighborMove(bonds)
    return move
end

function compute_ratio(chain::GutzwillerChain,move::SwapNeighborMove)
    """ Computes transition amplitude. Returns transition amplitude,
    det_ratio_up, det_ratio_down. First return argument must always be ratio"""
    state = get_State(chain)
    wavefunctions = state.wavefunctions

    # update vectors for up spin determinant (1) and down spin determinant (2)
    u1, v1, u2, v2 = get_update_vectors(state,move)

    det_ratio_up = det_ratio_factor(state.A_inv_up,u1,v1)
    det_ratio_down = det_ratio_factor(state.A_inv_down,u2,v2)

    # we also need the ratio of proposal factors
    pf_ratio = get_proposal_factor_ratio(chain,move)

    return det_ratio_up*det_ratio_down*pf_ratio, det_ratio_up, det_ratio_down
end

function get_proposal_factor_ratio(chain::GutzwillerChain,move::SwapNeighborMove)
    state = get_State(chain)
    bonds = state.bonds
    configuration_factor_forwards = 1.0/count_bonds(bonds)
    new_bonds = get_updated_bonds(chain,move)
    configuration_factor_backwards = 1.0/count_bonds(new_bonds)
    return configuration_factor_backwards/configuration_factor_forwards
end

function make_spin_config(R_up,R_down)::SpinConfiguration
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

function make_business_directory(r_up::Array,r_down::Array)::Array
    """ Makes the directory for the lattice. The i'th element
    is the index of the electron living on site i. indices <= N/2
    refer to up spins and indices > N/2 refer to down spins """
    r_full = [r_up ; r_down]
    sortkey = sortperm(r_full)
    return sortkey
end

end
