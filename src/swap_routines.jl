include("lattices.jl")

mutable struct ElectronTinder
    """ Database to match lonely up spins with lonely down spins by site"""
    lonely_up::Array{Int}
    lonely_down::Array{Int}
end

function addUser!(t::ElectronTinder,site::Int,spin::Int)
    spin==1 ? append!(t.lonely_up,site) : append!(t.lonely_down,site)
end


function SwipeRight(database::ElectronTinder)
    """ Creates list of matches (moves) generated by
    electron tinder's database """
    @assert length(database.lonely_up) == length(database.lonely_down)
    lu = database.lonely_up
    ld = database.lonely_down
    nmatches = length(lu)
    matches = Array{ExchangeMove}(undef,nmatches)
    for i in 1:nmatches
        matches[i] = ExchangeMove((lu[i],ld[i]))
    end
    return matches
end

using Statistics: mean
using LightGraphs: neighborhood
function Swap(chain::DoubleGutzwillerChain,regionsize::Int)

    """ regionsize ≈ r, the radius of the region selected """
    states = get_State(chain) # (state1, state2)
    lg = get_Model(chain)[1]["lattice"].graph # horrible
    ns = prod(get_Model(chain)[1]["dims"]) # also horrible

    swap_per_site = zeros(ns)

    for v in vertices(lg)
        sites = unique(neighborhood(lg,v,regionsize))
        swap_per_site[v] = compute_swapregion(chain,sites)
    end

    return mean(swap_per_site)
end

function compute_swapregion(chain::DoubleGutzwillerChain,
                            sites::Array)
    """ This function is somewhat nontrivial, the following is the
    layout of how it works. First, a list of exchange moves is created
    for each replica that maps the current configuration to the new
    one after swapping. This happens in a two stages: first the spin
    spin configurations inside the regions are compared to see where they
    disagree. Every site where they disagree, an entry is added the
    appropriate spin list in Electron Tinder. Finally, the completed
    ElectronTinder database is fed into SwipeRight to generate a list of moves.

    After this, the list of moves are submitted to preexisting functions
    to compute the determinant ratio produced by that sequence of moves.

    Finally what is returned is the expectation value of the swap operator
    acting on this region which is,
    det(swapped copy 1)*det(swapped copy 2) / (det(og copy 1) * det(og copy 2))
    """
    state = get_State(chain)
    (spinc1, spinc2) = (state[1].spin_config.sc,
                        state[2].spin_config.sc)

    if sum(spinc1[sites])!=sum(spinc2[sites]) # regions must have same quantum #
        return 0
    end
    if spinc1[sites] == spinc2[sites]
        return 1
    end


    Tinder = (ElectronTinder([],[]),ElectronTinder([],[]))
    for site in sites
        if spinc1[site] != spinc2[site]
            addUser!(Tinder[1],site,spinc1[site])
            addUser!(Tinder[2],site,spinc2[site])
        end
    end

    (MoveList1, MoveList2) = SwipeRight.(Tinder)
    @assert size(MoveList1) == size(MoveList2)



    # now because update_chain! modifies the state, make copies of
    # the two states to update
    free_chains = deepcopy(get_Replicas(chain))

    for i in 1:length(MoveList1)
        # process_move!(free_chains[1],MoveList1[i],forceaccept=true)
        # process_move!(free_chains[2],MoveList2[i],forceaccept=true)

        ratio, extras = compute_ratio(free_chains[1],MoveList1[i])
        update_state!(free_chains[1],MoveList1[i],extras=extras)

        ratio, extras = compute_ratio(free_chains[2],MoveList2[i])
        update_state!(free_chains[2],MoveList2[i],extras=extras)
    end

    s0 = state # shorter name to make the det product fit on one line
    ss = (get_State(free_chains[1]), get_State(free_chains[2])) # swapped states


    # determinant products in swapped and original wavefunctions
    sdp = ss[1].det_A_up * ss[1].det_A_down * ss[2].det_A_up * ss[2].det_A_down
    odp = s0[1].det_A_up * s0[1].det_A_down * s0[2].det_A_up * s0[2].det_A_down

    print("sdp: ",sdp," odp: ",odp,"\n")
    return sdp/odp
end