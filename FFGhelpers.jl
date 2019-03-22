""" Misc non-math functions used in FreeFermionGutzwiller.jl"""

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

function get_electron_index(site::Int64,state::GutzwillerState)::Int64
    ns = length(state.business_directory)
    return ((site-1) % (ns/2)) + 1
end

function get_conditioned_state(
    chain, n_up,n_down, conditioning_tol=10^6)
    sites = collect(1:n_up+n_down)
    cc_up = 2*conditioning_tol
    cc_down = 2*conditioning_tol
    while (cc_up > conditioning_tol  || cc_down > conditioning_tol)
        # vectors to store the location of the up and down electrons
        R_up = zeros(Int64, n_up)
        R_down = zeros(Int64, n_down)
        sample!(sites,R_up,replace=false)
        leftover_sites = collect(setdiff(Set(sites),Set(R_up)))
        sample!(leftover_sites,R_down,replace=false)

                # get occupied wavefunctions from lattice and fermi energy
        filled_states = get_Wavefunctions(chain)

        cc_up = cond(filled_states[R_up,:])
        cc_down = cond(filled_states[R_down,:])

        if (cc_up < conditioning_tol && cc_down < conditioning_tol)
            return R_up, R_down, filled_states
        end
    end
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

function update_Determinants!(chain::GutzwillerChain,move::SwapNeighborMove,extras)
    # extras contains det_ratio_up and det_ratio_down

    state = get_State(chain)

    if isnothing(extras)
        trash, extras = compute_ratio(chain,move)
    end
    det_ratio_up = extras[1]
    det_ratio_down = extras[2]
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

function make_business_directory(r_up::Array,r_down::Array)::Array
    """ Makes the directory for the lattice. The i'th element
    is the index of the electron living on site i. indices <= N/2
    refer to up spins and indices > N/2 refer to down spins """
    r_full = [r_up ; r_down]
    sortkey = sortperm(r_full)
    return sortkey
end

function get_SwapNeighborMove(bonds::SimpleGraph{Int64})::SwapNeighborMove
    """ draw a random bond and return the two sites attached to that bond """
    index = rand(1:ne(bonds))
    edge = nth(edges(bonds),index) #nth is from IterTools
    move_sites = (src(edge),dst(edge))
    return SwapNeighborMove(move_sites)
end
