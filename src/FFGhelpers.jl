""" Misc non-math functions used in FreeFermionGutzwiller.jl"""

count_bonds(graph::SimpleGraph{Int64}) = ne(graph)

function get_updated_bonds(chain::GutzwillerChain,move::ExchangeMove)::SimpleGraph{Int64}
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

function update_Rs!(state::GutzwillerState,move::ExchangeMove)
    r_up = state.r_up
    r_down = state.r_down
    bd = state.business_directory
    # use bd to find what electrons are being exchanged
    l1, l2 = move.sites[1], move.sites[2]
    electron_A = bd[l1]
    electron_B = bd[l2]

    if electron_A.spin=="up" # then A is spin up and B is spin down
        r_up[electron_A.label] = l2
        r_down[electron_B.label] = l1
    else # A is down and B is up
        r_down[electron_A.label] = l2
        r_up[electron_B.label] = l1
    end
end

function update_Spin_config!(state::GutzwillerState,move::ExchangeMove)
    sc = state.spin_config.sc
    l1 = move.sites[1]
    l2 = move.sites[2]
    sc[l1], sc[l2] = sc[l2], sc[l1]
end

function update_Bonds!(chain::GutzwillerChain,move::ExchangeMove)
    state = get_State(chain)
    tmp = get_updated_bonds(chain,move)
    state.bonds = tmp
end

function update_Determinants!(chain::GutzwillerChain,move::ExchangeMove,extras)
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

function update_Inverses!(state::GutzwillerState,move::ExchangeMove)
    u1,v1,u2,v2 = get_update_vectors(state,move)
    # recall u1 v1 are for up, u2 v2 are for down
    t = sherman_morrison_inverse_update(state.A_inv_up,u1,v1)
    r = sherman_morrison_inverse_update(state.A_inv_down,u2,v2)

    state.A_inv_up = t
    state.A_inv_down = r
end

function update_Business_directory!(state::GutzwillerState,move::ExchangeMove)
    bd = state.business_directory
    l1 = move.sites[1]
    l2 = move.sites[2]
    bd[l1], bd[l2] = bd[l2], bd[l1]
end

function make_business_directory(r_up::Array,r_down::Array)::Array
    """ Makes the directory for the lattice. The i'th element
    is the index of the electron living on site i. indices <= N/2
    refer to up spins and indices > N/2 refer to down spins """
    ns = length([r_up ; r_down])
    bd = Array{Electron}(undef, ns)

    for i in 1:length(r_up) # assume len(r_up) = len(r_down)
        bd[r_up[i]] = Electron(i,"up")
        bd[r_down[i]] = Electron(i,"down")
    end

    return bd
end

function get_ExchangeMove(bonds::SimpleGraph{Int64})::ExchangeMove
    """ draw a random bond and return the two sites attached to that bond """
    index = rand(1:ne(bonds))
    edge = nth(edges(bonds),index) #nth is from IterTools
    move_sites = (src(edge),dst(edge))
    return ExchangeMove(move_sites)
end

function get_update_vectors(state::GutzwillerState,move::ExchangeMove)
    """ Gets vectors which perform rank 1 row update to the slater determinant
    specifically, u is zero except at column i and v_i = ϕ_i(r') - ϕ_i(r)
    u1 and v1 are for up spin, and u2 and v2 are for down spin"""
    # get electron indices from business directory
    bd = state.business_directory

    site_A = move.sites[1]
    site_B = move.sites[2]

    electron_A = bd[site_A]
    electron_B = bd[site_B]

    u1 = zeros(length(state.r_up))
    v1 = zeros(length(state.r_up))
    u2 = zeros(length(state.r_down))
    v2 = zeros(length(state.r_down))

    if electron_A.spin=="up"
        # then A is up (1) and B is down (2)
        # up electron is moving from site A to site B
        u1[electron_A.label] = 1
        v1 = state.wavefunctions[site_B,:] - state.wavefunctions[site_A,:]

        # down electron is moving from site B to site A
        u2[electron_B.label] = 1
        v2 = state.wavefunctions[site_A,:] - state.wavefunctions[site_B,:]

    else
        # A is down (2) and B is up (1)
        u2[electron_A.label] = 1
        v2 = state.wavefunctions[site_B,:] - state.wavefunctions[site_A,:]

        u1[electron_B.label] = 1
        v1 = state.wavefunctions[site_A,:] - state.wavefunctions[site_B,:]
    end

    return u1,v1,u2,v2
end
