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
    n_up,n_down,
    model,conditioning_tol=10^6
    )
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
        filled_states = get_Wavefunctions(model)

        cc_up = cond(filled_states[R_up,:])
        cc_down = cond(filled_states[R_down,:])

        if (cc_up < conditioning_tol && cc_down < conditioning_tol)
            return R_up, R_down, filled_states
        end
    end
end
