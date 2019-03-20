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
