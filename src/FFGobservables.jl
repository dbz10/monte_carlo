""" Observables that can be measured on a gutzwiller state """

function Szi(chain::GutzwillerChain,i)
    state = get_State(chain)
    model = get_Model(chain)
    Sz = state.spin_config.sc[i]
    data = Dict("Szi" => Float64(Szi))
    return data
end

using Statistics: mean
function NeighborSzSz(chain::GutzwillerChain)
    state = get_State(chain)
    model = get_Model(chain)
    lg = model["lattice"].graph # latticegraph
    spinc = state.spin_config.sc # spin config array

    SzSz = 0

    for i in vertices(lg)
        SzSz += mean([spinc[i]*spinc[j]/4 for j in neighbors(lg,i)])
    end

    data = Dict("SzSz" => Float64(SzSz)/nv(lg))
    return data
end

function MeasureSwap(chain::DoubleGutzwillerChain)
    model = get_Model(chain)
    ns = prod(model["dims"])

    len = Int64(ns/2+2)

    sw = zeros(len)
    sw[1] = 1 # manually put in swap(nothing) = 1
    for i = 2:len
        sw[i] = Swap(chain,i-2) # regionsize 0 means just 1 site
    end

    data = Dict("swap" => sw)
    return data
end
