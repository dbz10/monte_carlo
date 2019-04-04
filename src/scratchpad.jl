include("lattices.jl")
include("FFGmain.jl")

using LinearAlgebra
using GraphPlot
using LightGraphs



# define a model
dims = (20,) # dimension of the lattice
lattice = Lattices.get_SquareLattice(dims,pbc=true) # make a square lattice
filling = 1 # setting filling ≢ 1 means there are holes.


sort(unique(neighborhood(lattice.graph,3,2)))

gplot(lattice.graph,nodelabel=collect(1:nv(lattice.graph)))

model = Dict(
    "dims" => dims,
    "lattice" => lattice,
    "filling" => filling,
    "hamiltonian" => Lattices.get_Tightbinding_Wavefunctions,
    "fermi_energy" => 0,
    )

function observable(chain::FFG.GutzwillerChain)
    state = FFG.get_State(chain)
    Sz1 = state.spin_config.sc[1]
    data = Dict("Sz1" => Sz1)
    return data
end


policy = FFG.SwapNeighborsPolicy()
mc_spec = Dict("mc_steps" => 1E3,
                "mc_warmup_steps" => 1E3,
                "sample_interval" => 1E2)

# gplot(lattice.graph, nodelabel = 1:nv(lattice.graph))
doublegutz = FFG.DoubleGutzwillerChain()
# FFG.init_Chain!(doublegutz,model=model,observable=observable,policy=policy,mc_spec=mc_spec)
typeof(doublegutz)

FFG.init_Chain!.(doublegutz.replicas,
        model=model,observable=observable,policy=policy,mc_spec=mc_spec
        )

(s1,s2) = FFG.get_State(doublegutz)

print(FFG.MeasureSwap(doublegutz))


et = FFG.ElectronTinder([1],[2])
et = FFG.ElectronTinder([],[])

FFG.addUser!(et,1,1)
FFG.addUser!(et,2,-1)
FFG.addUser!(et,3,1)
FFG.addUser!(et,4,-1)

et
m = FFG.SwipeRight(et)
