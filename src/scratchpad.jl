include("lattices.jl")
include("FFGmain.jl")

using LinearAlgebra
using GraphPlot
using LightGraphs
using Plots



# define a model
dims = (30,) # dimension of the lattice
lattice = Lattices.get_SquareLattice(dims) # make a square lattice
filling = 1 # setting filling ≢ 1 means there are holes.

# gplot(lattice.graph,nodelabel=collect(1:nv(lattice.graph)))

model = Dict(
    "dims" => dims,
    "lattice" => lattice,
    "filling" => filling,
    "hamiltonian" => Lattices.get_Tightbinding_Wavefunctions,
    "fermi_energy" => 0,
    )

# function observable(chain::FFG.GutzwillerChain)
#     state = FFG.get_State(chain)
#     Sz1 = state.spin_config.sc[1]
#     data = Dict("Sz1" => Sz1)
#     return data
# end



policy = FFG.SwapNeighborsPolicy()
mc_spec = Dict("mc_steps" => 1E6*prod(dims),
                "mc_warmup_steps" => 100*prod(dims),
                "sample_interval" => 4*prod(dims))

# gplot(lattice.graph, nodelabel = 1:nv(lattice.graph))
doublegutz = FFG.DoubleGutzwillerChain()
# FFG.init_Chain!(doublegutz,model=model,observable=observable,policy=policy,mc_spec=mc_spec)
# typeof(doublegutz)``

# gutz = FFG.GutzwillerChain()
# FFG.init_Chain!(gutz,model=model,observable=observable,policy=policy,mc_spec=mc_spec)
# FFG.get_MarkovChain(gutz)

# FFG.get_MarkovChain.(doublegutz.replicas)

FFG.init_Chain!(doublegutz,
        model=model,observable=FFG.MeasureSwap,policy=policy,mc_spec=mc_spec
        )




FFG.runMC!(doublegutz)

z = FFG.get_Data(doublegutz)["swap"]
plot(collect(1:2:dims[1]),-log.(z))
# print(-log.(FFG.MeasureSwap(doublegutz)["swap"]))
print(-log.(z))
print(FFG.get_Diagnostics(doublegutz))
