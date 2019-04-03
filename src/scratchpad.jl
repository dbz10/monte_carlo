include("lattices.jl")
include("FFGmain.jl")

using LinearAlgebra
using GraphPlot
using LightGraphs
using IterTools




# define a model
dims = (20,) # dimension of the lattice
lattice = Lattices.get_SquareLattice(dims,pbc=true) # make a square lattice
filling = 1 # setting filling ≢ 1 means there are holes.



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
gutz = FFG.GutzwillerChain()

FFG.init_Chain!(gutz,
        model=model,observable=observable,policy=policy,mc_spec=mc_spec
        )


state = FFG.get_State(gutz)

sc = state.spin_config.sc

for i in vertices(lattice.graph)
    print(sc[i])
end

dl = [1, 2, 3, 4]
zeros(size(dl))
