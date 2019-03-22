include("lattices.jl")
include("FreeFermionGutzwiller.jl")

using LinearAlgebra
using GraphPlot
using LightGraphs
using IterTools




# define a model
dims = (16,16) # dimension of the lattice
lattice = Lattices.get_SquareLattice(dims,true) # make a square lattice
filling = 1 # setting filling ≢ 1 means there are holes.


gplot(lattice.graph,nodelabel=collect(1:nv(lattice.graph)))

model = Dict(
    "dims" => dims,
    "lattice" => lattice,
    "filling" => filling,
    "hamiltonian" => Lattices.get_Tightbinding_Wavefunctions,
    "fermi_energy" => 0,
    )

function observable(chain::FreeFermionGutzwiller.GutzwillerChain)
    state = FreeFermionGutzwiller.get_State(chain)
    Sz1Sz2 = state.spin_config.sc[1]*state.spin_config.sc[2]
    data = Dict("Sz1Sz2" => Sz1Sz2)
    return data
end


policy = FreeFermionGutzwiller.SwapNeighborsPolicy()
mc_spec = Dict("mc_steps" => 1,
                "mc_warmup_steps" => 1E3)

# gplot(lattice.graph, nodelabel = 1:nv(lattice.graph))
gutz = FreeFermionGutzwiller.GutzwillerChain()

FreeFermionGutzwiller.init_Chain!(gutz,
        model=model,observable=observable,policy=policy,mc_spec=mc_spec
        )

@time FreeFermionGutzwiller.do_warmup!(gutz)

FreeFermionGutzwiller.get_Diagnostics(gutz)
