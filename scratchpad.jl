include("lattices.jl")
include("FreeFermionGutzwiller.jl")

using LinearAlgebra
using GraphPlot
using LightGraphs
using IterTools



# define a model
dims = (2,3) # dimension of the lattice
lattice = Lattices.get_SquareLattice(dims,false) # make a square lattice
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
    thingy = state.spin_config.sc[1]
    data = Dict("Sz1" => thingy)
    return data
end


policy = FreeFermionGutzwiller.SwapNeighborsPolicy()
mc_spec = Dict("mc_steps" => 1)

# gplot(lattice.graph, nodelabel = 1:nv(lattice.graph))
gutz = FreeFermionGutzwiller.GutzwillerChain()

FreeFermionGutzwiller.init_Chain!(gutz,
        model=model,observable=observable,policy=policy,mc_spec=mc_spec
        )

state = FreeFermionGutzwiller.get_State(gutz)
gplot(state.bonds,nodelabel=collect(1:nv(lattice.graph)))


move=FreeFermionGutzwiller.get_Move(gutz)
updated_bonds = FreeFermionGutzwiller.get_updated_bonds(gutz,move)
gplot(updated_bonds,nodelabel=collect(1:nv(lattice.graph)))

print(state.spin_config.sc)
