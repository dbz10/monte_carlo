include("FreeFeermionGutzwiller.jl")
include("lattices.jl")

lattice = Lattices.get_SquareLattice((4,4))

using GraphPlot
using LightGraphs

gplot(lattice.graph)
