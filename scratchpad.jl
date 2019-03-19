include("lattices.jl")
include("FreeFermionGutzwiller.jl")
using LinearAlgebra
using GraphPlot



# define a model
dims = (5,) # dimension of the lattice
lattice = Lattices.get_SquareLattice(dims,true) # make a square lattice
filling = 1 # setting filling ≢ 1 means there are holes.


model = Dict(
    "dims" => dims,
    "lattice" => lattice,
    "filling" => filling,
    )


gplot(lattice.graph, nodelabel = 1:nv(lattice.graph))

F = Lattices.get_Tightbinding_Wavefunctions(lattice)
F.values
F.vectors
