include("lattices.jl")
include("FreeFermionGutzwiller.jl")

using LinearAlgebra
using GraphPlot



# define a model
dims = (10,) # dimension of the lattice
lattice = Lattices.get_SquareLattice(dims,true) # make a square lattice
filling = 1 # setting filling ≢ 1 means there are holes.


model = Dict(
    "dims" => dims,
    "lattice" => lattice,
    "filling" => filling,
    "hamiltonian" => Lattices.get_Tightbinding_Wavefunctions,
    "fermi_energy" => 0,
    )


# gplot(lattice.graph, nodelabel = 1:nv(lattice.graph))
gutz = FreeFermionGutzwiller.GutzwillerChain()

R_up = [3,2,7,4,1]
R_down = [8,9,10,5,6]

states = FreeFermionGutzwiller.get_States(model)

FreeFermionGutzwiller.get_Determinant(R_up,states)
FreeFermionGutzwiller.get_Determinant(R_down,states)
