include("lattices.jl")
include("FreeFermionGutzwiller.jl")
using LinearAlgebra
using GraphPlot
using LightGraphs: adjacency_matrix



# define a model
dims = (5,) # dimension of the lattice
lattice = Lattices.get_SquareLattice(dims,true) # make a square lattice
filling = 1 # setting filling ≢ 1 means there are holes.


model = Dict(
    "dims" => dims,
    "lattice" => lattice,
    "filling" => filling,
    "hamiltonian" => FreeFermionGutzwiller.tight_binding_dispersion
    )


gplot(lattice.graph, nodelabel = 1:nv(lattice.graph))

v1 = FreeFermionGutzwiller.get_Wavefunctions(lattice)
v1 = round.(FreeFermionGutzwiller.get_Wavefunctions(lattice)*100)/100

k = 2*pi*collect(1:5)/5

r = collect(1:5)

v2 = round.(100*exp.([1im*ri*kj for kj in k, ri in r])./sqrt(5))/100
v2 = exp.([1im*ri*kj for kj in k, ri in r])./sqrt(5)


ham = Matrix(adjacency_matrix(lattice.graph))

ham
