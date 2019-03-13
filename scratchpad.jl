include("lattices.jl")
using LinearAlgebra

lattice = Lattices.get_SquareLattice((6,2))

using GraphPlot
using LightGraphs

gplot(lattice.graph)

lv = lattice.lattice_vectors


k = (pi/6, pi/20)
k_vector = collect(k)
lv[2,:,:]

terms = [exp(1im*dot(k_vector,lv[i,:,:])) for i in 1:length(k)]
sum(terms +conj(terms))


-2*cos(pi/6)-2cos(pi/20)
