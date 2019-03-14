include("lattices.jl")
include("FreeFermionGutzwiller.jl")
using LinearAlgebra

lattice = Lattices.get_SquareLattice((4,3))

using GraphPlot
using LightGraphs

nvertices = nv(lattice.graph)
nedges = ne(lattice.graph)
gplot(lattice.graph,nodelabel=1:nvertices,edgelabel=1:nedges)



lv = lattice.lattice_vectors

# for (ind,e) in zip(1:nedges,edges(lattice.graph))
#     u, v = src(e), dst(e)
#     println("edge $ind connects $u - $v")
# end


R_up = [1 4 5 7 8 9]
R_down = [2 3 6 10 11 12]

sco = FreeFermionGutzwiller.get_spin_config(R_up,R_down)



bonds = FreeFermionGutzwiller.get_init_bonds(lattice.graph,sco)
nvertices = nv(bonds)
nedges = ne(bonds)
gplot(bonds,nodelabel=1:nvertices,edgelabel=1:nedges)


# define a model
dims = (4,3) # dimension of the lattice
lattice = Lattices.get_SquareLattice(dims) # make a square lattice
filling = 1 # setting filling ≢ 1 means there are holes.

model = Dict(
    "dims" => dims,
    "lattice" => lattice,
    "filling" => filling,
    "hamiltonian" => FreeFermionGutzwiller.tight_binding_dispersion
    )

gutzwiller_chain = FreeFermionGutzwiller.GutzwillerChain()

FreeFermionGutzwiller.init_Chain!(gutzwiller_chain,model=model)



collect(1:4)

collect(setdiff(Set(collect(1:12)),Set(R_up)))

R_up = zeros(Int64,6)
R_down = zeros(Int64,6)

sites = collect(1:12)

using StatsBase: sample!

sample!(sites,R_up,replace=false)
R_up

leftover_sites=collect(setdiff(Set(sites),Set(R_up)))

sample!(leftover_sites,R_down,replace=false)
