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
