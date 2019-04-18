# Module for constructing lattices
# Things to add:
#   - support for unit cell
#   - support for generic lattices (so far only have square)

module Lattices


using LinearAlgebra
using LightGraphs
using SimpleWeightedGraphs
import LightGraphs.cartesian_product

abstract type Lattice end

export Lattice

struct SquareLattice <: Lattice
    dims::Tuple
    graph::LightGraphs.SimpleGraph{Int64}
    lattice_vectors::Array{Int64}
end

struct MixedBCSquareLattice <: Lattice
    dims::Tuple
    graph::SimpleWeightedGraphs.SimpleWeightedGraph
end

function get_MixedBoundarySquareLattice(dims::Tuple)::MixedBCSquareLattice
    """ This is an extension of the function Grid from LightGraphs to
    make a weighted graph with antiperiodic boundary conditions in one direction
    and periodic boundary conditions in the other, for the purposes
    of obtaining a nondegenerate fermi surface.
    https://journals.jps.jp/doi/pdf/10.1143/JPSJ.57.2482 """

   # make anti periodic cycle in x direction
   g = SimpleWeightedGraph(CycleGraph(dims[1]))
   add_edge!(g,1,dims[1],-1)

   for d in dims[2:end]
       g = cartesian_product(SimpleWeightedGraph(CycleGraph(d)), g)
   end
   return MixedBCSquareLattice(dims,g)
end

function get_SquareLattice(dims;pbc=true)::SquareLattice
    graph = LightGraphs.Grid([d for d in dims],periodic=pbc)
    # lattice vectors for a d dimensional square lattice are just
    # v_i = \hat e_i
    lattice_vectors = get_SquareLatticeVectors(dims)
    return SquareLattice(dims,graph,lattice_vectors)
end

function get_SquareLatticeVectors(dims::Tuple)::Array{Int64}
    """ returns array of lattice vectors
    v[i,:,:] is the i'th lattice vector"""
    ndims = length(dims)
    m0 = Matrix{Int}(I,ndims,ndims)
    return m0
end


struct TriangularLattice <: Lattice
    dims::Tuple
    graph::LightGraphs.SimpleGraph{Int64}
end


function get_Tightbinding_Wavefunctions(lattice::SquareLattice)::Eigen
    ham = -Matrix(adjacency_matrix(lattice.graph))
    F = eigen(ham) # F.values are evals, F.vectors[:,k] is k'th evec.
    return F
end

function get_Tightbinding_Wavefunctions(lattice::MixedBCSquareLattice)::Eigen
    ham = -Matrix(weights(lattice.graph))
    F = eigen(ham) # F.values are evals, F.vectors[:,k] is k'th evec.
    return F
end

function get_2DCSL_Wavefunctions(lattice::SquareLattice,Δ)::Eigen
    lg = lattice.graph
    dims = lattice.dims
    nsites = prod(dims)
    ham = zeros(ComplexF64,(nsites,nsites))
    Lx = dims[1]
    Ly = dims[2]

    for y = 1:dims[2]
        # loop over y layers. odd layers have hopping -t in the x dir
        # even layers have hopping +t in the x dir
        for x = 1:dims[1]
            id0 = sub2ind(dims,x,y)

            if (y % 2) == 1
                idp = sub2ind(dims,x,wrap(y+1,Ly))
                ham[idp,id0] = -1

                idp = sub2ind(dims,wrap(x+1,Lx),y)
                ham[idp,id0] = +1

                idp = sub2ind(dims,wrap(x+1,Lx),wrap(y+1,Ly))
                ham[idp,id0] = +1im*Δ

                idp = sub2ind(dims,wrap(x+1,Lx),wrap(y-1,Ly))
                ham[idp,id0] = -1im*Δ
            else
                idp = sub2ind(dims,x,wrap(y+1,Ly))
                ham[idp,id0] = -1

                idp = sub2ind(dims,wrap(x+1,Lx),y)
                ham[idp,id0] = -1

                idp = sub2ind(dims,wrap(x+1,Lx),wrap(y+1,Ly))
                ham[idp,id0] = -1im*Δ

                idp = sub2ind(dims,wrap(x+1,Lx),wrap(y-1,Ly))
                ham[idp,id0] = +1im*Δ
            end
        end
    end

    ham += ham'
    # print(ham)
    F = eigen(ham)
    return F
end

function wrap(x,a)
    """ wraps x to between 1 and a """
    return fldmod(x-1,a)[2]+1
end

"""transition to julia 1.0"""
function ind2sub(dims, ind)
    return Tuple(CartesianIndices(dims)[ind])
end

"""transition to julia 1.0"""
function sub2ind(dims, inds...)
    return (LinearIndices(dims))[inds...]
end



function cartesian_product(g::G, h::G) where G <: AbstractSimpleWeightedGraph
    """ Function overload for cartesian product for weighted graphs. For
    some reason, not correctly loading override functions from package..."""
    z = G(nv(g) * nv(h))
    id(i, j) = (i - 1) * nv(h) + j
    for e in edges(g)
        i1, i2 = Tuple(e)
        for j = 1:nv(h)
            add_edge!(z, id(i1, j), id(i2, j), weight(e))
        end
    end

    for e in edges(h)
        j1, j2 = Tuple(e)
        for i in vertices(g)
            add_edge!(z, id(i, j1), id(i, j2), weight(e))
        end
    end
    return z
end

end
