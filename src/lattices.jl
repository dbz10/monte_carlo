# Module for constructing lattices
# Things to add:
#   - support for unit cell
#   - support for generic lattices (so far only have square)

module Lattices

import LightGraphs
using LinearAlgebra
using LightGraphs: adjacency_matrix

abstract type Lattice end

export Lattice

struct SquareLattice <: Lattice
    dims::Tuple
    graph::LightGraphs.SimpleGraph{Int64}
    lattice_vectors::Array{Int64}
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


function get_Tightbinding_Wavefunctions(lattice)::Eigen
    ham = -Matrix(adjacency_matrix(lattice.graph))
    F = eigen(ham) # F.values are evals, F.vectors[:,k] is k'th evec.
    return F
end

end

function nth_nearest_neighbors(graph::SimpleGraph{Int64},center_site::Int,n::Int)
    """ returns indices of vertices that are within n links of the center site.
    algorithm works recursively.
    ex: on a graph representing a 1d chain,
     nth_nearest_neighbors(graph, 5,2) = (3,4,5,6,7)
     """
     neighbors_collection = [center_site]; # I don't think there's any way to preallocate this
     if n == 1
         neighbors_collection = [neighbors_collection ; neighbors(graph,center_site)]
     else
         for i in neighbors(graph,center_site)
             neighbors_collection = [neighbors_collection ; nth_nearest_neighbors(graph,i,n-1)]
         end
     end

    return unique(neighbors_collection)
end
