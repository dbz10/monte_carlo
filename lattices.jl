# Module for constructing lattices
# Things to add:
#   - support for unit cell
#   - support for generic lattices (so far only have square)

module Lattices

import LightGraphs
using LinearAlgebra


abstract type Lattice end

export Lattice

struct SquareLattice <: Lattice
    dims::Tuple
    graph::LightGraphs.SimpleGraph{Int64}
    lattice_vectors::Array{Int64}
end

function get_SquareLattice(dims)::SquareLattice
    graph = LightGraphs.Grid([d for d in dims],periodic=false)
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


end
