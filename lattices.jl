module Lattices

import LightGraphs

abstract type Lattice end

export Lattice

struct SquareLattice <: Lattice
    dims::Tuple
    graph::LightGraphs.SimpleGraph{Int64}
end

function get_SquareLattice(dims)
    graph = LightGraphs.Grid([d for d in dims],periodic=true)
    return SquareLattice(dims,graph)
end

struct TriangularLattice <: Lattice
    dims::Tuple
    graph::LightGraphs.SimpleGraph{Int64}
end

    
