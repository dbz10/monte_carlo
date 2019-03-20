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
observable = FreeFermionGutzwiller.GutzwillerObservable()
policy = FreeFermionGutzwiller.SwapNeighborsPolicy()
mc_spec = Dict("mc_steps" => 1)

# gplot(lattice.graph, nodelabel = 1:nv(lattice.graph))
gutz = FreeFermionGutzwiller.GutzwillerChain()

FreeFermionGutzwiller.init_Chain!(gutz,
        model=model,observable=observable,policy=policy,mc_spec=mc_spec
        )

R_up = [3,2,7,4,1]
R_down = [8,9,10,5,6]

states = FreeFermionGutzwiller.get_States(model)

FreeFermionGutzwiller.get_Determinant(R_up,states)
FreeFermionGutzwiller.get_Determinant(R_down,states)

abstract type dnum end

mutable struct num <: dnum
    val::Float64
end

function add(v::dnum,z::dnum)
    return v.val + z.val
end

a = num(1)

b=num(3)

add(a,b)
