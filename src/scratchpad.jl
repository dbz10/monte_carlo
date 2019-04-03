include("lattices.jl")
include("FFGmain.jl")

using LinearAlgebra
using GraphPlot
using LightGraphs
using IterTools
using RecursiveArrayTools: convert, VectorOfArray



# define a model
dims = (20,) # dimension of the lattice
lattice = Lattices.get_SquareLattice(dims,pbc=true) # make a square lattice
filling = 1 # setting filling ≢ 1 means there are holes.

VA =VectorOfArray([neighbors(lattice.graph,i)[:] for i in 4:5])
convert(Array,VA)[:]


gplot(lattice.graph,nodelabel=collect(1:nv(lattice.graph)))

model = Dict(
    "dims" => dims,
    "lattice" => lattice,
    "filling" => filling,
    "hamiltonian" => Lattices.get_Tightbinding_Wavefunctions,
    "fermi_energy" => 0,
    )

function observable(chain::FFG.GutzwillerChain)
    state = FFG.get_State(chain)
    Sz1 = state.spin_config.sc[1]
    data = Dict("Sz1" => Sz1)
    return data
end


policy = FFG.SwapNeighborsPolicy()
mc_spec = Dict("mc_steps" => 1E3,
                "mc_warmup_steps" => 1E3,
                "sample_interval" => 1E2)

# gplot(lattice.graph, nodelabel = 1:nv(lattice.graph))
gutz = FFG.GutzwillerChain()

FFG.init_Chain!(gutz,
        model=model,observable=observable,policy=policy,mc_spec=mc_spec
        )


state = FFG.get_State(gutz)

sc = state.spin_config.sc

for i in vertices(lattice.graph)
    print(sc[i])
end

function nth_nearest_neighbors(graph::SimpleGraph{Int64},center_site::Int{64},n::Int{64})
    """ returns indices of vertices that are within n links of the center site.
    algorithm works recursively.
    ex: on a graph representing a 1d chain,
     nth_nearest_neighbors(graph, 5,2) = (3,4,5,6,7)
     """
     neighbors_collection = [center_site]; # I don't think there's any way to preallocate this
     if n == 1
         neighbors_collection = [neighbors_collection ; neighbors(graph,center_site)]
     else
         # l = convert(Array,VectorOfArray([nth_nearest_neighbors(graph,j,n-1) for j in neighbors(center_site)]))[:]
         # this is some ridiculous series of commands to convert an array of lists into one list
         for i in neighbors(center_site)
             neighbors_collection = [neighbors_collection ; nth_nearest_neighbors(graph,i,n-1)]
         end
     end

    return neighbors_collection
end



unique([5; neighbors(lattice.graph,5); 5])
