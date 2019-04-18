include("lattices.jl")
include("FFGmain.jl")

using LinearAlgebra
using GraphPlot
using LightGraphs
using Plots
using JLD2, FileIO
using Statistics



# define a model
dims = (12,12) # dimension of the lattice
lattice = Lattices.get_SquareLattice(dims,pbc=true) # make a square lattice
filling = 1 # setting filling ≢ 1 means there are holes.

# gplot(SimpleGraph(adjacency_matrix(lattice.graph)),nodelabel=collect(1:nv(lattice.graph)))

F = Lattices.get_2DCSL_Wavefunctions(lattice,0.5)
histogram(F.values,bins=50)
plot(F.values)

Δ = 0.5
model = Dict(
    "dims" => dims,
    "lattice" => lattice,
    "filling" => filling,
    "hamiltonian" => Lattices.get_2DCSL_Wavefunctions(lattice,Δ),
    "fermi_energy" => 0,
    )

# function observable(chain::FFG.GutzwillerChain)
#     state = FFG.get_State(chain)
#     Sz1 = state.spin_config.sc[1]
#     data = Dict("Sz1" => Sz1)
#     return data
# end

n_runs = 10

policy = FFG.SwapNeighborsPolicy()
mc_spec = Dict("mc_steps" => 1E3*prod(dims),
                "mc_warmup_steps" => 100*prod(dims),
                "sample_interval" => 5*prod(dims))

x = collect(1:2:dims[1])
data = zeros(n_runs,length(x))
for j = 1:n_runs
    doublegutz = FFG.DoubleGutzwillerChain()
    FFG.init_Chain!(doublegutz,
            model=model,observable=FFG.MeasureSwap,policy=policy,mc_spec=mc_spec
            )
    FFG.runMC!(doublegutz)
    z = FFG.get_Data(doublegutz)["swap"]
    data[j,:] = z
end


d = -log.(data)

y = mean(d,dims=1)'
err = std(d,dims=1)'
plotly()
plot(x,y,yerror=err)
# sites = unique(neighborhood(lattice.graph,39,1))
# print(FFG.compute_swapregion(doublegutz,sites),"\n")
# #




FFG.runMC!(doublegutz)
# FFG.aggregate_measurements!(doublegutz)

z = FFG.get_Data(doublegutz)["swap"]
l = collect(1:2:dims[1])

plot(l,-log.(z))
# print(-log.(FFG.MeasureSwap(doublegutz)["swap"]))
print(-log.(z))
print(FFG.get_Diagnostics(doublegutz))

z = [1im 0.5im; 2 2im]

z+=z'
