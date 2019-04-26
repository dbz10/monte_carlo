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

# F = Lattices.get_2DCSL_Wavefunctions(lattice,0.25)
# histogram(F.values,bins=50)
# plot(F.values)

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

n_runs = 300

policy = FFG.SwapNeighborsPolicy()
mc_spec = Dict("mc_steps" => 3E3*prod(dims),
                "mc_warmup_steps" => 300*prod(dims),
                "sample_interval" => 10*prod(dims))

x = collect(1:2:dims[1])
data = zeros(ComplexF64,n_runs,length(x))
for j = 1:n_runs
    println("On run ",j)
    doublegutz = FFG.DoubleGutzwillerChain()
    FFG.init_Chain!(doublegutz,
            model=model,observable=FFG.MeasureSwap,policy=policy,mc_spec=mc_spec
            )
    FFG.runMC!(doublegutz)
    z = FFG.get_Data(doublegutz)["swap"]
    data[j,:] = z
end

real(1+2im)
print(mean(data,dims=1))



d = -log.(abs.(real.(data)))
y = -log.(abs.(mean(data,dims=1)))'

err = std(d,dims=1)'
plotly()
plot(x,y)

z = 2
# sites = unique(neighborhood(lattice.graph,39,1))
# print(FFG.compute_swapregion(doublegutz,sites),"\n")
# #


#
#
# FFG.runMC!(doublegutz)
# # FFG.aggregate_measurements!(doublegutz)
#
# z = FFG.get_Data(doublegutz)["swap"]
# l = collect(1:2:dims[1])
#
# plot(l,-log.(z))
# # print(-log.(FFG.MeasureSwap(doublegutz)["swap"]))
# print(-log.(z))
# print(FFG.get_Diagnostics(doublegutz))
#
# z = [1im 0.5im; 2 2im]
#
# z = zeros((2,2))
#
# typeof(z)
#
# z[1,1] = 1im
#
# z = 1.04 + 0.23im
#
# t = typeof(z)
#
# z = zeros(ComplexF64,(3,3))
#
# typeof(z) == Complex{Float64}
#
# z = zeros(t,(3,3))
#
