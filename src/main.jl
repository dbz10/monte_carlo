# main.jl

include("FFGmain.jl")
include("lattices.jl")

# define a model
dims = (dl[i],) # dimension of the lattice
lattice = Lattices.get_SquareLattice(dims) # make a square lattice
# observable = FFG.observe_Swap()
filling = 1 # setting filling ≢ 1 means there are holes.





model = Dict(
    "dims" => dims,
    "lattice" => lattice,
    "filling" => filling,
    "hamiltonian" => Lattices.get_Tightbinding_Wavefunctions,
    "fermi_energy" => 0,
    )

# Policy
policy = FFG.SwapNeighborsPolicy()

# MC specifications
mc_spec = Dict("mc_steps" => 1E4*prod(dims),
                "mc_warmup_steps" => 100*prod(dims),
                "sample_interval" => 4*prod(dims))



# create a markov chain object representing our gutzwiller projected
# free fermion model
dg = FFG.DoubleGutzwillerChain()

# initialize the chain: sets up the initial state
FFG.init_Chain!(dg,
    model=model,observable=FFG.MeasureSwap,
    policy=policy,mc_spec=mc_spec)



FFG.runMC!(dg)

# print(FFG.get_Data(gutzwiller_chain))
# print("\n")
# print("Acceptance Ratio: ", FFG.get_Diagnostics(gutzwiller_chain)["acceptance_ratio"])
# print("\n")
print(FFG.get_Data(dg), "\n")
