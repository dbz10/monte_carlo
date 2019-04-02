#main.jl

include("FreeFermionGutzwiller.jl")
include("lattices.jl")
# include("observables.jl")


# define a model
dims = (8,8) # dimension of the lattice
lattice = Lattices.get_SquareLattice(dims) # make a square lattice
# observable = FreeFermionGutzwiller.observe_Swap()
filling = 1 # setting filling ≢ 1 means there are holes.





model = Dict(
    "dims" => dims,
    "lattice" => lattice,
    "filling" => filling,
    "hamiltonian" => Lattices.get_Tightbinding_Wavefunctions,
    "fermi_energy" => 0,
    )

# Policy
policy = FreeFermionGutzwiller.SwapNeighborsPolicy()

# MC specifications
mc_spec = Dict(
    "mc_warmup_steps" => Int(1e2*prod(dims)),
    "mc_steps" => Int(1e4*prod(dims)),
    "sample_interval" => Int(prod(dims)*4),
    )



# create a markov chain object representing our gutzwiller projected
# free fermion model
gutzwiller_chain = FreeFermionGutzwiller.GutzwillerChain()

# initialize the chain: sets up the initial state
FreeFermionGutzwiller.init_Chain!(gutzwiller_chain,
    model=model,observable=FreeFermionGutzwiller.NeighborSzSz,
    policy=policy,mc_spec=mc_spec)

FreeFermionGutzwiller.get_Mc_Spec(gutzwiller_chain)


FreeFermionGutzwiller.runMC!(gutzwiller_chain)

# print(FreeFermionGutzwiller.get_Data(gutzwiller_chain))
# print("\n")
# print("Acceptance Ratio: ", FreeFermionGutzwiller.get_Diagnostics(gutzwiller_chain)["acceptance_ratio"])
# print("\n")
print(FreeFermionGutzwiller.get_Data(gutzwiller_chain)["SzSz"]*3)
