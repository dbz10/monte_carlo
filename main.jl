#main.jl

include("FreeFermionGutzwiller.jl")
include("lattices.jl")
# include("observables.jl")


# define a model
dims = (6) # dimension of the lattice
lattice = Lattices.get_SquareLattice(dims) # make a square lattice
observable = FreeFermionGutzwiller.observe_Swap()
filling = 1 # setting filling ≢ 1 means there are holes.

model = Dict(
    "dims" => dims,
    "lattice" => lattice,
    "observable" => observable,
    "filling" => filling,
)

# Policy
policy = FreeFermionGutzwiller.get_policy_swap_neighbors()

# MC specifications
mc_specs = Dict(
    "num_mc_warmup_steps" => Int(1e2),
    "num_mc_steps" => Int(1e2),
    "sample_interval" => Int(1e1),
    "num_times" => 5,
    )


# create a markov chain object representing our gutzwiller projected
# free fermion model
gutzwiller_chain = FreeFermionGutzwiller.Gutzwiller_Chain()

# initialize the chain: sets up the initial state
FreeFermionGutzwiller.init_Chain!(gutzwiller_chain,
    model=model,observable=observable,policy=policy,mc_spec=mc_spec)

FreeFermionGutzwiller.runMC!(gutzwiller_chain)
