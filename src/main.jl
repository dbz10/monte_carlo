#main.jl

include("FreeFermionGutzwiller.jl")
include("lattices.jl")
# include("observables.jl")


# define a model
dims = (400,) # dimension of the lattice
lattice = Lattices.get_SquareLattice(dims) # make a square lattice
# observable = FreeFermionGutzwiller.observe_Swap()
filling = 1 # setting filling ≢ 1 means there are holes.


function Sz1(chain::FreeFermionGutzwiller.GutzwillerChain)
    state = FreeFermionGutzwiller.get_State(chain)
    model = FreeFermionGutzwiller.get_Model(chain)
    Sz1 = state.spin_config.sc[1]
    data = Dict("Sz1" => Float64(Sz1))
    return data
end

using Statistics: mean
function Sz1Sz2(chain::FreeFermionGutzwiller.GutzwillerChain)
    state = FreeFermionGutzwiller.get_State(chain)
    model = FreeFermionGutzwiller.get_Model(chain)
    Sz1Sz2 = mean([state.spin_config.sc[i]*state.spin_config.sc[i+1] for i in 1:length(sc)])
    data = Dict("Sz1Sz2" => Float64(Sz1Sz2))
    return data
end


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
mc_specs = Dict(
    "mc_warmup_steps" => Int(1e5),
    "mc_steps" => Int(1e7),
    "sample_interval" => Int(1e2),
    )



# create a markov chain object representing our gutzwiller projected
# free fermion model
gutzwiller_chain = FreeFermionGutzwiller.GutzwillerChain()

# initialize the chain: sets up the initial state
FreeFermionGutzwiller.init_Chain!(gutzwiller_chain,
    model=model,observable=Sz1Sz2,policy=policy,mc_spec=mc_spec)

FreeFermionGutzwiller.get_Mc_Spec(gutzwiller_chain)


FreeFermionGutzwiller.runMC!(gutzwiller_chain)

print(FreeFermionGutzwiller.get_Data(gutzwiller_chain))
