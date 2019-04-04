# main.jl

include("FFGmain.jl")
include("lattices.jl")

dl = [6,10,14]
cl = zeros(size(dl))

for i in 1:length(dl)
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
mc_spec = Dict(
    "mc_warmup_steps" => Int(1e2*prod(dims)),
    "mc_steps" => Int(1e4*prod(dims)),
    "sample_interval" => Int(prod(dims)*4),
    )



# create a markov chain object representing our gutzwiller projected
# free fermion model
gutzwiller_chain = FFG.GutzwillerChain()

# initialize the chain: sets up the initial state
FFG.init_Chain!(gutzwiller_chain,
    model=model,observable=FFG.NeighborSzSz,
    policy=policy,mc_spec=mc_spec)

FFG.get_Mc_Spec(gutzwiller_chain)


FFG.runMC!(gutzwiller_chain)

# print(FFG.get_Data(gutzwiller_chain))
# print("\n")
# print("Acceptance Ratio: ", FFG.get_Diagnostics(gutzwiller_chain)["acceptance_ratio"])
# print("\n")
print(FFG.get_Data(gutzwiller_chain)["SzSz"]*3, "\n")
cl[i] = FFG.get_Data(gutzwiller_chain)["SzSz"]*3;
end

print(dl,"\n")
print(cl)
