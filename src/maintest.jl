# main.jl runs on small systems sizes, to test that code is fully run-able


function maintest()
    dl = [4,6,8]
    cl = zeros(size(dl))


    for i in 1:length(dl)
    # define a model
    dims = (dl[i],dl[i]) # dimension of the lattice
    lattice = Lattices.get_SquareLattice(dims) # make a square lattice
    filling = 1 # setting filling ≢ 1 means there are holes.





    model = Dict(
        "dims" => dims,
        "lattice" => lattice,
        "filling" => filling,
        "hamiltonian" => Lattices.get_Tightbinding_Wavefunctions,
        "fermi_energy" => 0,
        )

    policy = FFG.SwapNeighborsPolicy()


    mc_spec = Dict(
        "mc_warmup_steps" => Int(1e2*prod(dims)),
        "mc_steps" => Int(1e3*prod(dims)),
        "sample_interval" => Int(prod(dims)*4),
        )




    gutzwiller_chain = FFG.GutzwillerChain()


    FFG.init_Chain!(gutzwiller_chain,
        model=model,observable=FFG.NeighborSzSz,
        policy=policy,mc_spec=mc_spec)

    FFG.get_Mc_Spec(gutzwiller_chain)


    FFG.runMC!(gutzwiller_chain)


    cl[i] = FFG.get_Data(gutzwiller_chain)["SzSz"]*3;
    end


    return cl
end
