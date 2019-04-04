using Test
using LightGraphs
using LinearAlgebra
include("FFGmain.jl")
include("lattices.jl")


ffgmathtest = @testset "test-FFGMATH/LATTICES" begin
    ############### testing lattice routines #################################
    dims = (5,5);
    lg = Lattices.get_SquareLattice(dims,pbc=true).graph
    @test sort(unique(neighborhood(lg, 3,2))) ==
    sort([3,4,5,2,1,9,8,7,13,23,18,24,22])

    dims = (4,4)
    lattice = Lattices.get_MixedBoundarySquareLattice(dims)
    F = Lattices.get_Tightbinding_Wavefunctions(lattice)
    @test F.values ≈ [-3.414213562373095,-3.414213562373095,
    -1.4142135623730951,-1.4142135623730951,-1.4142135623730951,-1.4142135623730951,
    -0.5857864376269051,-0.5857864376269051,0.5857864376269051,0.5857864376269051,
    1.4142135623730951,1.4142135623730951,1.4142135623730951,1.4142135623730951,
    3.414213562373095,3.414213562373095]

    dims = (2,3);
    lattice = Lattices.get_SquareLattice(dims,pbc=true);
    @test Matrix(adjacency_matrix(lattice.graph)) == [0 1 1 0 1 0;
                                                     1 0 0 1 0 1;
                                                     1 0 0 1 1 0;
                                                     0 1 1 0 0 1;
                                                     1 0 1 0 0 1;
                                                     0 1 0 1 1 0]
    F = Lattices.get_Tightbinding_Wavefunctions(lattice);
    evals = F.values;
    @test evals ≈ [-3, -1, 0, 0, 2, 2]
    states = F.vectors[:,1:3]
    r = [1,2,3]
    ############### testing math routines #####################################
    @test FFG.get_Determinant(r,states) ≈ 0.3333333333
    tmat = [0 -1.2247448713915894 -1.2247448713915894;
            -1.2247448713915894 1.2247448713915894 0;
            1 0 -1]
    @test FFG.get_Inverse_Matrix(r,states) ≈ tmat

    ############### set up small matrix to test det ratio factor & inverse #####
    m1 = [1 2 4;
          4 8 12;
          8 12 16];
    m1inv = inv(m1);
    m2 = [1 4 4;
           4 12 12;
           8 15 16];

    u = [2 ;4 ;3 ];
    v = [0 ; 1 ; 0 ];

    @test FFG.det_ratio_factor(m1inv,u,v) ≈ det(m2)/det(m1)

    @test FFG.sherman_morrison_inverse_update(m1inv,u,v) ≈ inv(m2)
end

@testset "test-FFGMAIN/HELPERS" begin
    ############ make a test chain $$$$$$$$$$$$$$$$$$$$$$$$

    model = Dict(
        "dims" => (2,3),
        "lattice" => Lattices.get_SquareLattice((2,3),pbc=true),
        "filling" => 1,
        "hamiltonian" => Lattices.get_Tightbinding_Wavefunctions,
        "fermi_energy" => 0,
        )

    function observable(chain::FFG.GutzwillerChain)
        state = FFG.get_State(chain)
        Sz1 = state.spin_config.sc[1]
        return Dict("Sz1" => Sz1)
    end

    policy = FFG.SwapNeighborsPolicy()
    mc_spec = Dict("mc_steps" => 1,
                    "mc_warmup_steps" => 1)

    testchain = FFG.GutzwillerChain();
    testchain.basechain = FFG.MarkovChain()


    r_up = [3,1,5]
    r_down = [2,6,4]


    # STATE CREATION TESTS
    bd = [FFG.Electron(2,"up"), FFG.Electron(1,"down"), FFG.Electron(1,"up"),
        FFG.Electron(3,"down"), FFG.Electron(3,"up"), FFG.Electron(2,"down")]
    @test FFG.make_business_directory(r_up,r_down) == bd

    spinc = FFG.SpinConfiguration([1,-1,1,-1,1,-1])
    @test FFG.make_spin_config(r_up,r_down).sc == spinc.sc

    # testing construction of initial bond graph
    using LightGraphs: SimpleGraph, add_edge!
    ib = SimpleGraph(6);
    add_edge!(ib,1,2);
    add_edge!(ib,3,4);
    add_edge!(ib,5,6);

    @test FFG.get_init_bonds(model["lattice"].graph,spinc) ==
    ib

    ## STATE UPDATE TESTS
    teststate = FFG.GutzwillerState(
    r_up,r_down, bd, spinc,ib,
    [], 0, 0, [],[]
    )

    testchain.basechain.state = teststate
    testchain.basechain.model = model

    ## test a simple move
    move = FFG.ExchangeMove((1,2))

    @test FFG.get_proposal_factor_ratio(testchain,move) == 3/7
    FFG.update_Rs!(teststate,move)
    @test teststate.r_up == [3,2,5]
    @test teststate.r_down == [1,6,4]



    FFG.update_Bonds!(testchain,move)

    ub = SimpleGraph(6);
    add_edge!(ub,1,2);
    add_edge!(ub,1,3);
    add_edge!(ub,1,5);
    add_edge!(ub,2,4);
    add_edge!(ub,2,6);
    add_edge!(ub,3,4);
    add_edge!(ub,5,6);
    @test teststate.bonds == ub

    FFG.update_Spin_config!(teststate,move)
    @test teststate.spin_config.sc == [-1,1,1,-1,1,-1]

    bd = [FFG.Electron(1,"down"),FFG.Electron(2,"up"),FFG.Electron(1,"up"),
        FFG.Electron(3,"down"),FFG.Electron(3,"up"),FFG.Electron(2,"down")]
    FFG.update_Business_directory!(teststate,move)
    @test teststate.business_directory == bd

end


@testset "maintest" begin
    include("maintest.jl")
    @test size(maintest()) == (3,)
end

@testset "swaproutines" begin
    """ Test suite for swap routines. """
    dims = (6,) # dimension of the lattice
    lattice = Lattices.get_SquareLattice(dims,pbc=true) # make a square lattice
    filling = 1 # setting filling ≢ 1 means there are holes.
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


    dg = FFG.DoubleGutzwillerChain()
    c1 = dg.replicas[1]
    c2 = dg.replicas[2]
    c1.basechain = FFG.MarkovChain()
    c2.basechain = FFG.MarkovChain()
    c1.basechain.model = model
    c2.basechain.model = model

    r_up1 = [5,6,4]
    r_down1 = [2,1,3]
    r_up2 = [2,5,1]
    r_down2 = [6,4,3]

    teststate1 = FFG.get_test_state(c1,r_up1,r_down1)
    teststate2 = FFG.get_test_state(c2,r_up2,r_down2)


    dg.replicas[1].basechain.state = teststate1
    dg.replicas[2].basechain.state = teststate2


    # first lets test a region that has a single exchange
    swapsites = [1,4,5]

    (spinc1,spinc2) = (teststate1.spin_config.sc,teststate2.spin_config.sc)

    Tinder1 = FFG.ElectronTinder([],[])
    Tinder2 = FFG.ElectronTinder([],[])
    for site = swapsites

        if spinc1[site] != spinc2[site]
            FFG.addUser!(Tinder1,site,spinc1[site])
            FFG.addUser!(Tinder2,site,spinc2[site])
        end
    end
    @test FFG.SwipeRight(Tinder1) == FFG.SwipeRight(FFG.ElectronTinder([4],[1]))
    @test FFG.SwipeRight(Tinder2) == FFG.SwipeRight(FFG.ElectronTinder([1],[4]))

    # now comes the difficult test...

    #
    r_up1prime = [5,6,1]
    r_down1prime = [2,4,3]
    r_up2prime = [2,5,4]
    r_down2prime = [6,1,3]
    #
    # ts1p = FFG.get_test_state(c1,r_up1prime,r_down1prime)
    # ts2p = FFG.get_test_state(c2,r_up2prime,r_down2prime)

    # print(ts1p.det_A_up * ts1p.det_A_down * ts2p.det_A_up * ts2p.det_A_down,"\n")
    # print(teststate1.det_A_up * teststate1.det_A_down * teststate2.det_A_up * teststate2.det_A_down,"\n")

    # manual_eval = (ts1p.det_A_up * ts1p.det_A_down * ts2p.det_A_up * ts2p.det_A_down)/(teststate1.det_A_up * teststate1.det_A_down * teststate2.det_A_up * teststate2.det_A_down)
    print(FFG.compute_swapregion(dg,swapsites))

end
