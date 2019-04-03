using Test
include("FFGmain.jl")


ffgmathtest = @testset "test-FFGMATH" begin
    ############### testing lattice routines #################################
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
    move = FFG.SwapNeighborMove((1,2))

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
