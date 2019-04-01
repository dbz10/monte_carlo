using Test
include("FreeFermionGutzwiller.jl")


@testset "test-FFGMATH" begin
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
    @test FreeFermionGutzwiller.get_Determinant(r,states) ≈ 0.3333333333
    tmat = [0 -1.2247448713915894 -1.2247448713915894;
            -1.2247448713915894 1.2247448713915894 0;
            1 0 -1]
    @test FreeFermionGutzwiller.get_Inverse_Matrix(r,states) ≈ tmat

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

    @test FreeFermionGutzwiller.det_ratio_factor(m1inv,u,v) ≈ det(m2)/det(m1)

    @test FreeFermionGutzwiller.sherman_morrison_inverse_update(m1inv,u,v) ≈ inv(m2)
end

@testset "test-FFGMAIN" begin
    ############ make a test chain $$$$$$$$$$$$$$$$$$$$$$$$

    model = Dict(
        "dims" => (2,3),
        "lattice" => Lattices.get_SquareLattice((2,3),pbc=true),
        "filling" => 1,
        "hamiltonian" => Lattices.get_Tightbinding_Wavefunctions,
        "fermi_energy" => 0,
        )

    function observable(chain::FreeFermionGutzwiller.GutzwillerChain)
        state = FreeFermionGutzwiller.get_State(chain)
        Sz1 = state.spin_config.sc[1]
        return Dict("Sz1" => Sz1)
    end

    policy = FreeFermionGutzwiller.SwapNeighborsPolicy()
    mc_spec = Dict("mc_steps" => 1,
                    "mc_warmup_steps" => 1)

    testchain = FreeFermionGutzwiller.GutzwillerChain();
    FreeFermionGutzwiller.init_Chain!(testchain,
            model=model,observable=observable,policy=policy,mc_spec=mc_spec
            )

    r_up = [3,1,5]
    r_down = [2,6,4]

    # some trivial stuff
    bd = [2, 4, 1, 6, 3, 5]
    @test FreeFermionGutzwiller.make_business_directory(r_up,r_down) == bd

    spinc = FreeFermionGutzwiller.SpinConfiguration([1,-1,1,-1,1,-1])
    @test FreeFermionGutzwiller.make_spin_config(r_up,r_down).sc == spinc.sc

    # testing construction of initial bond graph
    using LightGraphs: SimpleGraph, add_edge!
    ib = SimpleGraph(6);
    add_edge!(ib,1,2);
    add_edge!(ib,3,4);
    add_edge!(ib,5,6);

    @test FreeFermionGutzwiller.get_init_bonds(model["lattice"].graph,spinc) ==
    ib



end
