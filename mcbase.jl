# collection of high level basic monte carlo tasks

# high level data structues

abstract type AbstractChain end
abstract type AbstractState end
abstract type AbstractPolicy end
abstract type AbstractMove end
abstract type AbstractObservable end

mutable struct MarkovChain <: AbstractChain
    model::Dict # lattice, temperature, hamiltonian, etc...
    state::AbstractState # spins, fermions, whatever
    policy::AbstractPolicy # random choice, neighbors, ...
    move::AbstractMove # flip a spin, exchange neighbors, ...
    mc_spec::Dict # MC steps, MC warmup, measurement interval
    observable::Any # some function
    data::Dict # outcomes of observable
    diagnostics::Dict # information on the MC run
    MarkovChain() = new()
end

# functions to access parameters of a generic markov chain
get_MarkovChain(c::MarkovChain) = c
get_Model(c::AbstractChain) = get_MarkovChain(c).model
get_State(c::AbstractChain) = get_MarkovChain(c).state
get_Policy(c::AbstractChain) = get_MarkovChain(c).policy
get_Observable(c::AbstractChain) = get_MarkovChain(c).observable
get_Diagnostics(c::AbstractChain) = get_MarkovChain(c).diagnostics
get_Mc_Spec(c::AbstractChain) = get_MarkovChain(c).mc_spec
get_Data(c::AbstractChain) = get_MarkovChain(c).data


""" Run Monte Carlo"""
function runMC!(chain::AbstractChain)
    mc_specs = get_mc_specs(chain)
    NUM_MC_STEPS = mc_specs["num_mc_steps"]
    SAMPLE_INTERVAL = mc_specs["sample_interval"]

    # run burn in steps
    do_warmup!(chain)

    # carry out monte carlo sampling of observable
    for i = 1:NUM_MC_STEPS
        do_move!(chain)
        if mod(i,SAMPLE_INTERVAL) == 0
            measure_observable!(chain)
        end
    end

    # condense measurement statistics and diagnostics from MC sampling
    aggregate_measurements!(chain)
    aggregate_diagnostics!(chain)
end

# measures a generic observable
function measure_observable!(chain::AbstractChain)
    observable = get_Observable(chain)
    result = observable(chain)
    update_data!(chain,result)
end

""" Initialize a markov chain """
function init_Chain!(
    chain::AbstractChain;
    model=model, observable=observable, policy=policy, mc_spec=mc_spec)
    chain.basechain = MarkovChain()
    chain.basechain.model = model
    chain.basechain.state = get_init_state(chain)
    chain.basechain.policy = policy
    chain.basechain.mc_spec = mc_spec
    chain.basechain.data = observable(chain)
    chain.basechain.diagnostics = get_init_diagnostics(chain)
    chain.basechain.observable = observable
end

function do_move!(chain)
    move = get_Move(chain)
    process_move!(chain,move)
end

# recursion for warmup

function do_warmup!(chain)
    mc_specs = get_Mc_Spec(chain)
    num_mc_warmup_steps = mc_specs["mc_warmup_steps"]
    for i = 1:num_mc_warmup_steps
        do_move!(chain)
    end
end




function process_move!(chain,move)
    """ Get the ratio for the proposed move
    and decide whether to accept or reject"""
    ratio, extras = compute_ratio(chain,move) # extras is an optional
    # model based extra return
    ACCEPT =  ratio > rand()
    update_chain!(ACCEPT,chain,move,extras=extras)
end
