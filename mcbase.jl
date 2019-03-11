# collection of high level basic monte carlo tasks

# high level data structues

abstract type AbstractChain end
abstract type AbstractState end
abstract type AbstractPolicy end
abstract type AbstractMove end
abstract type AbstractObservable end

mutable struct BaseChain <: AbstractChain
    model::Dict # lattice, temperature, hamiltonian, etc...
    state::AbstractState # spins, fermions, whatever
    policy::AbstractPolicy # random choice, neighbors, ...
    move::AbstractMove # flip a spin, exchange neighbors, ...
    mc_spec::Dict # MC steps, MC warmup, measurement interval
    observable::AbstractObservable # some function
    data::Dict # outcomes of observable
    diagnostics::Dict # information on the MC run
    BaseChain() = new()
end




function runMC(chain::AbstractChain)
    mc_specs = get_mc_specs(chain)
    NUM_MC_STEPS = mc_specs["num_mc_steps"]
    NUM_MC_WARMUP_STEPS = mc_specs["num_mc_warmup_steps"]
    SAMPLE_INTERVAL = mc_specs["sample_interval"]

    # run burn in steps
    do_warmup(chain,NUM_MC_WARMUP_STEPS)

    # carry out monte carlo sampling of observable
    for i = 1:NUM_MC_STEPS
        do_move(chain)
        if mod(i,SAMPLE_INTERVAL) == 0
            measure_observable!(chain)
        end
    end

    # condense measurement statistics and diagnostics from MC sampling
    aggregate_measurements!(chain)
    aggregate_diagnostics!(chain)
end

function measure_observable!(chain::AbstractChain)

end

function init_Chain(chain::AbstractChain
    model=model, observable=observable, policy=policy, mc_specs=mc_specs)
    
end
