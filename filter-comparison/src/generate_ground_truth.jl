
include("install.jl")

using DrWatson: datadir, produce_or_load
using Ensembles: Ensembles, NoisyObserver, get_state_keys, get_ensemble_matrix, split_clean_noisy, xor_seed!
using Random: Random

using Lorenz63: Lorenz63
ext = Ensembles.get_extension(Ensembles, :Lorenz63Ext)
using .ext: Lorenz63Model


# Generate synthetic ground-truth observations.
function generate_ground_truth(params::Dict)    
    observation_times = let
        step = params["observation"]["timestep_size"]
        length = params["observation"]["num_timesteps"] + 1
        range(; start=0, length, step)
    end

    ## Make operators.
    transitioner = Lorenz63Model(; params)
    observer = NoisyObserver(get_state_keys(transitioner); params);

    ## Set seed for ground-truth simulation.
    Random.seed!(0xfee55e45)
    xor_seed!(observer, UInt64(0x243ecae5));

    ground_truth = @time let
        state0 = Dict{Symbol,Any}(:state => randn(3))

        ## Set seed for ground-truth simulation.
        Random.seed!(0xfee55e45)
        xor_seed!(observer, UInt64(0x243ecae5))

        ## Generate states and observations.
        t0 = 0.0
        states = Vector{Dict{Symbol,Any}}(undef, length(observation_times))
        observations = Vector{Dict{Symbol,Any}}(undef, length(observation_times))
        let state = state0
            states[1] = state0
            observations[1] = observer(state0)
            for (i, t) in enumerate(observation_times[2:end])
                state = transitioner(state, t0, t)
                obs = observer(state)
                states[i+1] = state
                observations[i+1] = split_clean_noisy(observer, obs)[2]
                t0 = t
            end
        end
        (; states, observations)
    end
    println("  ^ timing for making ground truth data")
    data = Dict(
        "states" => ground_truth.states,
        "observations" => ground_truth.observations,
        "observation_times" => observation_times,
    )
end

function ground_truth_stem(params::Dict)
    string(hash(params["ground_truth"]), base=62)
end

function produce_or_load_ground_truth(params::Dict; kwargs...)
    params_gt = params["ground_truth"]
    filestem = ground_truth_stem(params)

    params_file = datadir("ground_truth", "params", "$filestem.jld2")
    wsave(params_file, params_gt)

    savedir = datadir("ground_truth", "data")
    data, filepath = produce_or_load(generate_ground_truth, params_gt, savedir;
        filename=filestem, verbose = false, tag = false,
        loadfile=false, kwargs...)
    return data, filepath, filestem
end

if abspath(PROGRAM_FILE) == @__FILE__
    params = include(ARGS[1])
    produce_or_load_ground_truth(params)
end