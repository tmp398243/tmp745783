
include("install.jl")

using DrWatson: datadir, produce_or_load
using Ensembles
using Random: Random

using Lorenz63: Lorenz63
ext = Ensembles.get_extension(Ensembles, :Lorenz63Ext)
using .ext: Lorenz63Model

include("filter.jl")
include("generate_ground_truth.jl")

function generate_initial_ensemble(params::Dict)
    seed = params["seed"]
    ensemble_size = params["size"]
    prior_type = params["prior"]

    members = Vector{Dict{Symbol, Any}}(undef, ensemble_size)
    if prior_type == "gaussian"
        rng = Random.MersenneTwister(seed)
        prior_mean, prior_std = params["prior_params"]
        for i in 1:ensemble_size
            data = prior_mean .+ prior_std .* randn(rng, 3)
            state = Dict{Symbol, Any}(:state => data)
            members[i] = state
        end
    else
        throw(ArgumentError("Invalid prior type: $prior_type"))
    end

    ensemble = Ensemble(members)
    if params["spinup"] == "none"
        return Dict("ensemble" => ensemble)
    end

    assimilator = get_filter(params["spinup"])

    params_gt = params["ground_truth"]
    path = datadir("ground_truth")
    data_gt, _ = produce_or_load(generate_ground_truth, params_gt, path;
        filename = hash, prefix = "ground_truth", verbose = false, tag = false)

    observations_gt = data_gt["observations"]
    observation_times = data_gt["observation_times"]

    transitioner = Lorenz63Model(; params = params["ground_truth"])
    observer = NoisyObserver(get_state_keys(transitioner); params = params["ground_truth"]);

    ensembles = let t_index_end = params["spinup"]["num_timesteps"],
            observation_times = observation_times[1:t_index_end],
            ground_truth_observations = observations_gt[1:t_index_end],
            t0 = 0.0,
            transition_noise = params["spinup"]["transition_noise_scale"],
            assimilation_type = params["spinup"]["assimilation_type"]

        Random.seed!(0x3289745)
        xor_seed!(observer, UInt64(0x375ef928))

        if assimilation_type == "monolithic"
            observers = [observer]
        elseif assimilation_type == "sequential"
        else
            error("Unknown assimilation type: $(assimilation_type)")
        end

        logs = []
        ensembles = []
        @time begin
            push!(ensembles, (; ensemble, t = t0))
            for (t, y_obs) in zip(observation_times, ground_truth_observations)
                ## Advance ensemble to time t.
                ensemble = transitioner(ensemble, t0, t; inplace = false)

                ## Keep ensemble separated.
                if transition_noise != 0
                    for em in ensemble.members
                        em[:state] .+= transition_noise .* Random.randn(3)
                    end
                end

                if assimilation_type == "sequential"
                    y_obs_vec = get_member_vector(ensemble, y_obs)
                    observers = [IndexObserver(observer, i) for i in 1:length(y_obs_vec)]
                    y_obs = [observer_i(y_obs) for observer_i in observers]
                    push!(ensembles, (; ensemble, t))
                else
                    y_obs = [y_obs]
                end
                for (observer_i, y_obs_i) in zip(observers, y_obs)
                    ## Take observation at time t.
                    ensemble_obs = observer_i(ensemble)
                    ensemble_obs_clean, ensemble_obs_noisy = split_clean_noisy(
                        observer, ensemble_obs)

                    ## Record.
                    if assimilation_type != "sequential"
                        push!(ensembles, (; ensemble, ensemble_obs_clean, ensemble_obs_noisy, t))
                    end

                    ## Assimilate observation
                    log_data = Dict{Symbol, Any}()
                    (posterior, timing...) = @timed assimilate_data(
                        assimilator, ensemble, ensemble_obs_clean, ensemble_obs_noisy, y_obs_i, log_data)
                    log_data[:timing] = timing
                    ensemble = posterior

                    push!(logs, log_data)
                end
                ## Record.
                push!(ensembles, (; ensemble, t))

                ## Let time pass.
                t0 = t
            end
        end
        println("  ^ timing for making initial ensemble")
        ensembles
    end

    data = Dict(
        "ensembles" => ensembles,
        "ensemble" => ensembles[end],
    )
end

function produce_or_load_initial_ensemble(params::Dict; kwargs...)
    params_ensemble = params["ensemble"]
    params_ensemble["ground_truth"] = params["ground_truth"]
    savedir = datadir("ensemble")
    filename = string(hash(params_ensemble))
    produce_or_load(generate_initial_ensemble, params_ensemble, savedir;
        filename, prefix = "initial_ensemble", verbose = false, tag = false,
        loadfile=false, kwargs...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    params = include("params.jl")
    produce_or_load_initial_ensemble(params)
end