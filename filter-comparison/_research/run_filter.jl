
include("install.jl")

using DrWatson: datadir, produce_or_load
using Ensembles
using Random: Random

using Lorenz63: Lorenz63
ext = Ensembles.get_extension(Ensembles, :Lorenz63Ext)
using .ext: Lorenz63Model

include("filter.jl")
include("generate_ground_truth.jl")
include("generate_initial_ensemble.jl")

function run_filter(params::Dict)
    data_gt, _ = produce_or_load_ground_truth(params; loadfile=true)

    data_initial, _ = produce_or_load_initial_ensemble(params; loadfile=true)

    ensemble = data_initial["ensembles"][end].ensemble
    t0 = data_initial["ensembles"][end].t

    states_gt = data_gt["states"]
    observations_gt = data_gt["observations"]
    ts_gt = data_gt["observation_times"]

    params_spinup = params["ensemble"]["spinup"]

    # TODO: check for off-by-one error in the time indexing.
    t0_idx = 1 + params_spinup["num_timesteps"]
    tf_idx = params_spinup["num_timesteps"] + params["estimator"]["num_timesteps"]

    ts_gt = ts_gt[t0_idx:tf_idx]
    observations_gt = observations_gt[t0_idx:tf_idx]

    transition_noise = params["estimator"]["transition_noise_scale"]
    assimilation_type = params["estimator"]["assimilation_type"]
    assimilator = get_filter(params["estimator"])

    transitioner = Lorenz63Model(; params = params["ground_truth"])
    observer = NoisyObserver(get_state_keys(transitioner); params = params["ground_truth"]);

    Random.seed!(0x02cc4823)
    xor_seed!(observer, UInt64(0x54847e5f))

    if assimilation_type == "monolithic"
        observers = [observer]
    elseif assimilation_type == "sequential"
        observer_by_index = function (ensemble, i)
            ensemble = observer(ensemble)
            ## Extract the i-th component from each observation.
            for em in get_ensemble_members(ensemble)
                for key in get_state_keys(observer)
                    em[key] = em[key][i]
                end
            end
            return ensemble
        end
    else
        error("Unknown assimilation type: $(assimilation_type)")
    end

    logs = []
    ensembles = []
    @time begin
        push!(ensembles, (; ensemble, t = t0))
        for (t, y_obs) in zip(ts_gt, observations_gt)
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
            else
                y_obs = [y_obs]
            end
            for (observer_i, y_obs_i) in zip(observers, y_obs)
                ## Take observation at time t.
                ensemble_obs = observer_i(ensemble)
                ensemble_obs_clean, ensemble_obs_noisy = split_clean_noisy(
                    observer, ensemble_obs)

                ## Record.
                push!(ensembles, (; ensemble, ensemble_obs_clean, ensemble_obs_noisy, t))

                ## Assimilate observation
                log_data = Dict{Symbol, Any}()
                (posterior, timing...) = @timed assimilate_data(
                    assimilator, ensemble, ensemble_obs_clean, ensemble_obs_noisy, y_obs_i, log_data)
                log_data[:timing] = timing
                ensemble = posterior

                ## Record.
                push!(ensembles, (; ensemble, t))
                push!(logs, log_data)
            end

            ## Let time pass.
            t0 = t
        end
    end
    println("  ^ timing for running $(params["estimator"]["algorithm"])")

    data = Dict(
        "ensembles" => ensembles,
        "ensemble" => ensembles[end],
    )
end

function produce_or_load_run_filter(params::Dict; kwargs...)
    savedir = datadir("ensemble")
    filename = string(hash(params))
    params_file = joinpath(savedir, filename)
    produce_or_load(run_filter, params, savedir;
        filename, prefix = "run_filter", verbose = false, tag = false,
        loadfile=false, kwargs...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    params = include("params.jl")
    produce_or_load_run_filter(params)
end
