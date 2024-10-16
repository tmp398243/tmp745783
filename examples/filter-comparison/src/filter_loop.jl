
using Random: Random
using ProgressLogging: @progress

function filter_loop(
    ensemble,
    t0,
    estimator,
    transitioner,
    observer,
    observations_gt,
    observation_times,
    params_estimator;
    name="Time",
)
    Random.seed!(0x3289745)
    xor_seed!(observer, UInt64(0x375ef928))

    transition_noise = params_estimator["transition_noise_scale"]
    assimilation_type = params_estimator["assimilation_type"]

    if assimilation_type == "monolithic"
        observers = [observer]
    elseif assimilation_type == "sequential"
        y_obs_vec = get_member_vector(ensemble, observations_gt[1])
        observers = [IndexObserver(observer, i) for i in 1:length(y_obs_vec)]
    else
        error("Unknown assimilation type: $(assimilation_type)")
    end

    logs = []
    ensembles = []
    progress_name = name * ": "
    @time begin
        push!(ensembles, (; ensemble, t=t0))
        @progress name = progress_name for (t, y_obs) in
                                           zip(observation_times, observations_gt)
            ## Advance ensemble to time t.
            ensemble = transitioner(ensemble, t0, t; inplace=false)

            ## Keep ensemble separated.
            if transition_noise != 0
                for em in ensemble.members
                    em[:state] .+= transition_noise .* Random.randn(3)
                end
            end

            if assimilation_type == "sequential"
                y_obs_vec = get_member_vector(ensemble, y_obs)
                append!(
                    observers,
                    IndexObserver(observer, i) for i in length(observers):length(y_obs_vec)
                )
                y_obs = [observer_i(y_obs) for observer_i in observers[1:length(y_obs_vec)]]
                push!(ensembles, (; ensemble, t))
            else
                y_obs = (y_obs,)
            end
            for (observer_i, y_obs_i) in zip(observers, y_obs)
                ## Take observation at time t.
                ensemble_obs = observer_i(ensemble)
                ensemble_obs_clean, ensemble_obs_noisy = split_clean_noisy(
                    observer, ensemble_obs
                )

                ## Record.
                if assimilation_type != "sequential"
                    push!(
                        ensembles, (; ensemble, ensemble_obs_clean, ensemble_obs_noisy, t)
                    )
                end

                ## Assimilate observation
                log_data = Dict{Symbol,Any}()
                (posterior, timing...) = @timed assimilate_data(
                    estimator,
                    ensemble,
                    ensemble_obs_clean,
                    ensemble_obs_noisy,
                    y_obs_i,
                    log_data,
                )
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
    println("  ^ timing for running filter loop ($name)")

    data = Dict("ensembles" => ensembles, "logs" => logs, "ensemble" => ensembles[end])
    return data
end
