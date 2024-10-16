
include("install.jl")

using TerminalLoggers: TerminalLogger
using Logging: global_logger
using ProgressLogging: @progress
global_logger(TerminalLogger())

using DrWatson: srcdir, datadir, produce_or_load
using Ensembles
using Random: Random

using Lorenz63: Lorenz63
ext = Ensembles.get_extension(Ensembles, :Lorenz63Ext)
using .ext: Lorenz63Model

using Distributed: addprocs, @everywhere, WorkerPool, rmprocs, myid

include(srcdir("filter.jl"))
include(srcdir("generate_ground_truth.jl"))

macro worker_imports()
    return esc(
        quote
            using Ensembles
            using Lorenz63: Lorenz63
            ext = Ensembles.get_extension(Ensembles, :Lorenz63Ext)
            using .ext: Lorenz63Model

            using Ensembles: AbstractEnsemble, get_ensemble_members, DistributedOperator

            function (M::IndexObserver)(ensemble::AbstractEnsemble)
                ensemble = M.op(ensemble)
                for em in get_ensemble_members(ensemble)
                    for key in get_state_keys(M.op)
                        em[key] = em[key][M.i]
                    end
                end
                return ensemble
            end


            function (M::DistributedOperator)(member::Dict{Symbol,Any})
                return M.op(member)
            end
        end,
    )
end

@worker_imports
worker_imports = @macroexpand1 @worker_imports

function get_distributed_operator(op, worker_pool, distributed_type, kwargs)
    if distributed_type == :none
        return op
    end
    distributed_type = Val(distributed_type)
    dist_op = DistributedOperator(op, worker_pool, distributed_type; kwargs...)
end

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
    params_estimator = params["spinup"]
    if params_estimator == "none"
        return Dict("ensemble" => ensemble)
    end

    estimator = get_filter(params_estimator)

    data_gt, _ = produce_or_load_ground_truth(params; loadfile=true)
    observations_gt = data_gt["observations"]
    observation_times = data_gt["observation_times"]

    transitioner = Lorenz63Model(; params = params["ground_truth"])
    observer = NoisyObserver(get_state_keys(transitioner); params = params["ground_truth"]);

    t_index_end = params_estimator["num_timesteps"]
    observation_times = observation_times[1:t_index_end]
    observations_gt = observations_gt[1:t_index_end]

    t0 = 0.0

    exec_params = params_estimator["exec"]
    workers = get(exec_params, "workers", 0)
    cleanup_workers = false
    if isa(workers, Int)
        if workers > 0
            workers = addprocs(workers; exeflags="--project=$(Base.active_project())")
            cleanup_workers = true
        else
            workers = ()
        end
    end
    try
        parallel = length(workers) > 0
        if parallel
            worker_pool = WorkerPool(workers)
            @everywhere workers $worker_imports

            t = exec_params["transitioner_distributed_type"]
            kwargs = get(exec_params, "transitioner_distributed_kwargs", ())
            transitioner = get_distributed_operator(transitioner, worker_pool, t, kwargs)

            t = exec_params["observer_distributed_type"]
            kwargs = get(exec_params, "observer_distributed_kwargs", ())
            observer = get_distributed_operator(observer, worker_pool, t, kwargs)
        end
        filter_loop(ensemble, t0, estimator, transitioner, observer, observations_gt, observation_times, params_estimator)
    finally
        if cleanup_workers
            rmprocs(workers)
        end
    end
end

function filter_loop(ensemble, t0, estimator, transitioner, observer, observations_gt, observation_times, params_estimator)
    Random.seed!(0x3289745)
    xor_seed!(observer, UInt64(0x375ef928))

    transition_noise = params_estimator["transition_noise_scale"]
    assimilation_type = params_estimator["assimilation_type"]

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
        @progress "Time: " for (t, y_obs) in zip(observation_times, observations_gt)
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
                    estimator, ensemble, ensemble_obs_clean, ensemble_obs_noisy, y_obs_i, log_data)
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

    data = Dict(
        "ensembles" => ensembles,
        "logs" => logs,
        "ensemble" => ensembles[end],
    )
    return data
end

function initial_ensemble_stem(params::Dict)
    ground_truth_stem(params)*"-"*string(hash(params["ensemble"]), base=62)
end

function produce_or_load_initial_ensemble(params::Dict; kwargs...)
    filestem = initial_ensemble_stem(params)
    params_ensemble = params["ensemble"]
    params_ensemble["ground_truth"] = params["ground_truth"]

    params_file = datadir("initial_ensemble", "params", "$filestem.jld2")
    wsave(params_file, params_ensemble)

    savedir = datadir("initial_ensemble", "data")
    data, filepath = produce_or_load(generate_initial_ensemble, params_ensemble, savedir;
        filename=filestem, verbose = false, tag = false,
        loadfile=false, kwargs...)
    return data, filepath, filestem
end

if abspath(PROGRAM_FILE) == @__FILE__
    params_file = abspath(ARGS[1])
    params = include(params_file)
    produce_or_load_initial_ensemble(params)
end