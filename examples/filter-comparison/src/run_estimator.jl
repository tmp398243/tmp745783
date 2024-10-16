
include("install.jl")

using TerminalLoggers: TerminalLogger
using Logging: global_logger
isinteractive() && global_logger(TerminalLogger())

using DrWatson: srcdir, datadir, produce_or_load, wsave
using Ensembles
using Random: Random

using Lorenz63: Lorenz63
ext = Ensembles.get_extension(Ensembles, :Lorenz63Ext)
using .ext: Lorenz63Model

include(srcdir("filter.jl"))
include(srcdir("generate_ground_truth.jl"))
include(srcdir("generate_initial_ensemble.jl"))
include(srcdir("filter_loop.jl"))

function run_estimator(params::Dict)
    data_gt, _ = produce_or_load_ground_truth(params; loadfile=true)

    data_initial, _ = produce_or_load_initial_ensemble(params; loadfile=true)

    states_gt = data_gt["states"]
    observations_gt = data_gt["observations"]
    ts_gt = data_gt["observation_times"]

    params_spinup = params["ensemble"]["spinup"]

    # TODO: check for off-by-one error in the time indexing.
    t0_idx = 1 + params_spinup["num_timesteps"]
    tf_idx = params_spinup["num_timesteps"] + params["estimator"]["num_timesteps"]

    ts_gt = ts_gt[t0_idx:tf_idx]
    observations_gt = observations_gt[t0_idx:tf_idx]

    params_estimator = params["estimator"]
    estimator = get_filter(params_estimator)

    transitioner = Lorenz63Model(; params=params["ground_truth"])
    observer = NoisyObserver(get_state_keys(transitioner); params=params["ground_truth"])

    ensemble = data_initial["ensemble"].ensemble
    t0 = data_initial["ensemble"].t

    Random.seed!(0x02cc4823)
    xor_seed!(observer, UInt64(0x54847e5f))

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
        filter_loop(
            ensemble,
            t0,
            estimator,
            transitioner,
            observer,
            observations_gt,
            ts_gt,
            params_estimator;
            name=params_estimator["algorithm"],
        )
    finally
        if cleanup_workers
            rmprocs(workers)
        end
    end
end

function filter_stem(params::Dict)
    # TODO: shouldn't include "exec" params in this hash.
    return ground_truth_stem(params) *
           "-" *
           initial_ensemble_stem(params) *
           "-" *
           string(hash(params); base=62)
end

function produce_or_load_run_estimator(params::Dict; kwargs...)
    filestem = filter_stem(params)

    params_file = datadir("estimator", "params", "$filestem.jld2")
    wsave(params_file, params)

    savedir = datadir("estimator", "data")
    data, filepath = produce_or_load(
        run_estimator,
        params,
        savedir;
        filename=filestem,
        verbose=false,
        tag=false,
        loadfile=false,
        kwargs...,
    )
    return data, filepath, filestem
end

if abspath(PROGRAM_FILE) == @__FILE__
    params_file = abspath(ARGS[1])
    params = include(params_file)
    produce_or_load_run_estimator(params)
end
