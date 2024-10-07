# # Description
# This runs good stuff with EnKF.
#
# # Environment setup
# First, we install the necessary packages.

refresh_cache = false
include("install.jl")

# Import necessary packages.

## Define a macro for doing imports to avoid duplicating it for remote processes later on.
macro initial_imports()
    return esc(
        quote
            using Ensembles
            using LinearAlgebra: norm
            using Distributed:
                addprocs, rmprocs, @everywhere, remotecall, fetch, WorkerPool
            using Test: @test
            using Random: Random
            using CairoMakie

            using Lorenz63: Lorenz63
            ext = Ensembles.get_extension(Ensembles, :Lorenz63Ext)
            using .ext

            using EnsembleKalmanFilters: EnsembleKalmanFilters
            ext = Ensembles.get_extension(Ensembles, :EnsembleKalmanFiltersExt)
            using .ext

            using Statistics: Statistics, mean, var
            using EnsembleKalmanFilters: EnKF
            using Lorenz63Filter
            using PairPlots: PairPlots
            using ImageFiltering: ImageFiltering, imfilter
        end,
    )
end

@initial_imports
worker_initial_imports = @macroexpand1 @initial_imports
include("filter.jl")

# Set parameters.

params = include("params.jl")

# Seed for reproducibility.
Random.seed!(1983745);

# Make operators.
transitioner = Lorenz63Model(; params)
observer = NoisyObserver(get_state_keys(transitioner); params);

# Set seed for ground-truth simulation.
Random.seed!(0xfee55e45)
xor_seed!(observer, UInt64(0x243ecae5));

# Define observation times
observation_times = let
    step = params["observation"]["timestep_size"]
    length = params["observation"]["num_timesteps"]
    range(; start=0, length, step)
end

# Generate synthetic ground-truth observations.
if !(@isdefined ground_truth) || isnothing(ground_truth) || refresh_cache
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
            for (i, t) in enumerate(observation_times)
                state = transitioner(state, t0, t)
                obs = observer(state)
                states[i] = state
                observations[i] = split_clean_noisy(observer, obs)[2]
                t0 = t
            end
        end
        (; states, observations)
    end
    println("  ^ timing for making ground truth observations")
    ground_truth_states_vec = get_ensemble_matrix([:state], ground_truth.states)
    ground_truth_obs_vec = get_ensemble_matrix([:state], ground_truth.observations)
end;

function rmse(ensemble, y_true)
    return sqrt(mean((ensemble .- y_true) .^ 2))
end

# Plot the ground-truth.

if !(@isdefined plot_ground_truth)
    plot_ground_truth = true
end
if plot_ground_truth
    figs = let
        ts = observation_times
        data = reduce(hcat, state[:state] for state in ground_truth.states)

        plot_kwargs = (; color="#7fc97f", marker='.', markersize=15, markercolor=:black)
        plot_state_over_time(ts, data; plot_kwargs...)
    end
    fig = figs[1]
    supertitle = Label(fig[0, :], "ground truth", fontsize = 20, tellwidth=false)
    show_interactive(fig)
end

# Define how to make the initial ensemble.
function generate_ensemble(params::Dict)
    seed = params["ensemble"]["seed"]
    ensemble_size = params["ensemble"]["size"]
    prior_type = params["ensemble"]["prior"]

    members = Vector{Dict{Symbol, Any}}(undef, ensemble_size)
    if prior_type == "gaussian"
        rng = Random.MersenneTwister(seed)
        prior_mean, prior_std = params["ensemble"]["prior_params"]
        for i in 1:ensemble_size
            data = prior_mean .+ prior_std .* randn(rng, 3)
            state = Dict{Symbol, Any}(:state => data)
            members[i] = state
        end
    else
        throw(ArgumentError("Invalid prior type: $prior_type"))
    end

    ensemble = Ensemble(members)
    return ensemble
end;


# Make initial ensemble.
if ! (@isdefined ensemble_initial_spinup) || isnothing(ensemble_initial_spinup) || refresh_cache
    ensemble_initial_spinup = generate_ensemble(params);
end

if !(@isdefined plot_initial_ensemble0)
    plot_initial_ensemble0 = true
end

if plot_initial_ensemble0
    fig = my_pairplot(ensemble_initial_spinup)
    supertitle = Label(fig[0, :], "spinup initial ensemble", fontsize = 20)
    show_interactive(fig)
end

# Run EnKF for 2000 steps just to get a baseline prior.
if ! (@isdefined assimilator_spinup) || isnothing(assimilator_spinup) || refresh_cache
    assimilator_spinup = get_filter(params["spinup"])
end

if ! (@isdefined ensembles_initial) || isnothing(ensembles_initial) || refresh_cache
    ensembles_initial = let t_index_end = params["spinup"]["num_timesteps"],
            observation_times = observation_times[1:t_index_end],
            ground_truth_observations = ground_truth.observations[1:t_index_end],
            ensemble = ensemble_initial_spinup,
            t0 = 0.0,
            transition_noise = params["spinup"]["transition_noise_scale"],
            assimilation_type = params["spinup"]["assimilation_type"],
            assimilator = assimilator_spinup
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
end

# Define some handy functions.

function get_ground_truth_iterator(ensembles_ts, observation_times)
    gt_index = 1
    gt_indices = []
    post_assim_indices = []
    for (i, t) in enumerate(ensembles_ts)
        while gt_index <= length(observation_times) && observation_times[gt_index] < t
            gt_index += 1
        end
        if gt_index > length(observation_times)
            error("Comparing at time $(t) is impossible because final ground-truth observation is at time $(observation_times[end])")
        end
        if observation_times[gt_index] != t
            error("No observation at time $(t). Closest are $(observation_times[gt_index-1]) and $(observation_times[gt_index])")
        end
        push!(gt_indices, gt_index)
        if i == length(ensembles_ts) || t < ensembles_ts[i+1]
            push!(post_assim_indices, i)
        end
    end
    return gt_indices, post_assim_indices
end

function compute_errors(gt_indices, ensembles_means_vec, ground_truth_states_vec)
    rmses = Vector{Float64}(undef, length(ensembles_means_vec))
    for (i, gt_index) in enumerate(gt_indices)
        rmses[i] = rmse(ensembles_means_vec[:, i], ground_truth_states_vec[:, gt_index])
    end
    return rmses
end


function compute_metrics(ensembles)
    means_vec = get_ensemble_matrix(ensembles[1].ensemble.state_keys, mean(e.ensemble) for e in ensembles)
    vars_vec = get_ensemble_matrix(ensembles[1].ensemble.state_keys, var(e.ensemble) for e in ensembles)
    ts = [e.t for e in ensembles]

    gt_indices, post_assim_indices = get_ground_truth_iterator(ts, observation_times)
    rmses = compute_errors(gt_indices, means_vec, ground_truth_states_vec)

    spread = sqrt.(mean(vars_vec; dims=1)[1, :])
    return (;
        ts,
        vars_vec,
        means_vec,
        gt_indices,
        post_assim_indices,
        rmses,
        spread,
        post_assim_rmses = (rmses[i] for i in post_assim_indices),
        post_assim_spread = (spread[i] for i in post_assim_indices),
    )
end;

# Plot ensemble mean error over time.
if !(@isdefined plot_spinup_mean_error)
    plot_spinup_mean_error = true
end
if plot_spinup_mean_error
    figs, ts, data, errors = let ensembles=ensembles_initial
        ts = [e.t for e in ensembles]
        data = get_ensemble_matrix(ensembles[1].ensemble.state_keys, mean(e.ensemble) for e in ensembles)

        errors = fill(NaN, size(data))
        gt_indices, post_assim_indices = get_ground_truth_iterator(ts, observation_times)
        for (i, gt_index) in enumerate(gt_indices)
            errors[:, i] = data[:, i] .- ground_truth_states_vec[:, gt_index]
        end
        plot_kwargs = (; color="#7fc97f", marker='.', markersize=15, markercolor=:black)
        figs = plot_state_over_time(ts, errors; plot_kwargs...)
        figs, ts, data, errors
    end
    fig = figs[1]
    supertitle = Label(fig[0, :], "spinup mean error", fontsize = 20, tellwidth=false)
    show_interactive(fig)
end

# Plot smoothed ensemble mean error over time.
if plot_spinup_mean_error
    figs = let
        smoothed_errors = imfilter(errors, ImageFiltering.Kernel.gaussian((0, length(ts) * 0.01)))
        plot_kwargs = (; color="#7fc97f", marker='.', markersize=15, markercolor=:black)
        figs = plot_state_over_time(ts, smoothed_errors; plot_kwargs...)
    end
    fig = figs[1]
    supertitle = Label(fig[0, :], "spinup smoothed mean error", fontsize = 20, tellwidth=false)
    show_interactive(fig)
end


# Plot smoothed ensemble root mean squared error over time.
if plot_spinup_mean_error
    figs = let
        smoothed_errors = imfilter(errors .^ 2, ImageFiltering.Kernel.gaussian((0, length(ts) * 0.01)))
        smoothed_errors .^= 0.5
        plot_kwargs = (; color="#7fc97f", marker='.', markersize=15, markercolor=:black)
        function handler(fig)
            for c in fig.content
                if isa(c, Axis)
                    ylims!(c; low=0)
                end
            end
        end
        figs = plot_state_over_time(ts, smoothed_errors; handler, plot_kwargs...)
    end
    fig = figs[1]
    supertitle = Label(fig[0, :], "spinup smoothed mean error", fontsize = 20, tellwidth=false)
    show_interactive(fig)
end

# Plot initial ensemble from spinup.
if !(@isdefined plot_initial_ensemble)
    plot_initial_ensemble = true
end

if plot_initial_ensemble
    fig = my_pairplot(ensembles_initial[end].ensemble)
    supertitle = Label(fig[0, :], "initial ensemble", fontsize = 20)
    show_interactive(fig)
end

