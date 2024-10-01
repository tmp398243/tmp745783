# # Description
# This example shows a simple use case for Lorenz63Filter.
# Specifically, we use the [`greeting`](@ref) function to print greetings.
#
# # Environment setup
# First, we install the necessary packages.

using Pkg: Pkg

using Lorenz63Filter

try
    using Ensembles
catch
    Pkg.add(; url="https://github.com/tmp398243/tmp32487543")
    using Ensembles
end

try
    using Lorenz63: Lorenz63
catch
    Ensembles.install(:Lorenz63)
end

try
    using EnsembleKalmanFilters: EnsembleKalmanFilters
catch
    Ensembles.install(:EnsembleKalmanFilters)
end

try
    using NormalizingFlowFilters: NormalizingFlowFilters
catch
    Ensembles.install(:NormalizingFlowFilters)
end

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
            using LinearAlgebra: Diagonal
            using EnsembleKalmanFilters: EnKF
        end,
    )
end

@initial_imports
worker_initial_imports = @macroexpand1 @initial_imports

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
if !(@isdefined ground_truth) || isnothing(ground_truth)
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
        plot_state_over_time(observation_times, data; plot_kwargs...)
    end
    figs[1]
end
