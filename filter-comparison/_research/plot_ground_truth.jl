
include("install.jl")

using DrWatson: datadir, plotsdir, produce_or_load, wsave
using CairoMakie: Label
using Format: cfmt
using Lorenz63Filter

include("generate_ground_truth.jl")
include("utils.jl")

# Set parameters.
params = include("params.jl")

params_gt = params["ground_truth"]
path = datadir("ground_truth")
data, file = produce_or_load(generate_ground_truth, params_gt, path; filename = hash, prefix = "ground_truth", verbose = false, tag = false)

states = data["states"]
observations = data["observations"]
observation_times = data["observation_times"]

figs = let
    ts = observation_times
    matrix = reduce(hcat, state[:state] for state in states)

    plot_kwargs = (; color="#7fc97f", marker='.', markersize=15, markercolor=:black)
    plot_state_over_time(ts, matrix; max_dt=50, plot_kwargs...)
end
savedir = plotsdir("ground_truth")
uniquename = string(hash(params_gt))
fileprefix = joinpath(savedir, "gt_" * uniquename)
for (i, fig) in enumerate(figs)
    wsave("$(fileprefix)_$(cfmt("%02d", i)).png", fig) 
end
