
include("install.jl")

using DrWatson: datadir, plotsdir, produce_or_load, wsave
using CairoMakie: Label
using Format: cfmt
using Lorenz63Filter

include("generate_ground_truth.jl")
include("utils.jl")

# Read data.
params = include("params.jl")
data_gt, _, filestem_gt = produce_or_load_ground_truth(params; loadfile=true)

states = data_gt["states"]
observations = data_gt["observations"]
observation_times = data_gt["observation_times"]

figs = let
    ts = observation_times
    matrix = reduce(hcat, state[:state] for state in states)

    plot_kwargs = (; color="#7fc97f", marker='.', markersize=15, markercolor=:black)
    plot_state_over_time(ts, matrix; max_dt=50, plot_kwargs...)
end

params_gt = params["ground_truth"]
savedir = plotsdir("ground_truth", "states", filestem_gt)
fileprefix = joinpath(filepath)
for (i, fig) in enumerate(figs)
    filepath = joinpath(savedir, "$(cfmt("%02d", i)).png")
    wsave(filepath, fig) 
end
