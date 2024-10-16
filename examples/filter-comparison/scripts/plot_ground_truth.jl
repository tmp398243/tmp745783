
params_file = abspath(ARGS[1])
include("../src/install.jl")

using DrWatson: srcdir, datadir, plotsdir, produce_or_load, wsave
using CairoMakie: Label
using Format: cfmt
using Lorenz63Filter

include(srcdir("generate_ground_truth.jl"))
include(srcdir("utils.jl"))

# Read data.
params = include(params_file)
@time data_gt, _, filestem_gt = produce_or_load_ground_truth(params; loadfile=true)

states = data_gt["states"]
observations = data_gt["observations"]
observation_times = data_gt["observation_times"]

figs = let
    ts = observation_times
    matrix = reduce(hcat, state[:state] for state in states)

    plot_kwargs = (; color="#7fc97f", marker='.', markersize=15, markercolor=:black)
    @time plot_state_over_time(ts, matrix; max_dt=50, plot_kwargs...)
end

savedir = plotsdir("ground_truth", "states", filestem_gt)
@time for (i, fig) in enumerate(figs)
    filepath = joinpath(savedir, "$(cfmt("%02d", i)).png")
    wsave(filepath, fig)
end
