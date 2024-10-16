
params_file = abspath(ARGS[1])
include("../src/install.jl")

using TerminalLoggers: TerminalLogger
using Logging: global_logger
isinteractive() && global_logger(TerminalLogger())
using ProgressLogging: @withprogress, @logprogress

using DrWatson: srcdir, datadir, plotsdir, produce_or_load, wsave
using CairoMakie: Label
using Format: cfmt
using Lorenz63Filter

include(srcdir("generate_ground_truth.jl"))
include(srcdir("utils.jl"))

# Read data.
params = include(params_file)
data_gt, _, filestem_gt = produce_or_load_ground_truth(params; loadfile=true)

states = data_gt["states"]
observations = data_gt["observations"]
observation_times = data_gt["observation_times"]
savedir = plotsdir("ground_truth", "states", filestem_gt)
@withprogress name = "state vs t" let
    ts = observation_times
    matrix = reduce(hcat, state[:state] for state in states)

    plot_kwargs = (; color="#7fc97f", marker='.', markersize=15, markercolor=:black)
    figs = plot_state_over_time(ts, matrix; max_dt=50, plot_kwargs...)

    for (i, fig) in enumerate(figs)
        filepath = joinpath(savedir, "$(cfmt("%02d", i)).png")
        wsave(filepath, fig)
        @logprogress i / length(figs)
    end
end
