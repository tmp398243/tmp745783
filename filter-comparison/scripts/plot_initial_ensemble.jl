
params_file = abspath(ARGS[1])
include("../src/install.jl")

using TerminalLoggers: TerminalLogger
using Logging: global_logger
global_logger(TerminalLogger())

using DrWatson: srcdir, datadir, plotsdir, produce_or_load, wsave
using CairoMakie: Label, @L_str, Axis, scatterlines!, ylims!, Legend
using Format: cfmt
using Ensembles
using Lorenz63Filter
using ImageFiltering: ImageFiltering, imfilter

include(srcdir("generate_initial_ensemble.jl"))
include(srcdir("utils.jl"))
include(srcdir("plotting.jl"))

# Read data.
params = include(params_file)
@time data_initial, _, filestem = produce_or_load_initial_ensemble(params; loadfile=true)
@time data_gt, _, filestem_gt = produce_or_load_ground_truth(params; loadfile=true)

savedir = plotsdir("ensemble_initial", filestem)
if haskey(data_initial, "ensembles")
    @time plot_ensemble_data(savedir, data_initial["ensembles"], data_gt)
end
