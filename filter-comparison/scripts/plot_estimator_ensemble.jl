
params_file = abspath(ARGS[1])
include("../src/install.jl")

using DrWatson: srcdir, datadir, plotsdir, produce_or_load, wsave
using CairoMakie: Label, @L_str, Axis, scatterlines!, ylims!, Legend
using Format: cfmt
using Ensembles
using Lorenz63Filter
using ImageFiltering: ImageFiltering, imfilter

include(srcdir("generate_initial_ensemble.jl"))
include(srcdir("run_estimator.jl"))
include(srcdir("utils.jl"))
include(srcdir("plotting.jl"))

# Read data.
params = include(params_file)
@time data_gt, _, filestem_gt = produce_or_load_ground_truth(params; loadfile=true)

@time data_final, _, filestem = produce_or_load_run_estimator(params; loadfile=true)
savedir = plotsdir("estimator", filestem)
if haskey(data_final, "ensembles")
    @time plot_ensemble_data(savedir, data_final["ensembles"], data_gt)
end
