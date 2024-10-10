
include("install.jl")

using DrWatson: datadir, plotsdir, produce_or_load, wsave
using CairoMakie: Label, @L_str, Axis, scatterlines!, ylims!, Legend
using Format: cfmt
using Ensembles
using Lorenz63Filter
using ImageFiltering: ImageFiltering, imfilter

include("generate_initial_ensemble.jl")
include("run_estimator.jl")
include("utils.jl")
include("plotting.jl")

# Read data.
params = include("params.jl")
data_gt, _, filestem_gt = produce_or_load_ground_truth(params; loadfile=true)

data_final, _, filestem = produce_or_load_run_estimator(params; loadfile=true)
savedir = plotsdir("estimator", filestem)
if haskey(data_final, "ensembles")
    plot_ensemble_data(savedir, data_final["ensembles"], data_gt)
end
