
include("install.jl")

using DrWatson: datadir, plotsdir, produce_or_load, wsave
using CairoMakie: Label, @L_str, Axis, scatterlines!, ylims!, Legend
using Format: cfmt
using Ensembles
using Lorenz63Filter
using ImageFiltering: ImageFiltering, imfilter

include("generate_initial_ensemble.jl")
include("utils.jl")
include("plotting.jl")

# Read data.
params = include("params.jl")
data_initial, _ = produce_or_load_initial_ensemble(params; loadfile=true)
data_gt, _ = produce_or_load_ground_truth(params; loadfile=true)

params_ensemble = params["ensemble"]
params_ensemble["ground_truth"] = params["ground_truth"]
savedir = plotsdir("ensemble_initial")
uniquename = string(hash(params_ensemble, hash(params["ground_truth"])))
fileprefix = joinpath(savedir, "spinup_" * uniquename)

if haskey(data_initial, "ensembles")
    plot_ensemble_data(fileprefix, data_initial["ensembles"], data_gt)
end

savedir = plotsdir("estimator")
uniquename = string(hash(params))
fileprefix = joinpath(savedir, "ensemble_" * uniquename)
data_final, _ = produce_or_load_run_filter(params; loadfile=true)
if haskey(data_final, "ensembles")
    plot_ensemble_data(fileprefix, data_final["ensembles"], data_gt)
end
