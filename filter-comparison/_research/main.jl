
include("install.jl")

include("utils.jl")
include("generate_ground_truth.jl")

params = include("params.jl")

params_gt = params["ground_truth"]
path = datadir("ground_truth")
data_gt, _ = produce_or_load(generate_ground_truth, params_gt, path;
    filename = hash, prefix = "gt", verbose = false, tag = false)

states_gt = data_gt["states"]
observations_gt = data_gt["observations"]
observation_times = data_gt["observation_times"]


include("generate_initial_ensemble.jl")

params_ensemble = params["ensemble"]
params_ensemble["ground_truth"] = params["ground_truth"]
savedir = datadir("ensemble")
filename = string(hash(params_ensemble))
data_initial, _ = produce_or_load(generate_initial_ensemble, params_ensemble, savedir;
    filename, prefix = "initial_ensemble", verbose = false, tag = false)

