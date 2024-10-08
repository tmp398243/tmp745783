
include("install.jl")

include("utils.jl")
include("generate_ground_truth.jl")
include("generate_initial_ensemble.jl")

params = include("params.jl")

data_gt, _ = produce_or_load_ground_truth(params; loadfile=true)

states_gt = data_gt["states"]
observations_gt = data_gt["observations"]
observation_times = data_gt["observation_times"]

data_initial, _ = produce_or_load_initial_ensemble(params; loadfile=true)


# Make a filter.

# Run the filter.