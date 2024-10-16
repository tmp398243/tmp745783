
# First, generate and plot the synthetic ground truth.
using Pkg

ENV["lorenz63filter_force_install"] = "true"

params_file = abspath(joinpath(@__DIR__, "small-params.jl"))
push!(ARGS, params_file)

include("scripts/plot_ground_truth.jl")

using Markdown

fig_path = relpath(joinpath(savedir, "01.png"))

Markdown.parse("![ground truth states]($fig_path)")

# Then, we can do other stuff.

include("scripts/plot_initial_ensemble.jl")

# Here's the initial ensemble.

fig_path = relpath(joinpath(savedir, "state", "initial.png"))
Markdown.parse("![initial state]($fig_path)")

# Here's the mean estimate over time.
fig_path = relpath(joinpath(savedir, "mean_state", "01.png"))
Markdown.parse("![mean state]($fig_path)")

# Here's the mean estimate over time with the ground truth.
fig_path = relpath(joinpath(savedir, "mean_state_gt", "01.png"))
Markdown.parse("![mean state with ground truth]($fig_path)")

# Here's the spread over time.
fig_path = relpath(joinpath(savedir, "spread", "01.png"))
Markdown.parse("![spread]($fig_path)")

# Here's the spread smoothed over time.
fig_path = relpath(joinpath(savedir, "spread_smoothed", "01.png"))
Markdown.parse("![spread smoothed]($fig_path)")

# Here's the standard deviation in each coordinate smoothed over time.
fig_path = relpath(joinpath(savedir, "std_smoothed", "01.png"))
Markdown.parse("![std smoothed]($fig_path)")

# Here's the absolute error of the mean in each coordinate smoothed over time.
fig_path = relpath(joinpath(savedir, "mean_abserror_smoothed", "01.png"))
Markdown.parse("![mean abserror smoothed]($fig_path)")

# Here's the rmse over time.
fig_path = relpath(joinpath(savedir, "mean_rmse", "01.png"))
Markdown.parse("![rmse]($fig_path)")

# Here's the rmse smoothed over time.
fig_path = relpath(joinpath(savedir, "mean_rmse_smoothed", "01.png"))
Markdown.parse("![rmse smoothed]($fig_path)")

# Here's the spread and rmse smoothed over time.
fig_path = relpath(joinpath(savedir, "spread_mean_rmse_smoothed", "01.png"))
Markdown.parse("![spread rmse smoothed]($fig_path)")
