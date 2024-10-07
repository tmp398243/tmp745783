
include("install.jl")

using DrWatson: datadir, plotsdir, produce_or_load, wsave
using CairoMakie: Label, @L_str, Axis, scatterlines!
using Format: cfmt
using Ensembles
using Lorenz63Filter

include("generate_initial_ensemble.jl")
include("utils.jl")

# Set parameters.
params = include("params.jl")

params_ensemble = params["ensemble"]
params_ensemble["ground_truth"] = params["ground_truth"]
savedir = datadir("ensemble")
filename = string(hash(params_ensemble))
data_initial, _ = produce_or_load(generate_initial_ensemble, params_ensemble, savedir;
    filename, prefix = "initial_ensemble", verbose = false, tag = false)

params_gt = params["ground_truth"]
path = datadir("ground_truth")
data_gt, _ = produce_or_load(generate_ground_truth, params_gt, path; filename = hash, prefix = "ground_truth", verbose = false, tag = false)


savedir = plotsdir("ensemble_initial")
uniquename = string(hash(params_ensemble))

if haskey(data_initial, "ensembles")
    ensembles = data_initial["ensembles"]
    metrics = compute_metrics(ensembles)
    figs = let
        plot_kwargs = (; color="#7fc97f", marker='.', markersize=15, markercolor=:black)
        plot_state_over_time(metrics.ts, metrics.means_vec; max_dt=50, plot_kwargs...)
    end
    fileprefix = joinpath(savedir, "spinup_" * uniquename)
    for (i, fig) in enumerate(figs)
        wsave("$(fileprefix)_$(cfmt("%02d", i)).png", fig) 
    end

    states_gt = data_gt["states"]
    observations_gt = data_gt["observations"]
    ts_gt = data_gt["observation_times"]

    states_vec_gt = get_ensemble_matrix([:state], states_gt)

    t0, tf = extrema(ts)
    start_gt = searchsortedfirst(ts_gt, t0) 
    finish_gt = searchsortedfirst(ts_gt, tf) 

    figs = let
        cut = Colon()

        xs_gt = view(states_vec_gt, 1, :)
        ys_gt = view(states_vec_gt, 2, :)
        zs_gt = view(states_vec_gt, 3, :)

        gt_kwargs = (;
            color = ("#d95f02", 0.5),
            marker = '.',
            markersize = 15,
            markercolor = (:yellow, 0.5),
        )

        handler = function (fig)
            for ax in fig.content
                if isa(ax, Axis)
                    if ax.ylabel[] == L"\text{x}"
                        scatterlines!(ax, ts_gt[start_gt:finish_gt], xs_gt[start_gt:finish_gt]; gt_kwargs...)
                    elseif ax.ylabel[] == L"\text{y}"
                        scatterlines!(ax, ts_gt[start_gt:finish_gt], ys_gt[start_gt:finish_gt]; gt_kwargs...)
                    elseif ax.ylabel[] == L"\text{z}"
                        scatterlines!(ax, ts_gt[start_gt:finish_gt], zs_gt[start_gt:finish_gt]; gt_kwargs...)
                    end
                end
            end
        end
        plot_kwargs = (;
            handler,
            max_dt = 50,
            color = "#7fc97f",
            marker = '.',
            markersize = 0,
            markercolor = :black,
            connect = (;
                linestyle = :dash,
                color = [1, 2],
                colormap = :BuGn,
                markersize = 0,
            ),
        )
        plot_state_over_time(metrics.ts[cut], metrics.means_vec[:, cut]; plot_kwargs...)
    end
    fileprefix = joinpath(savedir, "spinup_gt_" * uniquename)
    for (i, fig) in enumerate(figs)
        wsave("$(fileprefix)_$(cfmt("%02d", i)).png", fig) 
    end
end
