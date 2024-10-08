
include("install.jl")

using DrWatson: datadir, plotsdir, produce_or_load, wsave
using CairoMakie: Label, @L_str, Axis, scatterlines!, ylims!, Legend
using Format: cfmt
using Ensembles
using Lorenz63Filter
using ImageFiltering: ImageFiltering, imfilter

include("generate_initial_ensemble.jl")
include("utils.jl")

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
    ensembles = data_initial["ensembles"]
    states_gt = data_gt["states"]
    observations_gt = data_gt["observations"]
    ts_gt = data_gt["observation_times"]
    ground_truth_states_vec = get_ensemble_matrix([:state], states_gt)

    metrics = compute_metrics(ensembles; ts_gt, ground_truth_states_vec)
    N_states = length(metrics.ts)

    # Plot state over time.
    figs = let
        plot_kwargs = (; disjoint=false, color="#7fc97f", marker='.', markersize=15, markercolor=:black)
        plot_state_over_time(metrics.ts, metrics.means_vec; max_dt=50, plot_kwargs...)
    end
    for (i, fig) in enumerate(figs)
        wsave("$(fileprefix)_$(cfmt("%02d", i)).png", fig) 
    end

    t0, tf = extrema(metrics.ts)
    start_gt = searchsortedfirst(ts_gt, t0) 
    finish_gt = searchsortedfirst(ts_gt, tf) 

    # Plot state over time with ground truth state.
    figs = let
        xs_gt = view(ground_truth_states_vec, 1, :)
        ys_gt = view(ground_truth_states_vec, 2, :)
        zs_gt = view(ground_truth_states_vec, 3, :)

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
            disjoint=false,
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
        plot_state_over_time(metrics.ts, metrics.means_vec; plot_kwargs...)
    end
    for (i, fig) in enumerate(figs)
        wsave("$(fileprefix)_gt_$(cfmt("%02d", i)).png", fig) 
    end

    # Plot ensemble mean error over time.
    errors = fill(NaN, size(metrics.means_vec))
    for (i, gt_index) in enumerate(metrics.gt_indices)
        errors[:, i] = metrics.means_vec[:, i] .- ground_truth_states_vec[:, gt_index]
    end
    figs = let
        plot_kwargs = (;
            disjoint=false,
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
        figs = plot_state_over_time(metrics.ts, errors; plot_kwargs...)
    end
    for (i, fig) in enumerate(figs)
        wsave("$(fileprefix)_gt_error_$(cfmt("%02d", i)).png", fig) 
    end

    # Plot smoothed ensemble mean error over time.
    figs = let
        smoothed_errors = imfilter(errors, ImageFiltering.Kernel.gaussian((0, N_states * 0.01)))
        plot_kwargs = (;
            disjoint=false,
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
        figs = plot_state_over_time(metrics.ts, smoothed_errors; plot_kwargs...)
    end
    for (i, fig) in enumerate(figs)
        wsave("$(fileprefix)_gt_error_smoothed_$(cfmt("%02d", i)).png", fig) 
    end

    # Plot smoothed ensemble root mean squared error over time.
    figs = let
        smoothed_errors = imfilter(errors .^ 2, ImageFiltering.Kernel.gaussian((0, N_states * 0.01)))
        smoothed_errors .^= 0.5
        plot_kwargs = (;
            disjoint=false,
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
        function handler(fig)
            for c in fig.content
                if isa(c, Axis)
                    ylims!(c; low=0)
                end
            end
        end
        figs = plot_state_over_time(metrics.ts, smoothed_errors; handler, plot_kwargs...)
    end
    for (i, fig) in enumerate(figs)
        wsave("$(fileprefix)_gt_abserror_smoothed_$(cfmt("%02d", i)).png", fig) 
    end

    # Plot initial and final ensemble from spinup.
    fig = my_pairplot(ensembles[1].ensemble)
    wsave("$(fileprefix)_initial.png", fig) 

    fig = my_pairplot(ensembles[end].ensemble)
    wsave("$(fileprefix)_final.png", fig)

    # Plot standard deviation.
    let
        handler = function (fig)
            for ax in fig.content
                if isa(ax, Axis)
                    if ax.ylabel[] == L"\text{x}"
                        ax.ylabel = L"\sigma_x"
                        # ax.yscale = log10
                    elseif ax.ylabel[] == L"\text{y}"
                        ax.ylabel = L"\sigma_y"
                        # ax.yscale = log10
                    elseif ax.ylabel[] == L"\text{z}"
                        ax.ylabel = L"\sigma_z"
                        # ax.yscale = log10
                    end
                end
            end
        end
        plot_kwargs = (;
            handler,
            disjoint=false,
            max_dt = 50,
            make_positive = true,
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
        figs = plot_state_over_time(metrics.ts, metrics.vars_vec .^ 0.5; plot_kwargs...)
        for (i, fig) in enumerate(figs)
            wsave("$(fileprefix)_gt_std_$(cfmt("%02d", i)).png", fig) 
        end

        smoothed_stds = imfilter(metrics.vars_vec, ImageFiltering.Kernel.gaussian((0, N_states * 0.01))) .^ 0.5
        figs = plot_state_over_time(metrics.ts, smoothed_stds; plot_kwargs...)
        for (i, fig) in enumerate(figs)
            wsave("$(fileprefix)_gt_std_smoothed_$(cfmt("%02d", i)).png", fig) 
        end
    end

    # Plot ensemble spread.
    let
        handler = function (fig)
            for ax in fig.content
                if isa(ax, Axis)
                    if ax.ylabel[] == L"\text{metric}"
                        ax.ylabel = L"\text{spread}"
                        ax.yscale = log10
                    end
                end
            end
        end
        plot_kwargs = (;
            handler,
            disjoint=false,
            max_dt = 50,
            color = "#041a1c",
            marker = '.',
            markersize = 15,
            markercolor = :black,
            connect = (;
                linestyle = :dash,
                color = [1, 2],
                colormap = :BuGn,
                markersize = 0,
            ),
        )
        figs = plot_error_metric_over_time(metrics.ts, metrics.spread; plot_kwargs...)
        for (i, fig) in enumerate(figs)
            wsave("$(fileprefix)_gt_spread_$(cfmt("%02d", i)).png", fig) 
        end

        smoothed_spread = imfilter(metrics.spread, ImageFiltering.Kernel.gaussian((N_states * 0.01,)))
        figs = plot_error_metric_over_time(metrics.ts, smoothed_spread; plot_kwargs...)
        for (i, fig) in enumerate(figs)
            wsave("$(fileprefix)_gt_spread_smoothed_$(cfmt("%02d", i)).png", fig) 
        end
    end

    # Plot RMSE.
    let
        handler = function (fig)
            for ax in fig.content
                if isa(ax, Axis)
                    if ax.ylabel[] == L"\text{metric}"
                        ax.ylabel = L"\text{RMSE}"
                        ax.yscale = log10
                    end
                end
            end
        end
        plot_kwargs = (;
            handler,
            disjoint=false,
            max_dt = 50,
            color = "#e41a1c",
            marker = '.',
            markersize = 15,
            markercolor = :black,
            connect = (;
                linestyle = :dash,
                color = [1, 2],
                colormap = :BuGn,
                markersize = 0,
            ),
        )
        figs = plot_error_metric_over_time(metrics.ts, metrics.rmses; plot_kwargs...)
        for (i, fig) in enumerate(figs)
            wsave("$(fileprefix)_gt_rmse_$(cfmt("%02d", i)).png", fig) 
        end

        smoothed_rmses = imfilter(metrics.rmses .^ 2, ImageFiltering.Kernel.gaussian((N_states * 0.01,))) .^ 0.5
        figs = plot_error_metric_over_time(metrics.ts, smoothed_rmses; plot_kwargs...)
        for (i, fig) in enumerate(figs)
            wsave("$(fileprefix)_gt_rmse_smoothed_$(cfmt("%02d", i)).png", fig) 
        end
    end


    # Plot ensemble spread and RMSE.
    let
        handler = function (fig)
            for ax in fig.content
                if isa(ax, Axis)
                    if ax.ylabel[] == L"\text{metric}"
                        ax.yscale = log10
                    end
                end
            end
        end
        plot_kwargs_spread = (;
            disjoint=false,
            color = "#041a1c",
            marker = '.',
            markersize = 15,
            markercolor = :black,
            connect = (;
                linestyle = :dash,
                color = [1, 2],
                colormap = :BuGn,
                markersize = 0,
            ),
        )
        plot_kwargs_full = (;
            handler,
            max_dt = 50,
            plot_kwargs_spread...
        )
        figs = plot_error_metric_over_time(metrics.ts, metrics.spread; label="spread", plot_kwargs_full...)

        plot_kwargs_rmse = (;
            plot_kwargs_spread...,
            color = "#e41a1c"
        )
        for (i, fig) in enumerate(figs)
            ax = fig.content[1]
            plot_disjoint_lines!(ax, metrics.ts, metrics.rmses; label="rmse", plot_kwargs_rmse...)
            fig[1, 2] = Legend(fig, ax; unique=true)
            wsave("$(fileprefix)_gt_spread_rmse_$(cfmt("%02d", i)).png", fig) 
        end

        # Plot smoothed values.
        smoothed_spread = imfilter(metrics.spread, ImageFiltering.Kernel.gaussian((N_states * 0.01,)))
        figs = plot_error_metric_over_time(metrics.ts, smoothed_spread; label="spread", plot_kwargs_full...)

        smoothed_rmses = imfilter(metrics.rmses .^ 2, ImageFiltering.Kernel.gaussian((N_states * 0.01,))) .^ 0.5
        for (i, fig) in enumerate(figs)
            ax = fig.content[1]
            plot_disjoint_lines!(ax, metrics.ts, smoothed_rmses; label="rmse", plot_kwargs_rmse...)
            fig[1, 2] = Legend(fig, ax; unique=true)
            wsave("$(fileprefix)_gt_spread_rmse_smoothed_$(cfmt("%02d", i)).png", fig) 
        end
    end
end
