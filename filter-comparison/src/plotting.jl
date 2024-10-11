
using DrWatson: srcdir, datadir, plotsdir, produce_or_load, wsave
using CairoMakie: Label, @L_str, Axis, scatterlines!, ylims!, Legend
using Format: cfmt
using Ensembles
using Lorenz63Filter
using ImageFiltering: ImageFiltering, imfilter

include(srcdir("utils.jl"))

function plot_ensemble_data(savedir_root, ensembles, data_gt)
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

    savedir = joinpath(savedir_root, "mean_state")
    for (i, fig) in enumerate(figs)
        filepath = joinpath(savedir, "$(cfmt("%02d", i)).png")
        wsave(filepath, fig)
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

    savedir = joinpath(savedir_root, "mean_state_gt")
    for (i, fig) in enumerate(figs)
        filepath = joinpath(savedir, "$(cfmt("%02d", i)).png")
        wsave(filepath, fig)
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

    savedir = joinpath(savedir_root, "mean_error")
    for (i, fig) in enumerate(figs)
        filepath = joinpath(savedir, "$(cfmt("%02d", i)).png")
        wsave(filepath, fig)
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

    savedir = joinpath(savedir_root, "mean_error_smoothed")
    for (i, fig) in enumerate(figs)
        filepath = joinpath(savedir, "$(cfmt("%02d", i)).png")
        wsave(filepath, fig)
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

    savedir = joinpath(savedir_root, "mean_abserror_smoothed")
    for (i, fig) in enumerate(figs)
        filepath = joinpath(savedir, "$(cfmt("%02d", i)).png")
        wsave(filepath, fig)
    end

    # Plot initial and final ensemble from spinup.
    fig = my_pairplot(ensembles[1].ensemble)
    savedir = joinpath(savedir_root, "state")
    filepath = joinpath(savedir, "initial.png")
    wsave(filepath, fig)

    fig = my_pairplot(ensembles[end].ensemble)
    savedir = joinpath(savedir_root, "state")
    filepath = joinpath(savedir, "final.png")
    wsave(filepath, fig)

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
        savedir = joinpath(savedir_root, "std")
        for (i, fig) in enumerate(figs)
            filepath = joinpath(savedir, "$(cfmt("%02d", i)).png")
            wsave(filepath, fig)
        end

        smoothed_stds = imfilter(metrics.vars_vec, ImageFiltering.Kernel.gaussian((0, N_states * 0.01))) .^ 0.5
        figs = plot_state_over_time(metrics.ts, smoothed_stds; plot_kwargs...)
        savedir = joinpath(savedir_root, "std_smoothed")
        for (i, fig) in enumerate(figs)
            filepath = joinpath(savedir, "$(cfmt("%02d", i)).png")
            wsave(filepath, fig)
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
        savedir = joinpath(savedir_root, "spread")
        for (i, fig) in enumerate(figs)
            filepath = joinpath(savedir, "$(cfmt("%02d", i)).png")
            wsave(filepath, fig)
        end

        smoothed_spread = imfilter(metrics.spread, ImageFiltering.Kernel.gaussian((N_states * 0.01,)))
        figs = plot_error_metric_over_time(metrics.ts, smoothed_spread; plot_kwargs...)
        savedir = joinpath(savedir_root, "spread_smoothed")
        for (i, fig) in enumerate(figs)
            filepath = joinpath(savedir, "$(cfmt("%02d", i)).png")
            wsave(filepath, fig)
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
        savedir = joinpath(savedir_root, "mean_rmse")
        for (i, fig) in enumerate(figs)
            filepath = joinpath(savedir, "$(cfmt("%02d", i)).png")
            wsave(filepath, fig)
        end

        smoothed_rmses = imfilter(metrics.rmses .^ 2, ImageFiltering.Kernel.gaussian((N_states * 0.01,))) .^ 0.5
        figs = plot_error_metric_over_time(metrics.ts, smoothed_rmses; plot_kwargs...)
        savedir = joinpath(savedir_root, "mean_rmse_smoothed2")
        for (i, fig) in enumerate(figs)
            filepath = joinpath(savedir, "$(cfmt("%02d", i)).png")
            wsave(filepath, fig)
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
        savedir = joinpath(savedir_root, "spread_mean_rmse")
        for (i, fig) in enumerate(figs)
            ax = fig.content[1]
            plot_disjoint_lines!(ax, metrics.ts, metrics.rmses; label="rmse", plot_kwargs_rmse...)
            fig[1, 2] = Legend(fig, ax; unique=true)
            filepath = joinpath(savedir, "$(cfmt("%02d", i)).png")
            wsave(filepath, fig)
        end

        # Plot smoothed values.
        smoothed_spread = imfilter(metrics.spread, ImageFiltering.Kernel.gaussian((N_states * 0.01,)))
        figs = plot_error_metric_over_time(metrics.ts, smoothed_spread; label="spread", plot_kwargs_full...)

        smoothed_rmses = imfilter(metrics.rmses .^ 2, ImageFiltering.Kernel.gaussian((N_states * 0.01,))) .^ 0.5
        savedir = joinpath(savedir_root, "spread_mean_rmse_smoothed")
        for (i, fig) in enumerate(figs)
            ax = fig.content[1]
            plot_disjoint_lines!(ax, metrics.ts, smoothed_rmses; label="rmse", plot_kwargs_rmse...)
            fig[1, 2] = Legend(fig, ax; unique=true)
            filepath = joinpath(savedir, "$(cfmt("%02d", i)).png")
            wsave(filepath, fig)
        end
    end
end
