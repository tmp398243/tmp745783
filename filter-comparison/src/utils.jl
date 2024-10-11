using Ensembles: AbstractEnsemble
using PairPlots: PairPlots
using Statistics: mean, var

function my_pairplot(ensemble::AbstractEnsemble)
    data = get_ensemble_matrix(ensemble)'
    my_pairplot(data)
end

function my_pairplot(data::AbstractArray)
    N = size(data, 1)
    colorrange = (0, size(data, 1))
    markersize = 5
    alpha = 1
    PairPlots.pairplot(
        data => (
            PairPlots.HexBin(;colormap=:Blues, colorrange),
            PairPlots.Scatter(;color=:red, markersize, alpha),
            PairPlots.MarginDensity(;bandwidth=0.1),
            PairPlots.MarginHist(),
        )
    )
end

function show_interactive(fig)
    if Base.is_interactive
        display(fig)
    else
        fig
    end
end


function get_ground_truth_iterator(ensembles_ts, observation_times)
    gt_index = 1
    gt_indices = Int64[]
    post_assim_indices = Int64[]
    for (i, t) in enumerate(ensembles_ts)
        while gt_index <= length(observation_times) && observation_times[gt_index] < t
            gt_index += 1
        end
        if gt_index > length(observation_times)
            error("Comparing at time $(t) is impossible because final ground-truth observation is at time $(observation_times[end])")
        end
        if observation_times[gt_index] != t
            error("No observation at time $(t). Closest are $(observation_times[gt_index-1]) and $(observation_times[gt_index])")
        end
        push!(gt_indices, gt_index)
        if i == length(ensembles_ts) || t < ensembles_ts[i+1]
            push!(post_assim_indices, i)
        end
    end
    return gt_indices, post_assim_indices
end

function rmse(ensemble, y_true)
    return sqrt(mean((ensemble .- y_true) .^ 2))
end

function compute_errors(gt_indices, ensembles_means_vec, ground_truth_states_vec)
    rmses = Vector{Float64}(undef, length(gt_indices))
    for (i, gt_index) in enumerate(gt_indices)
        rmses[i] = rmse(ensembles_means_vec[:, i], ground_truth_states_vec[:, gt_index])
    end
    return rmses
end

function compute_metrics(ensembles; ts_gt=nothing, ground_truth_states_vec=nothing)
    state_keys = ensembles[1].ensemble.state_keys
    means_vec = get_ensemble_matrix(state_keys, mean(e.ensemble) for e in ensembles)
    vars_vec = get_ensemble_matrix(state_keys, var(e.ensemble) for e in ensembles)
    ts = [e.t for e in ensembles]

    if isnothing(ts_gt)
        return (;
        ts,
        vars_vec,
        means_vec,
    )
    end
    gt_indices, post_assim_indices = get_ground_truth_iterator(ts, ts_gt)
    rmses = compute_errors(gt_indices, means_vec, ground_truth_states_vec)

    spread = sqrt.(mean(vars_vec; dims=1)[1, :])
    return (;
        ts,
        vars_vec,
        means_vec,
        gt_indices,
        post_assim_indices,
        rmses,
        spread,
        post_assim_rmses = (rmses[i] for i in post_assim_indices),
        post_assim_spread = (spread[i] for i in post_assim_indices),
    )
end