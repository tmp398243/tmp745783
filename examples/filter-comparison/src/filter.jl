using LinearAlgebra: Diagonal
using EnsembleKalmanFilters: EnKF
using NormalizingFlowFilters:
    cpu,
    ConditionalGlowOptions,
    NetworkConditionalGlow,
    OptimizerOptions,
    create_optimizer,
    TrainingOptions,
    NormalizingFlowFilter
using Configurations: from_dict

function Ensembles.assimilate_data(
    filter::NormalizingFlowFilter,
    ensemble,
    ensemble_obs_clean,
    ensemble_obs_noisy,
    y_obs,
    log_data,
)
    X_matrix = NormalizingFlowFilters.assimilate_data(
        filter,
        Float64.(get_ensemble_matrix(ensemble)),
        Float64.(get_ensemble_matrix(ensemble_obs_noisy)),
        get_member_vector(ensemble_obs_clean, y_obs),
        log_data,
    )
    members = get_ensemble_dicts(ensemble, X_matrix)
    posterior = Ensemble(members, ensemble.state_keys)
    return posterior
end

function get_filter(params)
    alg = params["algorithm"]
    if alg == "enkf"
        obs_std = params["observation_noise_stddev"]
        noise_type = params["observation_noise_type"]
        n = params["assimilation_type"] == "sequential" ? 1 : 3
        if noise_type == "diagonal"
            R = Diagonal(fill(Float64(obs_std)^2, n))
        else
            throw(ArgumentError("Unknown observation noise type: $noise_type"))
        end
        filter = EnKF(R; params)
        return filter
    elseif alg == "nf"
        glow_config = from_dict(ConditionalGlowOptions, params["glow"]; chan_x=3, chan_y=3)
        network = NetworkConditionalGlow(2, glow_config)

        optimizer_config = from_dict(OptimizerOptions, params["optimizer"])
        optimizer = create_optimizer(optimizer_config)

        training_config = from_dict(TrainingOptions, params["training"])

        device = cpu
        filter = NormalizingFlowFilter(network, optimizer; device, training_config)
        return filter
    end
    throw(ArgumentError("Unknown filter algorithm: $alg"))
end
