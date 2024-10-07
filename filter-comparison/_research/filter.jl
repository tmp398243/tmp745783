using LinearAlgebra: Diagonal
using EnsembleKalmanFilters: EnKF

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
        filter = EnKF(R; params);
        return filter
    elseif alg == "nf"
        include("nf.jl")
        filter = NormalizingFlow(; params)
        return filter
    end
    throw(ArgumentError("Unknown filter algorithm: $alg"))
end
