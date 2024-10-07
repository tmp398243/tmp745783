using LinearAlgebra: Diagonal

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
