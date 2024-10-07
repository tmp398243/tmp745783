
using Pkg: Pkg

Pkg.add("DrWatson")

using DrWatson
# Pkg.activate()
# Pkg.activate(@__DIR__)
@assert basename(dirname(Base.active_project())) == "intense"

try
    using Lorenz63Filter
catch
    Pkg.develop(; path=joinpath(@__DIR__, "..", ".."))
end

try
    using Ensembles
catch
    Pkg.add(; url="https://github.com/tmp398243/tmp32487543")
    using Ensembles
end

try
    using Lorenz63: Lorenz63
catch
    Ensembles.install(:Lorenz63)
end

try
    using EnsembleKalmanFilters: EnsembleKalmanFilters
catch
    Ensembles.install(:EnsembleKalmanFilters)
end

try
    using NormalizingFlowFilters: NormalizingFlowFilters
catch
    Ensembles.install(:NormalizingFlowFilters)
end

Pkg.add(["LinearAlgebra", "Random", "CairoMakie", "Statistics", "PairPlots", "ImageFiltering"])

# @quickactivate "intense"
