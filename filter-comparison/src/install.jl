
if basename(dirname(Base.active_project())) in ["v1.10", "v1.9", "v1.8", "v1.7", "v1.6"]
    using Pkg: Pkg

    Pkg.add("DrWatson")

    using DrWatson: @quickactivate
    Pkg.activate(joinpath(@__DIR__, ".."))
    @assert basename(dirname(Base.active_project())) == "filter-comparison"

    try
        import Lorenz63Filter
    catch
        Pkg.develop(; path=joinpath(@__DIR__, "..", ".."))
    end

    try
        import Ensembles
    catch
        Pkg.add(; url="https://github.com/tmp398243/tmp32487543")
        import Ensembles
    end

    try
        import Lorenz63: Lorenz63
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

    Pkg.add(["LinearAlgebra", "Random", "CairoMakie", "Statistics", "PairPlots", "ImageFiltering", "JLD2", "Format", "Configurations", "TerminalLoggers", "ProgressLogging", "Logging"])

    Pkg.instantiate()
end
