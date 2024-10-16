
if get(ENV, "lorenz63filter_force_install", "false") == "true" ||
    basename(dirname(Base.active_project())) in ["v1.10", "v1.9", "v1.8", "v1.7", "v1.6"]
    using Pkg: Pkg

    Pkg.activate(joinpath(@__DIR__, ".."))
    @assert basename(dirname(Base.active_project())) == "filter-comparison"

    try
        using Lorenz63Filter: Lorenz63Filter
    catch
        path = get(ENV, "lorenz63filter_path", joinpath(@__DIR__, "..", "..", ".."))
        Pkg.develop(; path)
    end

    try
        using Ensembles: Ensembles
    catch
        Pkg.add(; url="https://github.com/tmp398243/tmp32487543")
        using Ensembles: Ensembles
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

    Pkg.add([
        "DrWatson",
        "LinearAlgebra",
        "Random",
        "CairoMakie",
        "Statistics",
        "PairPlots",
        "ImageFiltering",
        "JLD2",
        "Format",
        "Configurations",
        "TerminalLoggers",
        "ProgressLogging",
        "Logging",
        "Markdown",
        "Distributed",
    ])

    Pkg.instantiate()
end
