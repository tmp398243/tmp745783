function update_header(content, pth; build_notebooks, build_scripts)
    links = []
    if build_notebooks
        push!(links, "[Jupyter notebook](main.ipynb)")
    end
    if build_scripts
        push!(links, "[plain script](main.jl)")
    end
    if length(links) == 0
        return content
    end
    project_link = "[Project.toml](Project.toml)"
    return """
        # # Reproducing example
        # The packages for this example are documented in the $project_link.
        # # Accessing example
        # This can also be accessed as a $(join(links, ", a", ", or a ")).
    """ * content
end
