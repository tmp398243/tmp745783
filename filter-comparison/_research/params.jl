refresh_cache = false

function flatten_dict(d, prefix_delim = ".")
    new_d = empty(d, Any)
    for (key, value) in pairs(d)
        if isa(value, Dict)
             flattened_value = flatten_dict(value, prefix_delim)
             for (ikey, ivalue) in pairs(flattened_value)
                 new_d["$key.$ikey"] = ivalue
             end
        else
            new_d[key] = value
        end
    end
    return new_d
end

function expand_dict(d::Dict{String, <:Any}, prefix_delim::String = ".")
    result = Dict{String, Any}()
    for (key, value) in pairs(d)
        parts = split(key, prefix_delim)
        current = result
        for (i, part) in enumerate(parts)
            if i == length(parts)
                current[part] = value
            else
                if !haskey(current, part) || !isa(current[part], Dict)
                    current[part] = Dict{String, Any}()
                end
                current = current[part]
            end
        end
    end
    return result
end

get_params_iterator(params::Dict{String, <:Any}) = get_params_iterator(params, get(params, "params_list", nothing))
get_params_iterator(params::Dict{String, <:Any}, params_list::Nothing) = [params]
get_params_iterator(params::Dict{String, <:Any}, params_list::Dict{String, <:Any}) = get_params_iterator(params, [params_list])
function get_params_iterator(params::Dict{String, Any}, params_list::Vector{<:Dict{String, <:Any}})
    params_base = deepcopy(flatten_dict(params))
    delete!(params_base, "params_list")

    params_combs = []
    for p_mod in params_list
        p_mod = flatten_dict(p_mod)
        mod_keys = keys(p_mod)
        mod_values = [p_mod[k] for k in mod_keys]
        combinations = Iterators.product(values(mod_values)...)
        for combo in combinations
            new_params = deepcopy(params_base)
            for (key, value) in zip(mod_keys, combo)
                new_params[key] = value
            end
            push!(params_combs, new_params)
        end
    end
    return params_combs
end

# Define parameters.
DType = Dict{String, Any}
params_all = DType(
    "params_list" => [
        DType(
            "ensemble" => DType(
                "size" => [10, 60, 100, 200, 400, 600]
            ),
            "observation" => DType(
                "timestep_size" => [0.1]
            ),
        ),
        DType(
            "ensemble" => DType(
                "size" => [1000]
            ),
            "observation" => DType(
                "timestep_size" => [0.1, 0.2, 0.3, 0.4, 0.5]
            ),
        ),
        ## DType(
        ##     "estimator" => DType(
        ##         "multiplicative_prior_inflation" => [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
        ##     ),
        ## ),
    ],
    "format" => "v0.1",
    "transition" => DType(
        "sigma" => 10,
        "rho" => 28,
        "beta" => 8/3,
        "scaling" => 1,
        "ministep_nt" => missing,
        "ministep_dt" => 0.05,
    ),
    "observation" => DType(
        "noise_scale" => 2,
        "timestep_size" => missing,
        "num_timesteps" => 6000,
    ),
    "ensemble" => DType(
        "size" => missing,
        "seed" => 9347215,
        "prior" => "gaussian",
        "prior_params" => [0.0, 1.0],
    ),
    "spinup" => DType(
        "assimilation_type" => "sequential",
        "num_timesteps" => 2000,
        "transition_noise_scale" => 1.0,

        ## EnKF params
        "algorithm" => "enkf",
        "include_noise_in_y_covariance" => true,
        "multiplicative_prior_inflation" => 0.0,
        "observation_noise_stddev" => 2.0,
        "observation_noise_type" => "diagonal",
    ),
    "estimator" => DType(
        "assimilation_type" => "monolithic",
        "num_timesteps" => 4000,
        "transition_noise_scale" => 0.0,

        ## EnKF params
        "algorithm" => "enkf",
        "include_noise_in_y_covariance" => true,
        "multiplicative_prior_inflation" => 0.1,
        "observation_noise_stddev" => 2.0,
        "observation_noise_type" => "diagonal",
    )
);

params_iterator = get_params_iterator(params_all)
params_compat = params_iterator[1]
params = expand_dict(params_compat)

params_transition = DType(
        "sigma" => 10,
        "rho" => 28,
        "beta" => 8/3,
        "scaling" => 1,
        "ministep_nt" => missing,
        "ministep_dt" => 0.05,
    )
params = DType(
    "format" => "v0.2",
    "ground_truth" => DType(
        "format" => "v0.1",
        "transition" => params_transition,
        "observation" => DType(
            "noise_scale" => 2,
            "timestep_size" => 0.1,
            "num_timesteps" => 6000,
        ),
    ),
    "ensemble" => DType(
        "size" => 10,
        "seed" => 9347215,
        "prior" => "gaussian",
        "prior_params" => [0.0, 1.0],
        "spinup" => DType(
            "transition" => params_transition,
            "assimilation_type" => "sequential",
            "num_timesteps" => 2000,
            "transition_noise_scale" => 1.0,

            ## EnKF params
            "algorithm" => "enkf",
            "include_noise_in_y_covariance" => true,
            "multiplicative_prior_inflation" => 0.0,
            "observation_noise_stddev" => 2.0,
            "observation_noise_type" => "diagonal",
        ),
    ),
    "estimator" => DType(
        "assimilation_type" => "monolithic",
        "num_timesteps" => 4000,
        "transition_noise_scale" => 0.0,

        ## EnKF params
        "algorithm" => "enkf",
        "include_noise_in_y_covariance" => true,
        "multiplicative_prior_inflation" => 0.1,
        "observation_noise_stddev" => 2.0,
        "observation_noise_type" => "diagonal",
    )
);
