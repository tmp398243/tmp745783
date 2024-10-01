params = Dict(
    "format" => "v0.1",
    "transition" => Dict(
        "sigma" => 10,
        "rho" => 28,
        "beta" => 8 / 3,
        "scaling" => 1,
        "ministep_nt" => missing,
        "ministep_dt" => 0.05,
    ),
    "observation" =>
        Dict("noise_scale" => 2, "timestep_size" => 0.1, "num_timesteps" => 10),
    "ensemble" => Dict(
        "size" => 10,
        "seed" => 9347215,
        "prior" => "gaussian",
        "prior_params" => [0.0, 1.0],
    ),
    "spinup" => Dict(
        "num_timesteps" => 10,
        "transition_noise_scale" => 1.0,

        ## EnKF params
        "algorithm" => "enkf",
        "include_noise_in_y_covariance" => true,
        "multiplicative_prior_inflation" => 0.0,
        "observation_noise_stddev" => 2.0,
        "observation_noise_type" => "diagonal",
    ),
)
