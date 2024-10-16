# Define parameters.
DType = Dict{String,Any}

params_transition = DType(
    "sigma" => 10,
    "rho" => 28,
    "beta" => 8 / 3,
    "scaling" => 1,
    "ministep_nt" => missing,
    "ministep_dt" => 0.05,
)

params_exec = DType(
    "workers" => 0,
    "transitioner_distributed_type" => :none,
    "observer_distributed_type" => :none,
)

params = DType(
    "format" => "v0.2",
    "ground_truth" => DType(
        "format" => "v0.1",
        "transition" => params_transition,
        "observation" =>
            DType("noise_scale" => 2, "timestep_size" => 0.1, "num_timesteps" => 600),
    ),
    "ensemble" => DType(
        "size" => 100,
        "seed" => 9347215,
        "prior" => "gaussian",
        "prior_params" => [0.0, 1.0],
        "spinup" => DType(
            "transition" => params_transition,
            "exec" => params_exec,
            "assimilation_type" => "sequential",
            "num_timesteps" => 200,
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
        "num_timesteps" => 400,
        "transition_noise_scale" => 0.0,
        "exec" => params_exec,

        ## EnKF params
        "algorithm" => "enkf",
        "include_noise_in_y_covariance" => true,
        "multiplicative_prior_inflation" => 0.1,
        "observation_noise_stddev" => 2.0,
        "observation_noise_type" => "diagonal",

        # ## NF params
        # "algorithm" => "nf",
        # "glow" => DType(
        #     "L" => 3,
        #     "K" => 9,
        #     "n_hidden" => 8,
        #     "split_scales" => false,
        # ),
        # "training" => DType(
        #     "n_epochs" => 32,
        #     "batch_size" => 50,
        #     "noise_lev_x" => 0.005f0,
        #     "noise_lev_y" => 0.0f0,
        #     "num_post_samples" => 50,
        #     "validation_perc" => 0.5,
        #     "n_condmean" => 0,
        # ),
        # "optimizer" => DType(
        #     "lr" => 1.0f-3,
        #     "clipnorm_val" => 3.0f0,
        # ),
    ),
);
