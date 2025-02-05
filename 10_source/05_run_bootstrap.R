library(foreach)
library(doParallel)


####################
# Run bootstrap on simulated data
run_bootstrap <- function(simdata, sim_beta_hat, nboot, nboot_t, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Matrix to store results
  bs_results <- matrix(
    NA,
    nrow = nboot, ncol = 2,
    dimnames = list(NULL, c("beta_hat_b", "t_star"))
  )

  n_bs <- nrow(simdata)
  
  ####################
  # Start of Bootstrap
  start_of_bs <- proc.time()[3] # for timing
  t_loop_duration <- 0 # for timing
  for (b in 1:nboot) {
    ####################
    # sample with replacement from simdata
    bs_idx <- sample(1:n_bs, size = n_bs, replace = TRUE)
    bs_sample <- simdata[bs_idx, ]  # Direct indexing

    ####################
    # apply model to bs_sample
    bs_fit <- fit_model(bs_sample)

    ####################
    # calculate bootstrap estimates of beta
    bs_results[b, "beta_hat_b"] <- get_estimates(model_fit = bs_fit)

    ####################
    # Bootstrap with t estimates
    start_of_t_loop <- proc.time()[3]
    # Step 1: loop to run bootstrap with t estimates
    boot_mean_b <- numeric(nboot_t)
    for (k in 1:nboot_t) {
      ####################
      # sample with replacement from bs_sample
      bs_idx_k <- sample(1:n_bs, size = n_bs, replace = TRUE)
      bs_sample_k <- bs_sample[bs_idx_k, ]  # Direct indexing

      ####################
      # apply model on bs_sample_k
      bs_fit_k <- fit_model(bs_sample_k)

      ####################
      # calculate inner bootstrap estimates of beta
      boot_mean_b[k] <- get_estimates(model_fit = bs_fit_k)
    }

    # Step 2: calculate t_star
    se_star <- sd(boot_mean_b, na.rm = TRUE)
    bs_results[b, "t_star"] <- (bs_results[b, "beta_hat_b"] - sim_beta_hat) / se_star
    
    end_of_t_loop <- proc.time()[3]
    t_loop_duration <- t_loop_duration + end_of_t_loop - start_of_t_loop
  }
  end_of_bs <- proc.time()[3]
  bs_duration_total <- end_of_bs - start_of_bs
  bs_duration_no_t <- bs_duration_total - t_loop_duration
  
  bs_runtime <- c(bs_duration_no_t, bs_duration_no_t, bs_duration_total, t_loop_duration)
  names(bs_runtime) <- c("wald", "percentile", "t", "t_loop")
  
  # End of Bootstrap
  
  return(list(bs_results, bs_runtime))

  ####################
  # compute final bootstrap estimates
  # wald_results <- wald_estimates(
  #   sim_beta_hat, bs_results,
  #   beta_true = params$beta_true, alpha = 0.05
  # )
  # 
  # percentile_results <- percentile_estimates(
  #   sim_beta_hat, bs_results,
  #   beta_true = params$beta_true, alpha = 0.05
  # )
  # 
  # t_results <- t_estimates(
  #   sim_beta_hat, bs_results,
  #   beta_true = params$beta_true, alpha = 0.05
  # )
  # 
  # return(list(wald_results, percentile_results, t_results))
}


########################
# Parallelized Code ####
run_bootstrap_parallel <- function(simdata, sim_beta_hat, nboot, nboot_t,seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Set up parallel backend
  ncores <- detectCores() - 1 # Use all but one core
  cl <- makeCluster(ncores)
  registerDoParallel(cl)

  # Run outer bootstrap loop in parallel
  bs_results <- foreach(
    b = 1:nboot, .combine = rbind,
    .packages = c("stats", "dplyr", "broom", "foreach", "doParallel"), .export = c("get_estimates", "fit_model")
  ) %dopar% {
    # sample with replacement from simdata
    bs_sample <- simdata[
      sample(1:nrow(simdata),
        size = nrow(simdata),
        replace = TRUE
      ),
    ]

    # Fit model
    bs_fit <- fit_model(bs_sample)

    # Bootstrap estimate of beta
    beta_hat_b <- get_estimates(model_fit = bs_fit)$beta_hat

    t_star <- NA

    # Step 1: loop to run bootstrap with t estimates
    boot_mean_b <- numeric(nboot_t)
    for (k in 1:nboot_t) {
      ####################
      # sample with replacement from bs_sample
      bs_sample_k <- bs_sample[
        sample(1:nrow(simdata),
          size = nrow(simdata),
          replace = TRUE
        ),
      ]

      ####################
      # apply model on bs_sample_k
      bs_fit_k <- fit_model(bs_sample_k)

      ####################
      # calculate inner bootstrap estimates of beta
      boot_mean_b[k] <- get_estimates(
        model_fit = bs_fit_k
      )$beta_hat
    }

    # Compute t_star
    se_star <- sd(boot_mean_b, na.rm=TRUE)
    t_star <- (beta_hat_b - sim_beta_hat) / se_star

    c(beta_hat_b, t_star)
  }

  stopCluster(cl) # Shut down parallel cluster

  colnames(bs_results) <- c("beta_hat_b", "t_star")

  # End of Bootstrap

  ####################
  # compute final bootstrap estimates
  wald_results <- wald_estimates(
    sim_beta_hat, bs_results,
    beta_true = params$beta_true, alpha = 0.05
  )
  
  percentile_results <- percentile_estimates(
    sim_beta_hat, bs_results,
    beta_true = params$beta_true, alpha = 0.05
  )
  
  t_results <- t_estimates(
    sim_beta_hat, bs_results,
    beta_true = params$beta_true, alpha = 0.05
  )
  
  return(list(wald_results, percentile_results, t_results))

  # if (method == "wald") {
  #   bs_estimates = wald_estimates(
  #     sim_beta_hat, bs_results,
  #     beta_true = params$beta_true, alpha = 0.05
  #   )
  # } else if (method == "percentile") {
  #   bs_estimates = percentile_estimates(
  #     sim_beta_hat, bs_results,
  #     beta_true = params$beta_true, alpha = 0.05
  #   )
  # } else if (method == "t") {
  #   bs_estimates = t_estimates(
  #     sim_beta_hat, bs_results,
  #     beta_true = params$beta_true, alpha = 0.05
  #   )
  # }
}
