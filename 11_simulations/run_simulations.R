####################################################################
# Robert Wan
# January 2025
#
# This file produces simulations for linear regression under different data
# generating scenarios
####################################################################

library(stats)
library(foreach)
library(doParallel)

###############################################################
## define or source functions used in code below
###############################################################

source(here::here("10_source", "01_simulate_data.R"))
source(here::here("10_source", "02_apply_methods.R"))
source(here::here("10_source", "03_extract_estimates.R"))
source(here::here("10_source", "04_three_bootstrap_estimates.R"))
source(here::here("10_source", "05_run_bootstrap.R"))

###############################################################
## set simulation design elements
###############################################################

# nsim calculated based on desired coverage of 95%
# with Monte Carlo error of no more than 1%
nsim <- 475
nboot <- 200
nboot_t <- 100

n <- c(10, 50, 100)
beta_true <- c(0, 0.5, 2)
error_form <- c("normal", "lognormal")
error_sigma2 <- c(2)

params_grid <- expand.grid(
  n = n,
  nsim = nsim,
  beta_true = beta_true,
  error_form = error_form,
  error_sigma2 = error_sigma2
)

# Setup parallel backend
ncores <- detectCores() - 2 # Use available cores minus one
cl <- makeCluster(ncores, port=12000)
registerDoParallel(cl)

###############################################################
# start simulation code Parallelized Version ####
###############################################################

# Loop through all simulations
for (scenario in 1:nrow(params_grid)) {
  params <- params_grid[scenario, ]

  # Generate a random seed for each simulated dataset
  seed <- floor(runif(nsim, 1, 10000))

  print(paste("Running scenario: ", scenario))
  print(Sys.time())
  # Run simulation in parallel
  foreach_results <- foreach(
    i = 1:nsim, .combine = rbind, .multicombine = FALSE,
    .packages = c("here", "stats", "dplyr", "broom"),
    .export = c(
      "get_simdata", "fit_model", "get_estimates",
      "run_bootstrap", "wald_estimates",
      "percentile_estimates", "t_estimates"
    )
  ) %dopar% {
    # Set seed for reproducibility within each worker
    set.seed(seed[i])
    
    # if (i %% 50 == 0) {
    #   message(paste("Reached simulation: ", i))
    #   flush.console()
    # }

    ####################
    # Simulate data
    simdata <- get_simdata(
      n = params$n,
      beta_treat = params$beta_true,
      error_form = params$error_form,
      error_sigma2 = params$error_sigma2
    )

    ####################
    # Apply model on simulated data
    sim_fit <- fit_model(simdata)

    ####################
    # Calculate simulation estimates
    sim_beta_hat <- get_estimates(model_fit = sim_fit)

    #####################
    # Run bootstrap on simulated data using different methods
    loop_results <- run_bootstrap(simdata, sim_beta_hat, nboot, nboot_t)
    
    bs_results <- loop_results[[1]]
    bs_runtime <- loop_results[[2]]
    
    #####################
    # Calculate bootstrap estimates
    start_wald <- proc.time()[3]
    wald = wald_estimates(sim_beta_hat, bs_results, beta_true = params$beta_true, alpha = 0.05)
    end_wald <- proc.time()[3]
    bs_runtime['wald'] <- bs_runtime['wald'] + end_wald - start_wald
    
    start_percentile <- proc.time()[3]
    percentile = percentile_estimates(sim_beta_hat, bs_results, beta_true = params$beta_true, alpha = 0.05)
    end_percentile <- proc.time()[3]
    bs_runtime['percentile'] <- bs_runtime['percentile'] + end_percentile - start_percentile
    
    start_t <- proc.time()[3]
    t = t_estimates(sim_beta_hat, bs_results, beta_true = params$beta_true, alpha = 0.05)
    end_t <- proc.time()[3]
    bs_runtime['t'] <- bs_runtime['t'] + end_t - start_t

    # Store results in a structured format
    bs_runtime <- c(bs_runtime, 0)
    rbind(wald,percentile,t,bs_runtime)
  }

  # Store results into each of the matrices
  wald_results <- foreach_results[seq(1, 4*nsim, by=4),]
  percentile_results <- foreach_results[seq(2, 4*nsim, by=4),]
  t_results <- foreach_results[seq(3, 4*nsim, by=4),]
  bs_runtime_results <- foreach_results[seq(4, 4*nsim, by=4),]
  
  final_results <- list(
    wald_results <- foreach_results[seq(1, 4*nsim, by=4),],
    percentile_results <- foreach_results[seq(2, 4*nsim, by=4),],
    t_results <- foreach_results[seq(3, 4*nsim, by=4),],
    bs_runtime_results <- foreach_results[seq(4, 4*nsim, by=4),]
  )
  
  ####################
  # save results
  # note that I am saving results outside of the for loop. For slow simulations,
  # you may want to save each iteration separately
  filename <- paste0("20250204_scenario_", scenario, ".RDS")
  saveRDS(final_results,
    file = here::here("30_results", filename)
  )
  
  print(paste("Finished scenario: ", scenario))
  print(Sys.time())
}

stopCluster(cl) # Shut down cluster


