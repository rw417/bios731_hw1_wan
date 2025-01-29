####################################################################
# Robert Wan
# January 2025
#
# This file produces simulations for linear regression under different data
# generating scenarios
####################################################################


library(tidyverse)
library(stats)

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
nsim <- 50
nboot <- 50
nboot_t <- 25

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
cl <- makeCluster(ncores)
registerDoParallel(cl)

###############################################################
# start simulation code Parallelized Version ####
###############################################################

# Loop through all simulations
for (scenario in 1:nrow(params_grid)) {
  params <- params_grid[scenario, ]

  # Generate a random seed for each simulated dataset
  seed <- floor(runif(nsim, 1, 10000))

  # Genereate a list of matrices to store results
  sim_results <- replicate(3, matrix(NA,
    nrow = nsim, ncol = 5,
    dimnames = list(
      NULL, c("sim_beta_hat", "se_b", "ci_lower", "ci_upper", "coverage")
    )
  ), simplify = FALSE)

  names(sim_results) <- c("wald", "percentile", "t")

  print(paste("Running scenario: ", scenario))
  print(Sys.time())
  # Run simulation in parallel
  sim_results <- foreach(
    i = 1:nsim, .combine = list, .multicombine = TRUE,
    .packages = c("here", "stats", "dplyr", "broom"),
    .export = c(
      "get_simdata", "fit_model", "get_estimates",
      "run_bootstrap", "wald_estimates",
      "percentile_estimates", "t_estimates"
    )
  ) %dopar% {
    # Source files inside each worker
    # source(here::here("10_source", "01_simulate_data.R"))
    # source(here::here("10_source", "02_apply_methods.R"))
    # source(here::here("10_source", "03_extract_estimates.R"))
    # source(here::here("10_source", "04_three_bootstrap_estimates.R"))
    # source(here::here("10_source", "05_run_bootstrap.R"))

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
    sim_beta_hat <- get_estimates(model_fit = sim_fit)$beta_hat

    #####################
    # Run bootstrap on simulated data using different methods
    bs_results <- run_bootstrap(simdata, sim_beta_hat, nboot, nboot_t)

    # Store results in a structured format
    list(
      wald = wald_estimates(sim_beta_hat, bs_results, beta_true = params$beta_true, alpha = 0.05),
      percentile = percentile_estimates(sim_beta_hat, bs_results, beta_true = params$beta_true, alpha = 0.05),
      t = t_estimates(sim_beta_hat, bs_results, beta_true = params$beta_true, alpha = 0.05)
    )
  }

  # Convert results into structured lists
  sim_results <- list(
    wald = do.call(rbind, lapply(sim_results, function(x) x$wald)),
    percentile = do.call(rbind, lapply(sim_results, function(x) x$percentile)),
    t = do.call(rbind, lapply(sim_results, function(x) x$t))
  )


  ####################
  # save results
  # note that I am saving results outside of the for loop. For slow simulations,
  # you may want to save each iteration separately
  filename <- paste0("scenario_", scenario, ".RDS")
  saveRDS(sim_results,
    file = here::here("30_results", filename)
  )
  
  print(paste("Finished scenario: ", scenario))
  print(Sys.time())
}

stopCluster(cl) # Shut down cluster


