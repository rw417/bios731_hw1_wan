####################################################################
# Julia Wrobel
# January 2025
#
# This file produces simulations for linear regression under different data
# generating scenarios
####################################################################


library(tidyverse)


###############################################################
## define or source functions used in code below
###############################################################

source(here::here("source", "01_simulate_data.R"))
source(here::here("source", "02_apply_methods.R"))
source(here::here("source", "03_extract_estimates.R"))


###############################################################
## set simulation design elements
###############################################################

# how are you justifying nsim?
nsim = 900


n = c(10, 50, 100, 500)
beta_true = c(0.5, 2)
sigma2_true = c(1, 5)

params = expand.grid(n = n,
                     n_sim = nsim,
                     beta_true = beta_true,
                     sigma2_true = sigma2_true)



# define simulation scenario
scenario = 1
params = params[scenario,]




###############################################################
## start simulation code
###############################################################

# generate a random seed for each simulated dataset
seed = floor(runif(nsim, 1, 10000))
results = as.list(rep(NA, nsim))

for(i in 1:nsim){

  set.seed(seed[i])

  ####################
  # simulate data
  simdata = get_simdata(n = params$n,
                        beta_treat = params$beta_true,
                        sigma2 = params$sigma2_true)

  ####################
  # apply method(s)
  fit = fit_model(simdata)

  ####################
  # calculate estimates
  estimates = get_estimates(model_fit = fit,
                            true_beta = params$beta_true)

  ####################
  # store results, including estimates, speed, parameter scenarios
  estimates = estimates %>%
    mutate(true_beta = params$beta_true,
           n = params$n,
           sigma2_true = params$sigma2_true)

  results[[i]] = estimates

}

####################
# save results
# note that I am saving results outside of the for loop. For slow simulations,
  # you may want to save each iteration separately
filename = paste0("scenario_", scenario, ".RDA")
save(results,
     file = here::here("results", filename))

