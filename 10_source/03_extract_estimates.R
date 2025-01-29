library(broom)

####################
# Get estimates of one iteration in bootstrap
get_estimates = function(model_fit){
  
  tidy(model_fit, conf.int = TRUE) %>%
    filter(term == "x") %>%
    rename(beta_hat = estimate) %>%
    select(beta_hat)
}

