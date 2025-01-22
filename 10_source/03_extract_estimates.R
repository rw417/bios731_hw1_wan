library(broom)

get_estimates = function(model_fit, true_beta){

  tidy(model_fit, conf.int = TRUE) %>%
    filter(term == "x") %>%
    mutate(coverage = ifelse(true_beta >= conf.low & true_beta <= conf.high, 1, 0)) %>%
    rename(beta_hat = estimate) %>%
    select(term, beta_hat, std.error, coverage)
}

