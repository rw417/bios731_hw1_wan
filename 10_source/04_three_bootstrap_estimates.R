####################
# Calculate the Wald standard error and CI for all bootstrap runs
wald_estimates <- function(sim_beta_hat, bs_results, beta_true, alpha = 0.05) {
  # Estimate standard error of beta hat from bootstrap results
  se_b <- sd(bs_results[, "beta_hat_b"], na.rm=TRUE)
  
  # Estimate CI and coverage
  ci_lower <- sim_beta_hat + qt(alpha / 2, df = dim(bs_results)[1] - 1) * se_b
  ci_upper <- sim_beta_hat + qt(1 - alpha / 2, df = dim(bs_results)[1] - 1) * se_b
  coverage <- (beta_true > ci_lower & beta_true < ci_upper)

  bs_estimates <- c(sim_beta_hat, se_b, ci_lower, ci_upper, coverage)
  names(bs_estimates) <- c("sim_beta_hat", "se_b", "ci_lower", "ci_upper", "coverage")

  return(bs_estimates)
}

####################
# Calculate the Wald standard error and percentile-based CI for all bootstrap runs
percentile_estimates <- function(sim_beta_hat, bs_results, beta_true, alpha = 0.05) {
  # Estimate standard error of beta hat from bootstrap results
  se_b <- sd(bs_results[, "beta_hat_b"], na.rm=TRUE)
  
  # Estimate CI and coverage
  ci_lower <- quantile(bs_results[, "beta_hat_b"], alpha / 2, na.rm=TRUE)
  ci_upper <- quantile(bs_results[, "beta_hat_b"], 1 - alpha / 2, na.rm=TRUE)
  coverage <- (beta_true > ci_lower & beta_true < ci_upper)
  
  bs_estimates <- c(sim_beta_hat, se_b, ci_lower, ci_upper, coverage)
  names(bs_estimates) <- c("sim_beta_hat", "se_b", "ci_lower", "ci_upper", "coverage")

  return(bs_estimates)
}

####################
# Calculate the bootstrap t standard error and CI for all bootstrap runs
t_estimates <- function(sim_beta_hat, bs_results, beta_true, alpha = 0.05) {
  # Estimate standard error of beta hat from bootstrap results
  se_b <- sd(bs_results[, "beta_hat_b"], na.rm=TRUE)

  # Estimate CI and coverage
  t_quants <- quantile(bs_results[, "t_star"], probs = c(alpha / 2, 1 - (alpha / 2)), na.rm=TRUE)
  ci_lower <- sim_beta_hat - t_quants[2] * se_b
  ci_upper <- sim_beta_hat - t_quants[1] * se_b
  coverage <- (beta_true > ci_lower & beta_true < ci_upper)

  bs_estimates <- c(sim_beta_hat, se_b, ci_lower, ci_upper, coverage)
  names(bs_estimates) <- c("sim_beta_hat", "se_b", "ci_lower", "ci_upper", "coverage")

  return(bs_estimates)
}
