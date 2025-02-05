# load libraries
library(stats)

get_simdata <- function(n, beta_treat, error_form, error_sigma2) {
  beta0 <- 1
  x <- rbinom(n, 1, prob = 0.5)
  
  # if all x's are the same, run again
  while (length(unique(x)) == 1) {
    x <- rbinom(n, 1, prob = 0.5)
  }
  
  if (error_form == "normal") {
    epsilon <- rnorm(n, 0, sd = sqrt(error_sigma2))
  } else if (error_form == "lognormal") {
    epsilon <- rlnorm(n, 0, sdlog = sqrt(error_sigma2))
  } else {
    stop("Error form must be 'normal' or 'lognormal'")
  }
  y <- beta0 + beta_treat * x + epsilon

  # Put x and y into a matrix
  simdata <- cbind(y, x)
  return(simdata)
}
