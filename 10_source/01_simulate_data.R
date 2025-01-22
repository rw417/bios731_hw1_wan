# load libraries
library(tidyverse)

#
#
#
get_simdata = function(n, beta_treat, sigma2){
  beta0 = 1
  x = rbinom(n, 1, prob = 0.5)
  epsilon = rnorm(n, 0, sd = sqrt(sigma2))
  y = beta0 + beta_treat * x + epsilon

  tibble(
    x = x,
    y = y
  )

}



