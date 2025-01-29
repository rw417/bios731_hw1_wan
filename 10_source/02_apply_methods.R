fit_model = function(simdata){
  # apply linear regression
  lm(y ~ x, data = simdata)
}
