

fit_model = function(simulated_data){
  lm(y ~ x, data = simulated_data)
}
