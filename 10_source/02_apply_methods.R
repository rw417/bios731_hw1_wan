fit_model = function(simdata){
  # apply linear regression
  model_x = cbind(1, matrix(simdata[,'x']))
  model_y = matrix(simdata[,'y'])
  lm.fit(model_x, model_y)
}
