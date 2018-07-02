eval_metric_functions <- list()
eval_metric_functions[['MAPE']] <- function(Target, Prediccion, error){
  return(mean(100*error))
}
eval_metric_functions[['MEDAPE']] <- function(Target, Prediccion, error){
  return(median(100*error))
}
eval_metric_functions[['MSE']] <- function(Target, Prediccion, error){
  return(mean((Target-Prediccion)^2))
}
eval_metric_functions[['RMSE']] <- function(Target, Prediccion, error){
  return(sqrt(mean((Target-Prediccion)^2)))
}
eval_metric_functions[['MAE']] <- function(Target, Prediccion, error){
  return(mean(abs(Target-Prediccion)))
}
eval_metric_functions[['SMAPE']] <- function(Target, Prediccion, error){
  return(100*mean((abs(Prediccion-Target)/((abs(Prediccion)+abs(Target))/2))))
}










