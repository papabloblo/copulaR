genera_mejor_iter <- function(pasos_stepwise, pred_train, pred_valid, pred_test, modelo){
  errores_train <- data.frame()
  errores_valid <- data.frame()
  errores_test <- data.frame()
  mejor_num_iteraciones <- pasos_stepwise[which.min(pasos_stepwise$error_valid), 'iteracion']
  for (j in 1:mejor_num_iteraciones){
    if (j == 1){
      errores_train <- rbind(errores_train,
                             data.frame(iter = j,
                                        error = pasos_stepwise[1, 'error_train'],
                                        variable = pasos_stepwise[1,'variable']
                             )
      )
      errores_valid <- rbind(errores_valid,
                             data.frame(iter = j,
                                        error = pasos_stepwise[1, 'error_valid'],
                                        variable = pasos_stepwise[1,'variable']
                             )
      )
      errores_test <- rbind(errores_test,
                            data.frame(iter = j,
                                       error = pasos_stepwise[1, 'error_test'],
                                       variable = pasos_stepwise[1,'variable']
                            )
      )
      
    } else if (j < mejor_num_iteraciones) {
      aux <- pasos_stepwise[pasos_stepwise$iteracion == j ,]
      errores_train <- rbind(errores_train,
                             data.frame(iter = j,
                                        error = aux[nrow(aux), 'error_train'],
                                        variable = aux[nrow(aux),'variable']
                             )
      )
      errores_valid <- rbind(errores_valid,
                             data.frame(iter = j,
                                        error = aux[nrow(aux), 'error_valid'],
                                        variable = aux[nrow(aux),'variable']
                             )
      )
      errores_test <- rbind(errores_test,
                            data.frame(iter = j,
                                       error = aux[nrow(aux), 'error_test'],
                                       variable = aux[nrow(aux) ,'variable']
                            )
      )
      
    } else {
      aux <- pasos_stepwise[pasos_stepwise$iteracion == j ,]
      errores_train <- rbind(errores_train,
                             data.frame(iter = j,
                                        error = aux[1, 'error_train'],
                                        variable = aux[1,'variable']
                             )
      )
      errores_valid <- rbind(errores_valid,
                             data.frame(iter = j,
                                        error = aux[1, 'error_valid'],
                                        variable = aux[1,'variable']
                             )
      )
      errores_test <- rbind(errores_test,
                            data.frame(iter = j,
                                       error = aux[1, 'error_test'],
                                       variable = aux[1 ,'variable']
                            )
      )
    }
  }
  
  pasos_stepwise_aux <- pasos_stepwise[pasos_stepwise$iteracion < mejor_num_iteraciones,]
  aux <- pasos_stepwise[pasos_stepwise$iteracion ==  mejor_num_iteraciones,]
  pasos_stepwise <- rbind(pasos_stepwise_aux,
                          aux[1,])
  
  modelo_final <- modelo
  if (mejor_num_iteraciones > 1){
    modelo_final[['iteraciones']] <- lapply(modelo[['iteraciones']][1:(mejor_num_iteraciones - 1)], function(x){x$final})
    modelo_final[['iteraciones']][[mejor_num_iteraciones]] <- modelo[['iteraciones']][[mejor_num_iteraciones]][['original']]
  } else {
    modelo_final[['iteraciones']] <- modelo_final[['iteraciones']][[1]][['original']]
  }
 
  return(list(errores_train = errores_train,
              errores_valid = errores_valid,
              errores_test = errores_test,
              pasos_stepwise = pasos_stepwise,
              pred_train = pred_train[[mejor_num_iteraciones]],
              pred_valid = pred_valid[[mejor_num_iteraciones]],
              pred_test = pred_test[[mejor_num_iteraciones]],
              modelo = modelo_final))
}