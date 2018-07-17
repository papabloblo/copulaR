ajuste_var_cop <- function(datos_train,
                           datos_valid,
                           datos_test,
                           var_iter,
                           num_sim,
                           max_bins,
                           bin_target,
                           num_obs_fit){
  
  errores <- list()
  
  train2 <- unique(train[, var_iter, drop = FALSE])
  
  # ¿duplicated más rápido que unique? -> COMPROBAR
  # train2 <- data.frame(train2[!duplicated(train2),])
  colnames(train2) <- var_iter
  
  # ¿Por qué if solo si NO HAY NINGUNA VARIABLE UNARIA?
  if (sum(apply(train2,
                2,
                function(x){length(unique(x))}
                ) == 1
          ) == 0){
    
    valid2 <- unique(valid[, var_iter, drop = FALSE])
    test2 <- unique(test[, var_iter, drop = FALSE])
    
    
    train2 <- unique(rbind(train2, valid2, test2))
    
    
    
    train_var <- unique(cbind(train[, var_iter, drop = FALSE], error = datos_train$error))
    variables <- colnames(train_var)
    
    ########################## ¿?
    if (!is.null(max_bins)){
      num_valores <- max(c(1, 
                           floor(max_bins^(1/(length(var_iter))))
                           )
                         )
    } else {
      num_valores <- 0
    }
    ########################## ¿?
    
    aprox_variables <- list()
    for (i in var_iter){
      if (!is.null(max_bins)){
        # ----POR AQUÍ----------
        aprox_variables[[i]] <- quantile(train2[, i], probs = seq(0, 1, length.out = num_valores))
        train2[, paste0(colnames(train2)[i], '_hist')] <- dplyr::ntile(train2[, i], num_valores)
        train_var[,i] <- apply(as.matrix(train_var[,i]),
                               1,
                               function(x){
                                 aprox_variables[[i]][which.min(abs(aprox_variables[[i]] - x))]
                               })
      } else {
        train2 <- cbind(train2, train2[, i])
        colnames(train2)[ncol(train2)] <- paste0(i, '_hist')
      }
    }
    
    if (!is.null(max_bins)){
      if (bin_target){
        aprox_variables[[length(variables)]] <- quantile(train_var[,length(variables)], probs = seq(0,1, length.out = num_valores))
        train_var[,length(variables)] <- apply(as.matrix(train_var[,length(variables)]),
                                               1,
                                               function(x){
                                                 aprox_variables[[length(variables)]][which.min(abs(aprox_variables[[length(variables)]] - x))]
                                               })
      }
    }
    
    ini <- Sys.time()
    
    mejor_copula <- copula.optima(train_var, num_obs_fit)
    
    if (!mejor_copula$indep){
      resultados <- data.frame()
      if (is.null(max_bins)){
        n <- nrow(train2)
      } else {
        n <- min(c(max_bins, nrow(train2))) 
      }
      if ((n*num_sim)>=1000000){
        num_iter <- floor((n*num_sim)/1000000) + 1
        fila_ini <- 1
        fila_fin <- min(c(floor(fila_ini + (1000000/num_sim)),n))
        for (i in 1:num_iter){
          train_aux <- train2[fila_ini:fila_fin,]
          resultados_aux <- puntuacion_copula_opt(datos_iter = data.frame(train_aux),
                                                  n.ventas = num_sim,
                                                  copulaoptima=mejor_copula$copulaoptima,
                                                  train =  train_var)
          resultados <- rbind(resultados, resultados_aux)
          fila_ini <- fila_fin + 1
          fila_fin <- min(c(floor(fila_ini + (1000000/num_sim)),n))
        }
      } else {
        resultados <- puntuacion_copula_opt(datos_iter = train2,
                                            n.ventas = num_sim,
                                            copulaoptima=mejor_copula$copulaoptima,
                                            train =  train_var)
      }
      
      resultados2 <- resultados
      names(resultados2)[names(resultados2) == 'ERROR'] <- 'ERROR_COP'
      resultados3 <- datos_train %>% left_join(resultados2, by = var_iter)
      resultados3.5 <- datos_valid %>% left_join(resultados2, by = var_iter)
      resultados4 <- datos_test %>% left_join(resultados2, by = var_iter)
      resultados3$pred_nueva <- resultados3$PREDICCION/(1-resultados3$ERROR_COP)
      resultados3$nuevo_error <- (resultados3$Target-resultados3$pred_nueva)/resultados3$Target
      resultados3.5$pred_nueva <- resultados3.5$PREDICCION/(1-resultados3.5$ERROR_COP)
      resultados3.5$nuevo_error <- (resultados3.5$Target-resultados3.5$pred_nueva)/resultados3.5$Target
      resultados4$pred_nueva <- resultados4$PREDICCION/(1-resultados4$ERROR_COP)
      resultados4$nuevo_error <- (resultados4$Target-resultados4$pred_nueva)/resultados4$Target
      datos_train$PREDICCION <- resultados3$pred_nueva
      datos_train$ERROR <- resultados3$nuevo_error
      datos_valid$PREDICCION <- resultados3.5$pred_nueva
      datos_valid$ERROR <- resultados3.5$nuevo_error
      datos_test$PREDICCION <- resultados4$pred_nueva
      datos_test$ERROR <- resultados4$nuevo_error
      
      errores[[1]] <- datos_train
      errores[[2]] <- datos_valid
      errores[[3]] <- datos_test
      errores[[4]] <- data.frame(mejor_copula_var = as.character(mejor_copula$aic[1,1]),
                                 ind_indepCopula = 0)
      resultados2 <- resultados2 %>% left_join(train2, by = var_iter)
      errores[[5]] <- mejor_copula$copulaoptima
      info_iter <- resultados2[,c(paste0(var_iter, '_hist'), 'ERROR_COP')]
      info_iter <- info_iter[!duplicated(info_iter),]
      errores[[6]] <- info_iter
      errores[[7]] <- aprox_variables
      
    } else {
      errores[[1]] <- data.frame()
      errores[[2]] <- data.frame()
      errores[[3]] <- data.frame()
      errores[[4]] <- data.frame('',
                                 ind_indepCopula = 1
      )
      errores[[5]] <- data.frame()
      errores[[6]] <- data.frame()
      errores[[7]] <- list()
    }
  } else {
    errores[[1]] <- data.frame()
    errores[[2]] <- data.frame()
    errores[[3]] <- data.frame()
    errores[[4]] <- data.frame('',
                               ind_indepCopula = 1
    )
    errores[[5]] <- data.frame()
    errores[[6]] <- data.frame()
    errores[[7]] <- list()
  }
  
  return(errores)
}