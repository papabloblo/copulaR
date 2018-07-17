library(ks)
library(magrittr)
library(dplyr)
library(data.table)
library(zoo)
library(writexl)
library(BBmisc)
library(caTools)
library(VineCopula)
library(rapportools)


# DEPENDECIAS -------------------------------------------------------------

source('genera_combinaciones_variables.R')
source('copula_optima_BI.R')
source('puntuacion_copulas_comb.R')
source('genera_mejor_iter.R')
source('ajuste_var_cop.R')
source('eval_metric_functions.R')



#' Title
#'
#' @param train data.frame. 
#' @param target character
#' @param valid 
#' @param test 
#' @param num_iter 
#' @param early_stopping_round 
#' @param num_sim 
#' @param max_bins 
#' @param num_obs_fit 
#' @param bin_target 
#' @param eval_metric 
#' @param verbosity 
#'
#' @return
#' @export
#'
#' @examples
copula.model <- function(train,
                         target,
                         valid = train,
                         test = train,
                         num_iter = 10,
                         early_stopping_round = 0,
                         num_sim = 500,
                         max_bins = NULL,
                         num_obs_fit = NULL,
                         bin_target = FALSE,
                         eval_metric  = "MAPE",
                         verbosity = TRUE){
  
  
  

  # ERRORS -----------------------------------------------------------------
  
  if (length(num_iter) > 1) {
    stop('num_iter must be of length 1')
  } else if (!is.numeric(num_iter) || num_iter <= 0) {
    stop('num_iter must be a positive number')
  } else {
    num_iter <- floor(num_iter)
  }
  
  if (length(early_stopping_round) > 1) {
    stop('early_stopping_round must be of length 1')
  } else if (!is.numeric(early_stopping_round) || early_stopping_round <= 0) {
    stop('early_stopping_round must be a positive number')
  } else {
    early_stopping_round <- floor(early_stopping_round)
  }
  
  if (length(num_sim) > 1) {
    stop('num_sim must be of length 1')
  } else if (!is.numeric(num_sim) || num_sim <= 0) {
    stop('num_sim must be a positive number')
  } else {
    num_sim <- floor(num_sim)
  }
  
  if (length(max_bins) > 1) {
    stop('max_bins must be of length 1')
  } else if (!is.numeric(max_bins) || max_bins <= 0) {
    stop('max_bins must be a positive number')
  } else {
    max_bins <- floor(max_bins)
  }
  
  if (length(num_obs_fit) > 1) {
    stop('num_obs_fit must be of length 1')
  } else if (!is.numeric(num_obs_fit) || num_obs_fit <= 0) {
    stop('num_obs_fit must be a positive number')
  } else {
    num_obs_fit <- floor(num_obs_fit)
  }
  
  
  if (!is.logical(bin_target)){
    stop("bin_taget must be TRUE or FALSE")
  }
  
  if (!is.logical(verbosity)){
    stop("verbosity must be TRUE or FALSE")
  }
  
  if ( !(target %in% names(train))){
    stop("target must be in train")
  }
 
  
  if (!eval_metric %in% c("MAPE",
                          "MEDAPE",
                          "MSE",
                          "RMSE",
                          "MAE",
                          "SMAPE")){
    stop(paste("eval_metric", eval_metric, "not soported"))
  }

  if (ncol(test) != ncol(train) || all(sort(names(train)) == sort(names(test)))){
    stop("test and train don't have the same columns")
  } 
  
  if (ncol(valid) != ncol(train) || all(sort(names(valid)) == sort(names(test)))){
    stop("valid and train don't have the same columns")
  } 
  

  if (!is.matrix(train)) train <- as.matrix(train)
  if (!is.matrix(test)) train <- as.matrix(test)
  if (!is.matrix(valid)) train <- as.matrix(valid)
  
  colnames(train)[colnames(train) == target] <- "target"
  colnames(valid)[colnames(valid) == target] <- "target"
  if (target %in% colnames(test)){
    colnames(test)[colnames(test) == target] <- "target"
  } else {
    test <- cbind(test, NA)
    colnames(test)[dim(test)[2]] <- "target"
  }
  
  variables <- setdiff(colnames(train), "target")
  num_variables <- length(variables)
  
  ########
  max_dim_copulas <- 2 
  #######
  
  errors <- data_frame(iter = 0L,
                        error = 0,
                        var = character(1))
  
  errors_train <- errors
  errors_valid <- errors
  errors_test <- errors
  
  stepwise <- data_frame()
  
  iteracion <- 1
  
  combinaciones_variables <- combn(variables, max_dim_copulas - 1, simplify = FALSE)
  
  variables_cruce <- unlist(lapply(combinaciones_variables, 
                                   function(x) paste(x, collapse = ", ")
                                   )
                            )

  
  pred_train <- list()
  pred_valid <- list()
  pred_test <- list()
  
  modelo <- list(
    train = train,
    max_bins = max_bins,
    num_sim = num_sim,
    bin_target = bin_target,
    iteraciones = list()
  )
  
  i <- 1
  
  while (i <= num_iter){
    
    # comprobar!!!!
    if (i > 1){
      if (early_stopping_round > 0){
        
        if ((i - errores_valid[which.min(errores_valid$error), 'iter'] - 2) == 
            early_stopping_round){
          
          pasos_stepwise <- pasos_stepwise[pasos_stepwise$paso<(iteracion - 1),]
          tablas_output <- genera_mejor_iter(pasos_stepwise, pred_train, pred_valid, pred_test, modelo)
          
          if ((verbosity) & (nrow(pasos_stepwise[pasos_stepwise$paso==iteracion,])>0)){
            print(pasos_stepwise[pasos_stepwise$paso==iteracion,])
          }
          i <- num_iter + 1
          next
        }
      }
    }
    
    
    assign(paste0('errores_train_var_', i), data.frame())
    assign(paste0('errores_valid_var_', i), data.frame())
    assign(paste0('errores_test_var_', i), data.frame())
    
    for (j in 1:length(combinaciones_variables)){
      
      if (i == 1){
        
        pred <- mean(train[, "target"])
        
        datos_train <- list(target = train[, "target"], prediction = pred)
        datos_train$error <- (datos_train$target - datos_train$pred)/datos_train$target
        
        datos_valid <- list(target = valid[, "target"], prediction = pred)
        datos_valid$error <- (datos_valid$target - datos_valid$pred)/datos_valid$target
        
        datos_test  <- list(target = test[,  "target"], prediction = pred)
        datos_test$error <- (datos_test$target - datos_test$pred)/datos_test$target
        
        ###########################
        errores_train$error <- round(eval_metric_functions[[eval_metric]]
                                     (datos_train$Target,
                                       datos_train$PREDICCION,
                                       datos_train$ERROR),5)
        
        errores_valid$error <- round(eval_metric_functions[[eval_metric]]
                                     (datos_valid$Target,
                                       datos_valid$PREDICCION,
                                       datos_valid$ERROR),5)
        
        errores_test$error <- round(eval_metric_functions[[eval_metric]]
                                     (datos_test$Target,
                                       datos_test$PREDICCION,
                                       datos_test$ERROR),5)
        ###############  
      } else  {
        datos_train <- datos_train_fija
        datos_valid <- datos_valid_fija
        datos_test <- datos_test_fija
      }
      
      datos_train
      datos_valid
      datos_test
      var_iter = combinaciones_variables[[j]]
      num_sim
      max_bins
      bin_target
      num_obs_fit
      
      # Convertir en lista !!!!
      assign(paste0('errores_', i, '_', j), 
             ajuste_var_cop(datos_train,
                            datos_valid,
                            datos_test,
                            combinaciones_variables[[j]],
                            num_sim,
                            max_bins,
                            bin_target,
                            num_obs_fit)
             )
      
      assign(paste0('errores_train_var_', i), 
             rbind(get(paste0('errores_train_var_', i)),
                   data.frame(var = paste0(variables[combinaciones_variables[[j]]], collapse = ','),
                              error = ifelse(get(paste0('errores_', i, '_', j))[[4]]$ind_indepCopula==0, 
                                             round(eval_metric_functions[[eval_metric]]
                                             (get(paste0('errores_', i, '_', j))[[1]]$Target,
                                               get(paste0('errores_', i, '_', j))[[1]]$PREDICCION,
                                               get(paste0('errores_', i, '_', j))[[1]]$ERROR),5),
                                             Inf)
                   )
             )
      )
      
      assign(paste0('errores_valid_var_', i), 
             rbind(get(paste0('errores_valid_var_', i)),
                   data.frame(var = paste0(variables[combinaciones_variables[[j]]], collapse = ','),
                              error = ifelse(get(paste0('errores_', i, '_', j))[[4]]$ind_indepCopula==0, 
                                             round(eval_metric_functions[[eval_metric]]
                                             (get(paste0('errores_', i, '_', j))[[2]]$Target,
                                               get(paste0('errores_', i, '_', j))[[2]]$PREDICCION,
                                               get(paste0('errores_', i, '_', j))[[2]]$ERROR),5),
                                             Inf)
                   )
             )
      )
      
      assign(paste0('errores_test_var_', i),
             rbind(get(paste0('errores_test_var_', i)),
                   data.frame(var = paste0(variables[combinaciones_variables[[j]]], collapse = ', '),
                              error = ifelse(get(paste0('errores_', i, '_', j))[[4]]$ind_indepCopula==0, 
                                             round(eval_metric_functions[[eval_metric]]
                                             (get(paste0('errores_', i, '_', j))[[3]]$Target,
                                               get(paste0('errores_', i, '_', j))[[3]]$PREDICCION,
                                               get(paste0('errores_', i, '_', j))[[3]]$ERROR),5),
                                             Inf)
                   )
             )
      )
      
    }
    
    chequeo_errores <- get(paste0('errores_train_var_', i))$error[!duplicated(get(paste0('errores_train_var_', i))$error)]
    if (length(chequeo_errores)==1){
      if (chequeo_errores==Inf) { 
        if (i == 1){
          stop('Las variables son un ?ordo')
          i <- num_iter + 1
          next
        } else {
          var_quitadas <- c(var_quitadas, 
                            as.character(get(paste0('errores_train_var_', (i - 1)))[get(paste0('errores_train_var_', (i - 1)))$error==errores_train$error[i],'var']))
          if (length(var_quitadas)==length(combinaciones_variables)){
            tablas_output <- genera_mejor_iter(pasos_stepwise, pred_train, pred_valid, pred_test, modelo)
            if ((verbosity) & (nrow(pasos_stepwise[pasos_stepwise$paso==iteracion,])>0)){
              print(pasos_stepwise[pasos_stepwise$paso==iteracion,])
            }
            i <- num_iter + 1
            next
          }
          
          orden_var_ant <- get(paste0('errores_train_var_', (i-1)))[order(get(paste0('errores_train_var_', (i-1)))$error),]
          copula_stepwise_ant <- which(get(paste0('errores_train_var_', (i-1)))$var==orden_var_ant[length(var_quitadas) + 1 - num_inf,'var'])
          datos_train_fija <- get(paste0('errores_', (i-1), '_', copula_stepwise_ant))[[1]]
          datos_valid_fija <- get(paste0('errores_', (i-1), '_', copula_stepwise_ant))[[2]]
          datos_test_fija <- get(paste0('errores_', (i-1), '_', copula_stepwise_ant))[[3]]
          
          pasos_stepwise <- rbind(pasos_stepwise,
                                  data.frame(paso = iteracion,
                                             iteracion = (i - 1),
                                             estado = c('entrando', 
                                                        'saliendo',
                                                        'saliendo',
                                                        'entrando'),
                                             variable = c('',
                                                          '',
                                                          as.character(errores_train$var)[nrow(errores_train)],
                                                          paste0(variables[combinaciones_variables[[copula_stepwise_ant]]], collapse = ',')),
                                             copula = c('',
                                                        '',
                                                        as.character(get(paste0('errores_', (i - 1), '_', which(variables_cruce %in% as.character(errores_train$var)[nrow(errores_train)])))[[4]]$mejor_copula_var),
                                                        as.character(get(paste0('errores_', (i - 1), '_', copula_stepwise_ant))[[4]]$mejor_copula_var)),
                                             error_train = c(Inf,
                                                             NA,
                                                             NA,
                                                             get(paste0('errores_train_var_', (i - 1)))[ copula_stepwise_ant, 'error']),
                                             error_valid = c(Inf,
                                                             NA,
                                                             NA,
                                                             get(paste0('errores_valid_var_', (i - 1)))[ copula_stepwise_ant, 'error']),
                                             error_test = c(Inf,
                                                            NA,
                                                            NA,
                                                            get(paste0('errores_test_var_', (i - 1)))[ copula_stepwise_ant, 'error'])
                                             
                                  )
          )
          
          errores_train <- errores_train[-nrow(errores_train),]
          errores_valid <- errores_valid[-nrow(errores_valid),]
          errores_test <- errores_test[-nrow(errores_test),]
          errores_train <- rbind(errores_train,
                                 data.frame(iter = (i-1),
                                            error = get(paste0('errores_train_var_', (i-1)))[ copula_stepwise_ant, 'error'],
                                            var = paste0(variables[combinaciones_variables[[copula_stepwise_ant]]], collapse = ',')
                                 )
          )
          errores_valid <- rbind(errores_valid,
                                 data.frame(iter = (i-1),
                                            error = get(paste0('errores_valid_var_', (i-1)))[ copula_stepwise_ant, 'error'],
                                            var = paste0(variables[combinaciones_variables[[copula_stepwise_ant]]], collapse = ',')
                                 )
          )
          errores_test <- rbind(errores_test,
                                data.frame(iter = (i-1),
                                           error = get(paste0('errores_test_var_', (i-1)))[ copula_stepwise_ant, 'error'],
                                           var = paste0(variables[combinaciones_variables[[copula_stepwise_ant]]], collapse = ',')
                                )
          )
          
          modelo[['iteraciones']][[(i - 1)]][['final']][['copula']] <- get(paste0('errores_', (i - 1), '_', copula_stepwise_ant))[[5]]
          modelo[['iteraciones']][[(i - 1)]][['final']][['inf_iter']] <- get(paste0('errores_', (i - 1), '_', copula_stepwise_ant))[[6]]
          modelo[['iteraciones']][[(i - 1)]][['final']][['aprox_variables']] <- get(paste0('errores_', (i - 1), '_', copula_stepwise_ant))[[7]]
          
          
          if ((errores_train$error[nrow(errores_train)-1] <= errores_train$error[nrow(errores_train)]) & (i > 2)){
            tablas_output <- genera_mejor_iter(pasos_stepwise, pred_train, pred_valid, pred_test, modelo)
            if ((verbosity) & (nrow(pasos_stepwise[pasos_stepwise$paso==iteracion,])>0)){
              print(pasos_stepwise[pasos_stepwise$paso==iteracion,])
            }
            i <- num_iter + 1
            next
          }
          
          if ((verbosity) & (nrow(pasos_stepwise[pasos_stepwise$paso==iteracion,])>0)){
            print(pasos_stepwise[pasos_stepwise$paso==iteracion,])
          }
          iteracion <- iteracion + 1
          next
        }
      }
    }
    
    copula_stepwise <- which.min(get(paste0('errores_train_var_', i))$error)
    
    error_anterior <- errores_train[errores_train$iter == (i - 1), 'error']
    
    if ((get(paste0('errores_train_var_', i))[ copula_stepwise, 'error'] >= error_anterior) & (i > 1)){
      var_quitadas <- c(var_quitadas, 
                        as.character(get(paste0('errores_train_var_', (i - 1)))[get(paste0('errores_train_var_', (i - 1)))$error==errores_train$error[i],'var']))
      if (length(var_quitadas)==length(combinaciones_variables)){
        tablas_output <- genera_mejor_iter(pasos_stepwise, pred_train, pred_valid, pred_test, modelo)
        if ((verbosity) & (nrow(pasos_stepwise[pasos_stepwise$paso==iteracion,])>0)){
          print(pasos_stepwise[pasos_stepwise$paso==iteracion,])
        }
        i <- num_iter + 1
        next
      }
      
      orden_var_ant <- get(paste0('errores_train_var_', (i-1)))[order(get(paste0('errores_train_var_', (i-1)))$error),]
      copula_stepwise_ant <- which(get(paste0('errores_train_var_', (i-1)))$var==orden_var_ant[length(var_quitadas) + 1 - num_inf,'var'])
      datos_train_fija <- get(paste0('errores_', (i-1), '_', copula_stepwise_ant))[[1]]
      datos_valid_fija <- get(paste0('errores_', (i-1), '_', copula_stepwise_ant))[[2]]
      datos_test_fija <- get(paste0('errores_', (i-1), '_', copula_stepwise_ant))[[3]]
      
      pasos_stepwise <- rbind(pasos_stepwise,
                              data.frame(paso = iteracion,
                                         iteracion = (i - 1),
                                         estado = c('entrando', 
                                                    'saliendo',
                                                    'saliendo',
                                                    'entrando'),
                                         variable = c(paste0(variables[combinaciones_variables[[copula_stepwise]]], collapse = ', '),
                                                      paste0(variables[combinaciones_variables[[copula_stepwise]]], collapse = ', '),
                                                      as.character(errores_train$var)[nrow(errores_train)],
                                                      paste0(variables[combinaciones_variables[[copula_stepwise_ant]]], collapse = ', ')),
                                         copula = c(as.character(get(paste0('errores_', i, '_', copula_stepwise))[[4]]$mejor_copula_var),
                                                    as.character(get(paste0('errores_', i, '_', copula_stepwise))[[4]]$mejor_copula_var),
                                                    as.character(get(paste0('errores_', (i - 1), '_', which(variables_cruce %in% as.character(errores_train$var)[nrow(errores_train)])))[[4]]$mejor_copula_var),
                                                    as.character(get(paste0('errores_', (i - 1), '_', copula_stepwise_ant))[[4]]$mejor_copula_var)),
                                         error_train = c(get(paste0('errores_train_var_', i))[ copula_stepwise, 'error'],
                                                         NA,
                                                         NA,
                                                         get(paste0('errores_train_var_', (i - 1)))[ copula_stepwise_ant, 'error']),
                                         error_valid = c(get(paste0('errores_valid_var_', i))[ copula_stepwise, 'error'],
                                                         NA,
                                                         NA,
                                                         get(paste0('errores_valid_var_', (i - 1)))[ copula_stepwise_ant, 'error']),
                                         error_test = c(get(paste0('errores_test_var_', i))[ copula_stepwise, 'error'],
                                                        NA,
                                                        NA,
                                                        get(paste0('errores_test_var_', (i - 1)))[ copula_stepwise_ant, 'error'])
                                         
                              )
      )
      
      errores_train <- errores_train[-nrow(errores_train),]
      errores_valid <- errores_valid[-nrow(errores_valid),]
      errores_test <- errores_test[-nrow(errores_test),]
      errores_train <- rbind(errores_train,
                             data.frame(iter = (i-1),
                                        error = get(paste0('errores_train_var_', (i-1)))[ copula_stepwise_ant, 'error'],
                                        var = paste0(variables[combinaciones_variables[[copula_stepwise_ant]]], collapse = ', ')
                             )
      )
      errores_valid <- rbind(errores_valid,
                             data.frame(iter = (i-1),
                                        error = get(paste0('errores_valid_var_', (i-1)))[ copula_stepwise_ant, 'error'],
                                        var = paste0(variables[combinaciones_variables[[copula_stepwise_ant]]], collapse = ', ')
                             )
      )
      errores_test <- rbind(errores_test,
                            data.frame(iter = (i-1),
                                       error = get(paste0('errores_test_var_', (i-1)))[ copula_stepwise_ant, 'error'],
                                       var = paste0(variables[combinaciones_variables[[copula_stepwise_ant]]], collapse = ', ')
                            )
      )
      
      modelo[['iteraciones']][[(i - 1)]][['final']][['copula']] <- get(paste0('errores_', (i - 1), '_', copula_stepwise_ant))[[5]]
      modelo[['iteraciones']][[(i - 1)]][['final']][['inf_iter']] <- get(paste0('errores_', (i - 1), '_', copula_stepwise_ant))[[6]]
      modelo[['iteraciones']][[(i - 1)]][['final']][['aprox_variables']] <- get(paste0('errores_', (i - 1), '_', copula_stepwise_ant))[[7]]
      
      if ((errores_train$error[nrow(errores_train)-1] <= errores_train$error[nrow(errores_train)]) & (i > 2)){
        tablas_output <- genera_mejor_iter(pasos_stepwise, pred_train, pred_valid, pred_test, modelo)
        if ((verbosity) & (nrow(pasos_stepwise[pasos_stepwise$paso==iteracion,])>0)){
          print(pasos_stepwise[pasos_stepwise$paso==iteracion,])
        }
        i <- num_iter + 1
        next
      }
      if ((verbosity) & (nrow(pasos_stepwise[pasos_stepwise$paso==iteracion,])>0)){
        print(pasos_stepwise[pasos_stepwise$paso==iteracion,])
      }
      iteracion <- iteracion + 1
      next
    } else {
      pasos_stepwise <- rbind(pasos_stepwise,
                              data.frame(paso = iteracion,
                                         iteracion = i,
                                         estado = 'entrando',
                                         variable = paste0(variables[combinaciones_variables[[copula_stepwise]]], collapse = ', '),
                                         copula = as.character(get(paste0('errores_', i, '_', copula_stepwise))[[4]]$mejor_copula_var),
                                         error_train = get(paste0('errores_train_var_', i))[ copula_stepwise, 'error'],
                                         error_valid = get(paste0('errores_valid_var_', i))[ copula_stepwise, 'error'],
                                         error_test = get(paste0('errores_test_var_', i))[ copula_stepwise, 'error']
                              )
      )
      var_quitadas <- as.character(get(paste0('errores_train_var_', i))[get(paste0('errores_train_var_', i))$error==Inf, 'var'])
      num_inf <- length(var_quitadas)
    
      datos_train_fija <- get(paste0('errores_', i, '_', copula_stepwise))[[1]]
      datos_valid_fija <- get(paste0('errores_', i, '_', copula_stepwise))[[2]]
      datos_test_fija <- get(paste0('errores_', i, '_', copula_stepwise))[[3]]
      pred_train[[i]] <- datos_train_fija$PREDICCION
      pred_valid[[i]] <- datos_valid_fija$PREDICCION
      pred_test[[i]] <- datos_test_fija$PREDICCION
      modelo[['iteraciones']][[i]] <- list()
      modelo[['iteraciones']][[i]][['original']][['copula']] <- get(paste0('errores_', i, '_', copula_stepwise))[[5]]
      modelo[['iteraciones']][[i]][['original']][['inf_iter']] <- get(paste0('errores_', i, '_', copula_stepwise))[[6]]
      modelo[['iteraciones']][[i]][['original']][['aprox_variables']] <- get(paste0('errores_', i, '_', copula_stepwise))[[7]]
      modelo[['iteraciones']][[i]][['final']][['copula']] <- get(paste0('errores_', i, '_', copula_stepwise))[[5]]
      modelo[['iteraciones']][[i]][['final']][['inf_iter']] <- get(paste0('errores_', i, '_', copula_stepwise))[[6]]
      modelo[['iteraciones']][[i]][['final']][['aprox_variables']] <- get(paste0('errores_', i, '_', copula_stepwise))[[7]]
      
      errores_train <- rbind(errores_train,
                             data.frame(iter = i,
                                        error = get(paste0('errores_train_var_', i))[ copula_stepwise, 'error'],
                                        var = paste0(variables[combinaciones_variables[[copula_stepwise]]], collapse = ', ')
                             )
      )
      errores_valid <- rbind(errores_valid,
                             data.frame(iter = i,
                                        error = get(paste0('errores_valid_var_', i))[ copula_stepwise, 'error'],
                                        var = paste0(variables[combinaciones_variables[[copula_stepwise]]], collapse = ', ')
                             )
      )
      errores_test <- rbind(errores_test,
                            data.frame(iter = i,
                                       error = get(paste0('errores_test_var_', i))[ copula_stepwise, 'error'],
                                       var = paste0(variables[combinaciones_variables[[copula_stepwise]]], collapse = ', ')
                            )
      )
      
      if (i == num_iter){
        tablas_output <- genera_mejor_iter(pasos_stepwise, pred_train, pred_valid, pred_test, modelo)
        if ((verbosity) & (nrow(pasos_stepwise[pasos_stepwise$paso==iteracion,])>0)){
          print(pasos_stepwise[pasos_stepwise$paso==iteracion,])
        }
        i <- num_iter + 1
        next
      }
      
      if ((verbosity) & (nrow(pasos_stepwise[pasos_stepwise$paso==iteracion,])>0)){
        print(pasos_stepwise[pasos_stepwise$paso==iteracion,])
      }
      i <- i + 1
      iteracion <- iteracion + 1
    }
  }
  return(tablas_output)
}

