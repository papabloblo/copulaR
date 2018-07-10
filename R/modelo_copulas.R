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
  
  
  

  # ERRORES -----------------------------------------------------------------
  
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
  

  colnames(train)[which(colnames(train) == target)] <- "Target"
  colnames(valid)[which(colnames(valid) == target)] <- "Target"
  if (target %in% colnames(test)){
    colnames(test)[which(colnames(test) == target)] <- "Target"
  } else {
    test$Target <- NA
  }
  
  variables <- names(train)[which(names(train)!="Target")]
  num_variables <- length(variables)
  max_dim_copulas <- 2 
  
  errores_train <- data.frame(iter = 0,
                              error = 0,
                              var = "")
  errores_valid <- data.frame(iter = 0,
                              error = 0,
                              var = "")
  errores_test <- data.frame(iter=0,
                             error = 0,
                             var = "")
  pasos_stepwise <- data.frame()
  iteracion <- 1
  combinaciones_variables <- genera_combinaciones_variables(num_variables,
                                                            max_dim_copulas)
  variables_cruce <- c()
  for (i in 1:length(combinaciones_variables)){
    variables_cruce <- c(variables_cruce,
                         paste0(variables[combinaciones_variables[[i]]], collapse = ', '))
  }
  i <- 1
  pred_train <- list()
  pred_valid <- list()
  pred_test <- list()
  
  modelo <- list()
  modelo[['train']] <- train
  modelo[['max_bins']] <- max_bins
  modelo[['num_sim']] <- num_sim
  modelo[['bin_target']] <- bin_target
  modelo[['iteraciones']] <- list()
  
  while (i <= num_iter){
    
    if (i > 1){
      if (early_stopping_round>0){
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
        datos_train <- train
        datos_valid <- valid
        datos_test <- test
        
        datos_train$PREDICCION <- mean(datos_train$Target)
        datos_valid$PREDICCION <- mean(datos_train$Target)
        datos_test$PREDICCION <- mean(datos_train$Target)
        datos_train$ERROR <- (datos_train$Target - datos_train$PREDICCION)/datos_train$Target
        datos_valid$ERROR <- (datos_valid$Target - datos_valid$PREDICCION)/datos_valid$Target
        datos_test$ERROR <- (datos_test$Target - datos_test$PREDICCION)/datos_test$Target
        
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
        
      } else  {
        
        datos_train <- datos_train_fija
        datos_valid <- datos_valid_fija
        datos_test <- datos_test_fija
      }
      
      assign(paste0('errores_', i, '_', j), 
             ajuste_var_cop(datos_train,
                            datos_valid,
                            datos_test,
                            variables[combinaciones_variables[[j]]],
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

###KDD98##
train <- read.csv('targetContinuoKDDTrain.csv', 
                  header = TRUE,
                  stringsAsFactors = FALSE)

test <- read.csv('Kdd1988targetcontinuoTest.csv', 
                 header = TRUE,
                 stringsAsFactors = FALSE)

train$TargetD <- as.numeric(substr(train$TargetD,2, nchar(train$TargetD)))
train$GiftAvgLast <- as.numeric(substr(train$GiftAvgLast,2, nchar(train$GiftAvgLast)))
train$GiftAvg36 <- as.numeric(substr(train$GiftAvg36,2, nchar(train$GiftAvg36)))
train$GiftAvgAll <- as.numeric(substr(train$GiftAvgAll,2, nchar(train$GiftAvgAll)))
var_aux <- substr(train$DemMedHomeValue,2, nchar(train$DemMedHomeValue))
train$DemMedHomeValue <- as.numeric(gsub(",", "", var_aux))
var_aux <- substr(train$DemMedIncome,2, nchar(train$DemMedIncome))
train$DemMedIncome <- as.numeric(gsub(",", "", var_aux))

test$TargetD <- as.numeric(substr(test$TargetD,2, nchar(test$TargetD)))
test$GiftAvgLast <- as.numeric(substr(test$GiftAvgLast,2, nchar(test$GiftAvgLast)))
test$GiftAvg36 <- as.numeric(substr(test$GiftAvg36,2, nchar(test$GiftAvg36)))
test$GiftAvgAll <- as.numeric(substr(test$GiftAvgAll,2, nchar(test$GiftAvgAll)))
var_aux <- substr(test$DemMedHomeValue,2, nchar(test$DemMedHomeValue))
test$DemMedHomeValue <- as.numeric(gsub(",", "", var_aux))
var_aux <- substr(test$DemMedIncome,2, nchar(test$DemMedIncome))
test$DemMedIncome <- as.numeric(gsub(",", "", var_aux))

ind_train <- sample(1:nrow(train), 0.7*nrow(train))
train_tab <- train[ind_train,]
valid_tab <- train[-ind_train,]

### tablas communities##
train <- read.csv('communitiestest.csv', 
                  header = TRUE,
                  stringsAsFactors = FALSE)
train <- train[-90,]

ind_train <- sample(1:nrow(train), 0.7*nrow(train))
train_tab <- train[ind_train,]
valid_tab <- train[-ind_train,]

test <- read.csv('communitiesTraining.csv', 
                 header = TRUE,
                 stringsAsFactors = FALSE)

colnames(train_tab)[ncol(train_tab)] <- 'TargetD'
colnames(valid_tab)[ncol(valid_tab)] <- 'TargetD'
colnames(test)[ncol(test)] <- 'TargetD'

###ailerons##
train_tab <- read.csv('aileronsTrain.csv', 
                      header = TRUE,
                      stringsAsFactors = FALSE)

valid_tab <- read.csv('aileronsValidate.csv', 
                      header = TRUE,
                      stringsAsFactors = FALSE)

test <- read.csv('aileronsTestSinTarget.csv', 
                 header = TRUE,
                 sep = ";",
                 stringsAsFactors = FALSE)

colnames(train_tab)[ncol(train_tab)] <- 'TargetD'
colnames(valid_tab)[ncol(valid_tab)] <- 'TargetD'
test$TargetD <- NA

### elevators ####
train <- read.csv('elevators.csv', 
                  header = TRUE,
                  sep = ";",
                  stringsAsFactors = FALSE)

ind_train <- sample(1:nrow(train), 0.6*nrow(train))
train_tab <- train[ind_train,]
valid_tab <- train[-ind_train,]

ind_test <- sample(1:nrow(valid_tab), 0.5*nrow(valid_tab))
test <- valid_tab[ind_test,]
valid_tab <- valid_tab[-ind_test,] 

colnames(train_tab)[ncol(train_tab)] <- 'TargetD'
colnames(valid_tab)[ncol(valid_tab)] <- 'TargetD'
colnames(test)[ncol(test)] <- 'TargetD'

modelo <- copula.model(train = train_tab,
                  valid = valid_tab,
                  target = "TargetD",
                  test = NULL,
                  num_iter = 350,
                  num_sim = 550,
                  max_bins = 30,
                  bin_target = FALSE,
                  eval_metric = "SMAPE",
                  num_obs_fit = 350,
                  early_stopping_round = 2)

predict.copula <- function(score = NULL,
                           modelo = NULL){
  if (is.null(score)){
    stop('sigue asi')
  } else if (!is.data.frame(score)){
    stop('cambia eso hombre')
  }
  
  names_modelo <- c("errores_train",
                    "errores_valid",
                    "errores_test",
                    "pasos_stepwise",
                    "pred_train",
                    "pred_valid",
                    "pred_test",
                    "modelo")
  
  if (is.null(modelo)){
    stop('un modelo hombre')
  } else if (!(all(length(sort(names(modelo)))==
                   length(sort(names_modelo))) &
               all(sort(names(modelo))==
                   sort(names_modelo)))){
    stop('un modelo decente')
  }
  
  if (!(all(length(sort(names(modelo$modelo$train)[names(modelo$modelo$train)!="Target"]))==
            length(sort(names(score)))) &
        all(sort(names(modelo$modelo$train)[names(modelo$modelo$train)!="Target"])==
            sort(names(score))))){
    stop('que variables pasas')
  }
  
  modelo$modelo$train$PREDICCION <- mean(modelo$modelo$train$Target)
  modelo$modelo$train$ERROR <- (modelo$modelo$train$Target - modelo$modelo$train$PREDICCION)/
    modelo$modelo$train$Target
  score$PREDICCION <- mean(modelo$modelo$train$Target)
  
  for (i in 1:length(modelo$modelo$iteraciones)){
    
    dim_iter <- ncol(modelo$modelo$iteraciones[[i]]$inf_iter) - 1
    var_iter <- substr(colnames(modelo$modelo$iteraciones[[i]]$inf_iter)[1:dim_iter],
                       1, 
                       nchar(colnames(modelo$modelo$iteraciones[[i]]$inf_iter)[1]) - 5)
    
    valores_scores <- data.frame(score[!duplicated(score[,var_iter]),var_iter])
    colnames(valores_scores) <- var_iter
    
    if (is.null(modelo$modelo$max_bins)){
      
      for (j in 1:dim_iter){
        if (j == 1){
          coincidencias <- data.frame(apply(as.matrix(valores_scores[,j]),
                                 1,
                                 function(x){x %in% modelo$modelo$iteraciones[[i]]$inf_iter[,j]}))
        } else {
          coincidencias <- cbind(coincidencias,
                                 apply(as.matrix(valores_scores[,j]),
                                       1,
                                       function(x){x %in% modelo$modelo$iteraciones[[i]]$inf_iter[,j]}))
        }
      }
      
      valores_var_nuevos <- data.frame(valores_scores[apply(coincidencias,
                                  1,
                                  function(x){sum(x)!=dim_iter}),])
      
      if (nrow(valores_var_nuevos)>0){
        colnames(valores_var_nuevos) <- var_iter
        
        variables <- c(var_iter, 'ERROR')
        
        train_var <- modelo$modelo$train %>% 
          select_(.dots = variables) 
        train_var <- train_var[!duplicated(train_var),]
        
        valores_var_nuevos[,paste0(colnames(valores_var_nuevos)[1:dim_iter], '_hist')] <- 
          valores_var_nuevos[,colnames(valores_var_nuevos)[1:dim_iter]]
          
        resultados <- data.frame()
        
        n <- nrow(valores_var_nuevos)
        
        if ((n*modelo$modelo$num_sim)>=1000000){
          num_iter <- floor((n*modelo$modelo$num_sim)/1000000) + 1
          fila_ini <- 1
          fila_fin <- min(c(floor(fila_ini + (1000000/modelo$modelo$num_sim)),n))
          for (j in 1:num_iter){
            train_aux <- valores_var_nuevos[fila_ini:fila_fin,]
            resultados_aux <- puntuacion_copula_opt(datos_iter = data.frame(train_aux),
                                                    n.ventas = modelo$modelo$num_sim,
                                                    copulaoptima=modelo$modelo$iteraciones[[i]]$copula,
                                                    train =  train_var)
            resultados <- rbind(resultados, resultados_aux)
            fila_ini <- fila_fin + 1
            fila_fin <- min(c(floor(fila_ini + (1000000/modelo$modelo$num_sim)),n))
          }
        } else {
          resultados <- puntuacion_copula_opt(datos_iter = valores_var_nuevos,
                                              n.ventas = modelo$modelo$num_sim,
                                              copulaoptima=modelo$modelo$iteraciones[[i]]$copula,
                                              train =  train_var)
        }
        
        info_iter <-modelo$modelo$iteraciones[[i]]$inf_iter
        colnames(info_iter)[1:dim_iter] <- var_iter
        names(resultados)[names(resultados) == 'ERROR'] <- 'ERROR_COP'
        
        resultados2 <- rbind(resultados[,c(var_iter, "ERROR_COP")],
                             info_iter)
        
        resultados3 <- modelo$modelo$train %>% left_join(resultados2, by = var_iter)
        resultados3.5 <- score %>% left_join(resultados2, by = var_iter)
        resultados3$pred_nueva <- resultados3$PREDICCION/(1-resultados3$ERROR_COP)
        resultados3$nuevo_error <- (resultados3$Target-resultados3$pred_nueva)/resultados3$Target
        resultados3.5$pred_nueva <- resultados3.5$PREDICCION/(1-resultados3.5$ERROR_COP)
        modelo$modelo$train$PREDICCION <- resultados3$pred_nueva
        modelo$modelo$train$ERROR <- resultados3$nuevo_error
        score$PREDICCION <- resultados3.5$pred_nueva
      } else {
        
        info_iter <- modelo$modelo$iteraciones[[i]]$inf_iter
        colnames(info_iter)[1:dim_iter] <- var_iter
        
        resultados2 <- info_iter
        
        resultados3 <- modelo$modelo$train %>% left_join(resultados2, by = var_iter)
        resultados3.5 <- score %>% left_join(resultados2, by = var_iter)
        resultados3$pred_nueva <- resultados3$PREDICCION/(1-resultados3$ERROR_COP)
        resultados3$nuevo_error <- (resultados3$Target-resultados3$pred_nueva)/resultados3$Target
        resultados3.5$pred_nueva <- resultados3.5$PREDICCION/(1-resultados3.5$ERROR_COP)
        modelo$modelo$train$PREDICCION <- resultados3$pred_nueva
        modelo$modelo$train$ERROR <- resultados3$nuevo_error
        score$PREDICCION <- resultados3.5$pred_nueva
        
      }
    } else {
      
      bins_var <- modelo$modelo$iteraciones[[i]]$aprox_variables
      
      variables <- c(var_iter, 'ERROR')
      
      train_var <- modelo$modelo$train %>% 
        select_(.dots = variables) 
      train_var <- train_var[!duplicated(train_var),]
      
      for (j in 1:(length(var_iter))){
        valores_scores[,paste0(colnames(valores_scores)[j], '_hist')] <- apply(as.matrix(valores_scores[,j]),
                                                                 1,
                                                                 function(x){
                                                                   bins_var[[j]][which.min(abs(bins_var[[j]] - x))]
                                                                 })
        train_var[,j] <- apply(as.matrix(train_var[,j]),
                               1,
                               function(x){
                                 bins_var[[j]][which.min(abs(bins_var[[j]] - x))]
                               })
      }
      
      if (modelo$modelo$bin_target){
        train_var[,length(variables)] <- apply(as.matrix(train_var[,length(variables)]),
                                               1,
                                               function(x){
                                                 bins_var[[length(variables)]][which.min(abs(bins_var[[length(variables)]] - x))]
                                               })
        
      }
      
      for (j in 1:dim_iter){
        if (j == 1){
          coincidencias <- data.frame(apply(as.matrix(valores_scores[,paste0(colnames(valores_scores)[j], '_hist')]),
                                            1,
                                            function(x){x %in% modelo$modelo$iteraciones[[i]]$inf_iter[,j]}))
        } else {
          coincidencias <- cbind(coincidencias,
                                 apply(as.matrix(valores_scores[,paste0(colnames(valores_scores)[j], '_hist')]),
                                       1,
                                       function(x){x %in% modelo$modelo$iteraciones[[i]]$inf_iter[,j]}))
        }
      }
      
      valores_var_nuevos <- data.frame(valores_scores[apply(coincidencias,
                                                   1,
                                                   function(x){sum(x)!=dim_iter}),])
      
      if (nrow(valores_var_nuevos)>0){
        
        resultados <- data.frame()
        
        n <- nrow(valores_var_nuevos)
        
        if ((n*modelo$modelo$num_sim)>=1000000){
          num_iter <- floor((n*modelo$modelo$num_sim)/1000000) + 1
          fila_ini <- 1
          fila_fin <- min(c(floor(fila_ini + (1000000/modelo$modelo$num_sim)),n))
          for (j in 1:num_iter){
            train_aux <- valores_var_nuevos[fila_ini:fila_fin,]
            resultados_aux <- puntuacion_copula_opt(datos_iter = data.frame(train_aux),
                                                    n.ventas = modelo$modelo$num_sim,
                                                    copulaoptima=modelo$modelo$iteraciones[[i]]$copula,
                                                    train =  train_var)
            resultados <- rbind(resultados, resultados_aux)
            fila_ini <- fila_fin + 1
            fila_fin <- min(c(floor(fila_ini + (1000000/modelo$modelo$num_sim)),n))
          }
        } else {
          resultados <- puntuacion_copula_opt(datos_iter = valores_var_nuevos,
                                              n.ventas = modelo$modelo$num_sim,
                                              copulaoptima=modelo$modelo$iteraciones[[i]]$copula,
                                              train =  train_var)
        }
        
        info_iter <-modelo$modelo$iteraciones[[i]]$inf_iter
        colnames(train_var)[1:dim_iter] <- paste0(colnames(train_var)[1:dim_iter], '_hist')
        train_var[,var_iter] <- modelo$modelo$train[!duplicated(modelo$modelo$train[,variables]),var_iter]
        info_iter <- rbind(train_var[!duplicated(train_var[,colnames(valores_scores)]),colnames(valores_scores)],
                           valores_scores[apply(coincidencias,
                                          1,
                                          function(x){sum(x)==dim_iter}),]) %>% left_join(info_iter, colnames(info_iter)[1:dim_iter])
        info_iter <- info_iter[!duplicated(info_iter),]
        names(resultados)[names(resultados) == 'ERROR'] <- 'ERROR_COP'
        
        resultados2 <- rbind(resultados[,c(var_iter, "ERROR_COP")],
                             info_iter[,c(var_iter, "ERROR_COP")])
        
        resultados3 <- modelo$modelo$train %>% left_join(resultados2, by = var_iter)
        resultados3.5 <- score %>% left_join(resultados2, by = var_iter)
        resultados3$pred_nueva <- resultados3$PREDICCION/(1-resultados3$ERROR_COP)
        resultados3$nuevo_error <- (resultados3$Target-resultados3$pred_nueva)/resultados3$Target
        resultados3.5$pred_nueva <- resultados3.5$PREDICCION/(1-resultados3.5$ERROR_COP)
        modelo$modelo$train$PREDICCION <- resultados3$pred_nueva
        modelo$modelo$train$ERROR <- resultados3$nuevo_error
        score$PREDICCION <- resultados3.5$pred_nueva
        
      } else {
        
        info_iter <-modelo$modelo$iteraciones[[i]]$inf_iter
        colnames(train_var)[1:dim_iter] <- paste0(colnames(train_var)[1:dim_iter], '_hist')
        train_var[,var_iter] <- modelo$modelo$train[!duplicated(modelo$modelo$train[,variables]),var_iter]
        info_iter <- rbind(train_var[!duplicated(train_var[,colnames(valores_scores)]),colnames(valores_scores)],
                           valores_scores[apply(coincidencias,
                                                1,
                                                function(x){sum(x)==dim_iter}),]) %>% left_join(info_iter, colnames(info_iter)[1:dim_iter])
        info_iter <- info_iter[!duplicated(info_iter),]
        
        resultados2 <- info_iter
        
        resultados3 <- modelo$modelo$train %>% left_join(resultados2, by = var_iter)
        resultados3.5 <- score %>% left_join(resultados2, by = var_iter)
        resultados3$pred_nueva <- resultados3$PREDICCION/(1-resultados3$ERROR_COP)
        resultados3$nuevo_error <- (resultados3$Target-resultados3$pred_nueva)/resultados3$Target
        resultados3.5$pred_nueva <- resultados3.5$PREDICCION/(1-resultados3.5$ERROR_COP)
        modelo$modelo$train$PREDICCION <- resultados3$pred_nueva
        modelo$modelo$train$ERROR <- resultados3$nuevo_error
        score$PREDICCION <- resultados3.5$pred_nueva
        
      }
    }
  }
  return(score$PREDICCION)
}

score <- train_tab
score$TargetD <- NULL
punt <- predict.copula(score = score,
                       modelo = modelo_aux)
