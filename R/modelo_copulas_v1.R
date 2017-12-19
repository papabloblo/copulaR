library(ggplot2)
#library(plotly)
library(DiagTest3Grp)
library(copula)
library(magrittr)
library(dplyr)
library(data.table)
library(zoo)
library(writexl)
library(BBmisc)

setwd('C:/Users/juan.carrillo/Desktop/Articulos/articulo_copulas/')
#Sys.setlocale("LC_ALL","English")
source('copula_optima.R')
source('puntuacion_copulas.R')
source('genera_combinaciones_variables.R')

###KDD98##
train <- read.csv('targetContinuoKDDTrain.csv', 
                  header = TRUE,
                  stringsAsFactors = FALSE)

test <- read.csv('Kdd1988targetcontinuoTest.csv', 
                 header = TRUE,
                 stringsAsFactors = FALSE)


### tablas communities##
train <- read.csv('communitiestest.csv', 
                  header = TRUE,
                  stringsAsFactors = FALSE)
train <- train[-90,]

test <- read.csv('communitiesTraining.csv', 
                 header = TRUE,
                 stringsAsFactors = FALSE)

colnames(train)[ncol(train)] <- 'TargetD'
colnames(test)[ncol(test)] <- 'TargetD'

head(train)
names(train)
head(test)

#Tratamineto_datos

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

# #Gráfico de la serie
# p <- ggplot(data = train, aes(x = FECHA)) + 
#   geom_line(aes(y = DEMANDA_GAS), color = 'steelblue') +  
#   geom_line(aes(y = PREDICCION_ARIMA), color = 'firebrick')+
#   scale_x_date(date_breaks = '1 month',
#                #date_minor_breaks='1 months',
#                date_labels = "%b")
# ggplotly(p)
# 
# ggplot(data = train, aes(x = ERROR_ARIMA,
#                          y = INC_TEMPERATURA_MAXIMA_ABS)) + 
#   geom_point(color = 'steelblue',
#              alpha = 0.5,
#              size = 4)


#Ajuste cópula

puntua_cop <- function(datos_train, datos_test, variables2, precision){
  
  errores <- list()
  
  train2 <- datos_train %>% 
    select_(.dots = variables2) %>% unique()
  
  train3 <- datos_test %>%
    select_(.dots = variables2) %>% unique()
  
  train2 <- unique(rbind(train2,train3))
  
  variables <- c(variables2, 'ERROR')
  
  train_var <- datos_train %>% 
    select_(.dots = variables) %>% unique()
  
  ini <- Sys.time()
  
  mejor_copula <- train_var %>% 
    copula.optima()
  
   #datos_iter <- train2[1,]
#     precision <- 600
# copulaoptima=mejor_copula$copulaoptima
#    train =  train_var
#print(mejor_copula$indep)
  
  # resultados <- apply(train2, 
  #                     1,
  #                     function(x){puntuacion_copula_opt(datos_iter = x,
  #                                                       n.ventas = precision,
  #                                                       copulaoptima=mejor_copula$copulaoptima,
  #                                                       train =  train_var)})
  
  if (mejor_copula$indep<=0.05){
    ini <- Sys.time()
    resultados <- data.frame()
    n <- nrow(train2)
    if ((n*precision)>=1000000){
      num_iter <- floor((n*precision)/1000000) + 1
      fila_ini <- 1
      fila_fin <- min(floor(fila_ini + (1000000/precision)),n) 
      for (i in 1:num_iter){
        train_aux <- train2[fila_ini:fila_fin,]
        resultados_aux <- puntuacion_copula_opt(datos_iter = data.frame(train_aux),
                                                n.ventas = precision,
                                                copulaoptima=mejor_copula$copulaoptima,
                                                train =  train_var)
        resultados <- rbind(resultados, resultados_aux)
        fila_ini <- fila_fin + 1
        fila_fin <- min(floor(fila_ini + (1000000/precision)),n)
      }
    } else {
      resultados <- puntuacion_copula_opt(datos_iter = train2,
                                          n.ventas = precision,
                                          copulaoptima=mejor_copula$copulaoptima,
                                          train =  train_var)
    }
    print(Sys.time() - ini)
    
    #resultados2 <- do.call(rbind, resultados)
    resultados2 <- resultados
    names(resultados2)[names(resultados2) == 'ERROR'] <- 'ERROR_COP'
    
    resultados3 <- datos_train %>% left_join(resultados2)
    resultados4 <- datos_test %>% left_join(resultados2)
    resultados3$pred_nueva <- resultados3$PREDICCION/(1-resultados3$ERROR_COP)
    resultados3$nuevo_error <- (resultados3$TargetD-resultados3$pred_nueva)/resultados3$TargetD
    resultados4$pred_nueva <- resultados4$PREDICCION/(1-resultados4$ERROR_COP)
    resultados4$nuevo_error <- (resultados4$TargetD-resultados4$pred_nueva)/resultados4$TargetD
    datos_train$PREDICCION <- resultados3$pred_nueva
    datos_train$ERROR <- resultados3$nuevo_error
    datos_test$PREDICCION <- resultados4$pred_nueva
    datos_test$ERROR <- resultados4$nuevo_error
    
    errores[[1]] <- datos_train
    errores[[2]] <- datos_test
    errores[[3]] <- data.frame(mejor_copula_var = as.character(mejor_copula$aic[1,1]),
                               ind_indepCopula = 0
    )
  } else {
    errores[[1]] <- data.frame()
    errores[[2]] <- data.frame()
    errores[[3]] <- data.frame('',
                               ind_indepCopula = 1
                               )
    
  }
  
  
  return(errores)
}

variables <- names(train)[c(3:20)]
num_variables <- length(variables)
max_dim_copulas <- 3 

errores_train <- data.frame(iter = 0,
                            error = 0,
                            var = "")
errores_test <- data.frame(iter=0,
                           error = 0,
                           var = "")

##############################################################################

i <- 1
num_iter <- 30
pasos_stepwise <- data.frame()
iteracion <- 1
combianciones_variables <- genera_combinaciones_variables(num_variables,
                                                          max_dim_copulas)

while (i <= num_iter){
  
  assign(paste0('errores_train_var_', i), data.frame())
  assign(paste0('errores_test_var_', i), data.frame())
  
  for (j in 1:length(combianciones_variables)){
    
    if (i == 1){
      datos_train <- train
      datos_test <- test
      
      datos_train$PREDICCION <- mean(datos_train$TargetD)
      datos_test$PREDICCION <- mean(datos_train$TargetD)
      datos_train$ERROR <- (datos_train$TargetD - datos_train$PREDICCION)/datos_train$TargetD
      datos_test$ERROR <- (datos_test$TargetD - datos_test$PREDICCION)/datos_test$TargetD
    } else  {
      
      datos_train <- datos_train_fija
      datos_test <- datos_test_fija
    }
    
    assign(paste0('errores_', i, '_', j), 
           puntua_cop(datos_train,
                      datos_test,
                      variables[combianciones_variables[[j]]],
                      750)
    )
    assign(paste0('errores_train_var_', i), 
           rbind(get(paste0('errores_train_var_', i)),
                 data.frame(var = paste(variables[combianciones_variables[[j]]], sep = ', '),
                            error = ifelse(get(paste0('errores_', i, '_', j))[[3]]$ind_indepCopula==0, 
                                           mean(abs(get(paste0('errores_', i, '_', j))[[1]]$ERROR)),
                                           Inf)
                 )
           )
    )
    assign(paste0('errores_test_var_', i),
           rbind(get(paste0('errores_test_var_', i)),
                 data.frame(var = paste(variables[combianciones_variables[[j]]], sep = ', '),
                            error = ifelse(get(paste0('errores_', i, '_', j))[[3]]$ind_indepCopula==0, 
                                           mean(abs(get(paste0('errores_', i, '_', j))[[2]]$ERROR)),
                                           Inf)
                 )
           )
    )
    
    # print(errores_train_var)
    # print(errores_test_var)
  }
  
  chequeo_errores <- unique(get(paste0('errores_train_var_', i))$error)
  if (length(chequeo_errores)==1){
    if (chequeo_errores==Inf) { 
      if (i == 1){
        print('Las variables son un �ordo')
        i <- num_iter + 1
        next
      } else {
        var_quitadas <- c(var_quitadas, as.character(errores_train$var[i]))
        if (length(var_quitadas)==length(combianciones_variables)){
          errores_train <- data.frame()
          errores_test <- data.frame()
          mejor_num_iteraciones <- pasos_stepwise[which.min(pasos_stepwise$error_test), 'iteracion']
          for (j in 1:mejor_num_iteraciones){
            if (j == 1){
              errores_train <- rbind(errores_train,
                                     data.frame(iter = j,
                                                error = pasos_stepwise[1, 'error_train'],
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
          
          print(pasos_stepwise)
          print(errores_train)
          print(errores_test)
          i <- num_iter + 1
          next
        }
        
        orden_var_ant <- get(paste0('errores_train_var_', (i-1)))[order(get(paste0('errores_train_var_', (i-1)))$error),]
        copula_stepwise_ant <- which(get(paste0('errores_train_var_', (i-1)))$var==orden_var_ant[length(var_quitadas) + 1 - num_inf,'var'])
        datos_train_fija <- get(paste0('errores_', (i-1), '_', copula_stepwise_ant))[[1]]
        datos_test_fija <- get(paste0('errores_', (i-1), '_', copula_stepwise_ant))[[2]]
        
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
                                                        variables[combianciones_variables[[copula_stepwise_ant]]]),
                                           copula = c('',
                                                      '',
                                                      as.character(get(paste0('errores_', (i - 1), '_', which(variables %in% as.character(errores_train$var)[nrow(errores_train)])))[[3]]$mejor_copula_var),
                                                      as.character(get(paste0('errores_', (i - 1), '_', copula_stepwise_ant))[[3]]$mejor_copula_var)),
                                           error_train = c(Inf,
                                                           NA,
                                                           NA,
                                                           get(paste0('errores_train_var_', (i - 1)))[ copula_stepwise_ant, 'error']),
                                           error_test = c(Inf,
                                                          NA,
                                                          NA,
                                                          get(paste0('errores_test_var_', (i - 1)))[ copula_stepwise_ant, 'error'])
                                           
                                )
        )
        
        errores_train <- errores_train[-nrow(errores_train),]
        errores_test <- errores_test[-nrow(errores_test),]
        errores_train <- rbind(errores_train,
                               data.frame(iter = (i-1),
                                          error = get(paste0('errores_train_var_', (i-1)))[ copula_stepwise_ant, 'error'],
                                          var = variables[combinaciones_variables[[copula_stepwise_ant]]]
                               )
        )
        errores_test <- rbind(errores_test,
                              data.frame(iter = (i-1),
                                         error = get(paste0('errores_test_var_', (i-1)))[ copula_stepwise_ant, 'error'],
                                         var = variables[combinaciones_variables[[copula_stepwise_ant]]]
                              )
        )
        
        if ((errores_train$error[nrow(errores_train)-1] <= errores_train$error[nrow(errores_train)]) & (i > 2)){
          errores_train <- data.frame()
          errores_test <- data.frame()
          mejor_num_iteraciones <- pasos_stepwise[which.min(pasos_stepwise$error_test), 'iteracion']
          for (j in 1:mejor_num_iteraciones){
            if (j == 1){
              errores_train <- rbind(errores_train,
                                     data.frame(iter = j,
                                                error = pasos_stepwise[1, 'error_train'],
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
          
          print(pasos_stepwise)
          print(errores_train)
          print(errores_test)
        }
        
        iteracion <- iteracion + 1
        print(var_quitadas)
        print(errores_train)
        print(errores_test)
        print(pasos_stepwise)
        next
       
      }
    }
  }
  
  copula_stepwise <- which.min(get(paste0('errores_train_var_', i))$error)
  
  error_anterior <- errores_train[errores_train$iter==(i - 1), 'error']
  
  if ((round(get(paste0('errores_train_var_', i))[ copula_stepwise, 'error'],5) >= round(error_anterior,5)) & (i > 1)){
    var_quitadas <- c(var_quitadas, as.character(errores_train$var[i]))
    if (length(var_quitadas)==length(combianciones_variables)){
      errores_train <- data.frame()
      errores_test <- data.frame()
      mejor_num_iteraciones <- pasos_stepwise[which.min(pasos_stepwise$error_test), 'iteracion']
      for (j in 1:mejor_num_iteraciones){
        if (j == 1){
          errores_train <- rbind(errores_train,
                                 data.frame(iter = j,
                                            error = pasos_stepwise[1, 'error_train'],
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
      
      print(pasos_stepwise)
      print(errores_train)
      print(errores_test)
      i <- num_iter + 1
      next
    }
    
    orden_var_ant <- get(paste0('errores_train_var_', (i-1)))[order(get(paste0('errores_train_var_', (i-1)))$error),]
    copula_stepwise_ant <- which(get(paste0('errores_train_var_', (i-1)))$var==orden_var_ant[length(var_quitadas) + 1 - num_inf,'var'])
    datos_train_fija <- get(paste0('errores_', (i-1), '_', copula_stepwise_ant))[[1]]
    datos_test_fija <- get(paste0('errores_', (i-1), '_', copula_stepwise_ant))[[2]]
    
    pasos_stepwise <- rbind(pasos_stepwise,
                            data.frame(paso = iteracion,
                                       iteracion = (i - 1),
                                       estado = c('entrando', 
                                                  'saliendo',
                                                  'saliendo',
                                                  'entrando'),
                                       variable = c(variables[copula_stepwise],
                                                    variables[copula_stepwise],
                                                    as.character(errores_train$var)[nrow(errores_train)],
                                                    variables[copula_stepwise_ant]),
                                       copula = c(as.character(get(paste0('errores_', i, '_', copula_stepwise))[[3]]$mejor_copula_var),
                                                  as.character(get(paste0('errores_', i, '_', copula_stepwise))[[3]]$mejor_copula_var),
                                                  as.character(get(paste0('errores_', (i - 1), '_', which(variables %in% as.character(errores_train$var)[nrow(errores_train)])))[[3]]$mejor_copula_var),
                                                  as.character(get(paste0('errores_', (i - 1), '_', copula_stepwise_ant))[[3]]$mejor_copula_var)),
                                       error_train = c(get(paste0('errores_train_var_', i))[ copula_stepwise, 'error'],
                                                       NA,
                                                       NA,
                                                       get(paste0('errores_train_var_', (i - 1)))[ copula_stepwise_ant, 'error']),
                                       error_test = c(get(paste0('errores_test_var_', i))[ copula_stepwise, 'error'],
                                                      NA,
                                                      NA,
                                                      get(paste0('errores_test_var_', (i - 1)))[ copula_stepwise_ant, 'error'])
                                       
                            )
    )
    
    errores_train <- errores_train[-nrow(errores_train),]
    errores_test <- errores_test[-nrow(errores_test),]
    errores_train <- rbind(errores_train,
                           data.frame(iter = (i-1),
                                      error = get(paste0('errores_train_var_', (i-1)))[ copula_stepwise_ant, 'error'],
                                      var = variables[copula_stepwise_ant]
                           )
    )
    errores_test <- rbind(errores_test,
                          data.frame(iter = (i-1),
                                     error = get(paste0('errores_test_var_', (i-1)))[ copula_stepwise_ant, 'error'],
                                     var = variables[copula_stepwise_ant]
                          )
    )
    
    if ((errores_train$error[nrow(errores_train)-1] <= errores_train$error[nrow(errores_train)]) & (i > 2)){
      errores_train <- data.frame()
      errores_test <- data.frame()
      mejor_num_iteraciones <- pasos_stepwise[which.min(pasos_stepwise$error_test), 'iteracion']
      for (j in 1:mejor_num_iteraciones){
        if (j == 1){
          errores_train <- rbind(errores_train,
                                 data.frame(iter = j,
                                            error = pasos_stepwise[1, 'error_train'],
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
      
      print(pasos_stepwise)
      print(errores_train)
      print(errores_test)
      i <- num_iter + 1
      next
    }
    iteracion <- iteracion + 1
    print(var_quitadas)
    print(errores_train)
    print(errores_test)
    print(pasos_stepwise)
    next
  } else {
    pasos_stepwise <- rbind(pasos_stepwise,
                            data.frame(paso = iteracion,
                                       iteracion = i,
                                       estado = 'entrando',
                                       variable = variables[copula_stepwise],
                                       copula = as.character(get(paste0('errores_', i, '_', copula_stepwise))[[3]]$mejor_copula_var),
                                       error_train = get(paste0('errores_train_var_', i))[ copula_stepwise, 'error'],
                                       error_test = get(paste0('errores_test_var_', i))[ copula_stepwise, 'error']
                            )
    )
    var_quitadas <- as.character(get(paste0('errores_train_var_', i))[get(paste0('errores_train_var_', i))$error==Inf, 'var'])
    num_inf <- length(var_quitadas)
    print(pasos_stepwise)
  }
  
  datos_train_fija <- get(paste0('errores_', i, '_', copula_stepwise))[[1]]
  datos_test_fija <- get(paste0('errores_', i, '_', copula_stepwise))[[2]]
  
  
  
  errores_train <- rbind(errores_train,
                         data.frame(iter = i,
                                    error = get(paste0('errores_train_var_', i))[ copula_stepwise, 'error'],
                                    var = variables[copula_stepwise]
                         )
  )
  errores_test <- rbind(errores_test,
                        data.frame(iter = i,
                                   error = get(paste0('errores_test_var_', i))[ copula_stepwise, 'error'],
                                   var = variables[copula_stepwise]
                        )
  )
  
  if (i == num_iter){
    errores_train <- data.frame()
    errores_test <- data.frame()
    mejor_num_iteraciones <- pasos_stepwise[which.min(pasos_stepwise$error_test), 'iteracion']
    for (j in 1:mejor_num_iteraciones){
      if (j == 1){
        errores_train <- rbind(errores_train,
                               data.frame(iter = j,
                                          error = pasos_stepwise[1, 'error_train'],
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
    
    print(pasos_stepwise)
    print(errores_train)
    print(errores_test)
  }
  
  i = i + 1
  iteracion <- iteracion + 1
  
  print(errores_train)
  print(errores_test)
}

iteraciones <- errores_test$iter
MAPE <- errores_test$error

plot(iteraciones, 
     MAPE, 
     type = "l")

write_xlsx( pasos_stepwise , 'pasos_stepwise_communities.xlsx')
write_xlsx( errores_train , 'errores_train_communities.xlsx')
write_xlsx( errores_test , 'errores_test_communities.xlsx')
