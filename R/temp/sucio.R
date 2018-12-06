

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
