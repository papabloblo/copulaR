puntuacion_copula_opt <- function(datos_iter, n.ventas, copulaoptima, train){
  
  #ini <- Sys.time()
  d <- ncol(train)
  
  Ventas <-  vector.datos(train[,d], n.ptos = n.ventas)
  #Barrido de train simulados
  names_hist <- paste0(colnames(train)[1:(d-1)], '_hist')
  parte1 <- do.call(rbind, replicate(n.ventas, as.matrix(datos_iter[,names_hist]),
                                     simplify=FALSE))
  rownames(parte1) <- NULL
  parte2 <- data.frame(Ventas = sort(unlist(replicate(nrow(data.frame(datos_iter[,names_hist])),
                                                      as.matrix(Ventas), 
                                                      simplify=FALSE))))
  rownames(parte2) <- NULL
  ptos.sim <- data.frame(parte1,
                         parte2)
  
  #ptos.sim <- data.frame(ptos.sim)
  
  #C?lculo de CDF para los train simulados
  datos.orig <- list()
  
  grid.cdf.datos <- matrix(0,nrow = nrow(ptos.sim), ncol = d)
  grid.datos.orig <- matrix(0,nrow = nrow(ptos.sim), ncol = d)
  
  for (i in 1:d){
    
    uni_ptos_sim <- as.matrix(ptos.sim[!duplicated(ptos.sim[,i]), i])
    
    # cdf <- apply(uni_ptos_sim,
    #              1,
    #              function(x){
    #                predict(ks::kcde(train[,i]), x = x)
    #              }
    # )
    
    cdf <- predict(ks::kcde(train[,i]), x = uni_ptos_sim)
    
    cdf <- data.frame(uni_ptos_sim,
                      cdf = cdf)
    if (i == d){
      cdf_aux <- cdf
      cdf_aux$redondeo <- round(cdf_aux$cdf, 7)
      tabla_agrup <- cdf_aux %>% 
        group_by(redondeo) %>%
        summarise(media = mean(uni_ptos_sim))
      y2 <- tabla_agrup$media
    }
    
    colnames(cdf)[1] <- colnames(ptos.sim)[i]
    
    ptos.sim_aux <- ptos.sim %>% left_join(cdf, by = colnames(ptos.sim)[i]) 
    
    grid.cdf.datos[, i] <- ptos.sim_aux[,d + 1]
    grid.datos.orig[, i] <- ptos.sim[,i]
    
  }
  
  #Se crea un grid con todas las posibles combinaciones de datos
  
  
  #grid.cdf.datos <- as.data.table(grid.cdf.datos)
  #names(grid.cdf.datos) <- names(train)
  
  #var.indep <- paste('x', 1:(d-1), sep ='')
  
  #Se crea un grid con todas las posibles combinaciones de datos
  
  #grid.cdf.datos <- expand.grid(cdf.datos)
  
  #grid.datos.orig <- as.data.table(grid.datos.orig)
  #names(grid.datos.orig) <- paste(names(train), '.orig', sep = '')
  
  #names(grid.datos.orig)[1:(d-1)] <- paste(var.indep, '.orig', sep = '')
  #names(grid.datos.orig)[d] <- 'y.orig'
  #grid.cdf.datos <- data.frame(grid.cdf.datos)
  # for (z in 1:d){
  #   col_aux <- grid.cdf.datos[,z]
  #   col_aux2 <- ifelse(col_aux<0.999,col_aux,0.999)
  #   col_aux3 <- ifelse(col_aux2>0.001,col_aux2,0.001)
  #   grid.cdf.datos[,z] <- col_aux3
  # }
  #  
  
  grid.cdf.datos <- ifelse(grid.cdf.datos<0.99999,grid.cdf.datos,0.99999)
  grid.cdf.datos <- ifelse(grid.cdf.datos>0.00001,grid.cdf.datos,0.00001)
  
  grid.cdf.datos <- data.table(grid.cdf.datos)
  var.indep <- names(train)[1:(ncol(train)-1)]
  names(grid.cdf.datos)[1:(d-1)] <- var.indep
  names(grid.cdf.datos)[d] <- 'y'
  
  grid.datos.orig <- data.table(grid.datos.orig)
  names(grid.datos.orig) <- paste(names(train), '.orig', sep = '')
  
  names(grid.datos.orig)[1:(d-1)] <- paste(var.indep, '.orig', sep = '')
  names(grid.datos.orig)[d] <- 'y.orig'
  
  grid.cdf.datos <- round(grid.cdf.datos, 7)
  cruce_final <- cbind(grid.datos.orig[,1:(d-1)], grid.cdf.datos[,1:(d-1)])
  cruce_final <- cruce_final[!duplicated(cruce_final),]
  ind_dup <- duplicated(grid.cdf.datos)
  grid.cdf.datos <- grid.cdf.datos[!ind_dup, ]
  grid.datos.orig <- grid.datos.orig[!ind_dup, ]
  #uniq.cdf.datos <- unique(grid.cdf.datos)
  
  # grid.cdf.datos$aux <- cCopula(as.matrix(grid.cdf.datos), copula = copulaoptima@copula, indices = d, inverse = T)
  # 
  # grid.esp.condic <- grid.cdf.datos %>%
  #   group_by_at(vars(one_of(var.indep))) %>%
  #   summarise(esp.condic = mean(aux))
  
  #distr.cop <- dCopula(as.matrix(grid.cdf.datos), copulaoptima@copula)
  distr.cop <- BiCopPDF(unlist(grid.cdf.datos[,1]), unlist(grid.cdf.datos[,2]), copulaoptima)
  # grid.cdf.datos_aux <- data.frame(grid.cdf.datos_aux,
  #                                  distr.cop = distr.cop)
  # grid.cdf.datos <- grid.cdf.datos %>% left_join(grid.cdf.datos_aux)
  
  #print(Sys.time() - ini)
  #uniq.cdf.datos <- cbind(uniq.cdf.datos, distr.cop)
  #grid.cdf.datos <- grid.cdf.datos %>% left_join(uniq.cdf.datos)
  grid.cdf.datos <- cbind(grid.cdf.datos, distr.cop)
  grid.cdf.datos <- data.table(grid.cdf.datos)
  
  var.agrup <- paste(var.indep, collapse = ',')
  grid.cdf.datos <- grid.cdf.datos[, prob.margin := estim.area(y, distr.cop), by = var.agrup]
  
  # grid.cdf.datos <- grid.cdf.datos %>%
  #   group_by_at(vars(one_of(var.indep))) %>%
  #   mutate( prob.margin = sum(abs(diff(y)*rollmean(distr.cop, 2))))
  
  
  grid.cdf.datos$distr.condic <- grid.cdf.datos$distr.cop/grid.cdf.datos$prob.margin
  
  
  grid.cdf.datos$esperanza <- grid.cdf.datos$distr.condic*grid.cdf.datos$y
  
  
  grid.cdf.datos <- grid.cdf.datos[, esp.condic := estim.area(y, esperanza), by = var.agrup]
  
  # grid.cdf.datos <- grid.cdf.datos %>%
  #   group_by_at(vars(one_of(var.indep))) %>%
  #   mutate( esp.condic = sum(abs(diff(y)*rollmean(esperanza, 2))))
  # grid.cdf.datos$esperanza.ord2 <- grid.cdf.datos$distr.condic*(grid.cdf.datos[, y]^2)
  
  # grid.cdf.datos <- grid.cdf.datos[, momen.ord2 := estim.area(y, esperanza.ord2), by = var.agrup]
  
  # grid.cdf.datos$std.condic <- sqrt(grid.cdf.datos$momen.ord2-(grid.cdf.datos$esp.condic)^2)
  
  # grid.esp.condic <- grid.cdf.datos[, c(var.indep, 'esp.condic','std.condic'), with = FALSE]
  
  grid.esp.condic <- grid.cdf.datos[, c(var.indep, 'esp.condic'), with = FALSE]
  
  #grid.esp.condic <- unique(grid.cdf.datos[, c(var.indep, 'esp.condic'), with = FALSE])
  
  grid.esp.condic <- cbind(grid.esp.condic, grid.datos.orig)
  
  grid.esp.condic <- grid.esp.condic[!duplicated(grid.esp.condic[,1:(ncol(grid.esp.condic)-1)]),1:(ncol(grid.esp.condic)-1)]
  
  dato_real <- grid.cdf.datos[!duplicated(grid.cdf.datos[,d, with = FALSE]), d, with = FALSE]
  
  grid.esp.condic$estim.copula <- sapply(grid.esp.condic$esp.condic,
                                         function(x){y2[min(which(x < dato_real))]})
  
  # banda.superior <- grid.esp.condic$esp.condic+grid.esp.condic$std.condic
  # banda.superior[banda.superior>max(cdf.datos[[d]])] <- max(cdf.datos[[d]])
  # 
  # banda.inferior <- grid.esp.condic$esp.condic-grid.esp.condic$std.condic
  # banda.inferior[banda.inferior<min(cdf.datos[[d]])] <- min(cdf.datos[[d]])
  # 
  # grid.esp.condic$int.inferior <- sapply(banda.inferior,
  #                                        function(x){y2[min(which(x <= cdf.datos[[d]]))]})
  # 
  # grid.esp.condic$int.superior <- sapply(banda.superior,
  #                                        function(x){y2[min(which(x <= cdf.datos[[d]]))]})
  
  grid.esp.condic$int.inferior <- 0
  grid.esp.condic$int.superior <- 0
  
  cruce_final <- cruce_final %>% left_join(grid.esp.condic[,colnames(grid.esp.condic)[!grepl('.orig', colnames(grid.esp.condic))], with = FALSE], by = colnames(datos_iter)[1:(d-1)])
  #var.orig <- paste(var.indep, '.orig', sep = '')
  nombres <- c(paste(names(train)[-d], '.orig', sep = ''), 'estim.copula', 'int.inferior', 'int.superior')
  
  final <- cruce_final[, nombres]
  colnames(final)[1:(d-1)] <- paste0(colnames(datos_iter)[1:(d-1)], '_hist')
  final <- final  %>% left_join(datos_iter, by = names_hist)
  final <- final[,colnames(final)[!colnames(final) %in% names_hist]]
  colnames(final)[1] <- "ERROR"
  #print(Sys.time() - ini)
  return(final)
}

