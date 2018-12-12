
KernelSmoothing.cdf <- function(xx,c0,bw)
  {
    # https://github.com/cran/DiagTest3Grp/blob/master/R/KernelSmoothing.cdf.R
    new.xx <- (c0-xx)/bw
    ks.prob <- pnorm(new.xx,mean=0,sd=1,lower.tail=TRUE)
    mean(ks.prob,na.rm=TRUE)
  }



puntuacion_copula_opt <- function(datos_iter, 
                                  n.ventas, 
                                  copulaoptima, 
                                  train){
  #ini <- Sys.time()
  d <- ncol(train)
  ventas <- vector.datos(train[, d], n.ptos = n.ventas)
  
  #Barrido de train simulados
  ptos.sim <- expand.grid(datos_iter[, 1:(d - 1)],
                          ventas)
  
  #C?lculo de CDF para los train simulados
  datos.orig <- list()
  
  grid.cdf.datos <- matrix(0, nrow = nrow(ptos.sim), ncol = d)
  grid.datos.orig <- matrix(0, nrow = nrow(ptos.sim), ncol = d)
   
  for (i in 1:d){
    
    uni_ptos_sim <- as.matrix(unique(ptos.sim[, i]))
    
    cdf <- apply(uni_ptos_sim,
                 1,
                 function(x){KernelSmoothing.cdf(x,
                                                 xx = train[, i],
                                                 bw = 0.1)}
    )
    
    cdf <- data.frame(uni_ptos_sim,
                      cdf = cdf)
    
    colnames(cdf)[1] <- colnames(ptos.sim)[i]
    
    ptos.sim_aux <- ptos.sim %>% left_join(cdf) 
  
    grid.cdf.datos[, i] <- ptos.sim_aux[,3]
    grid.datos.orig [, i] <- ptos.sim[,i]
  
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
  for (z in 1:(d - 1)){
    col_aux <- grid.cdf.datos[,z]
    col_aux2 <- ifelse(col_aux<0.99999,col_aux,0.99999)
    col_aux3 <- ifelse(col_aux2>0.00001,col_aux2,0.00001)
    grid.cdf.datos[,z] <- col_aux3
  }
   
  grid.cdf.datos <- data.table(grid.cdf.datos)
  var.indep <- colnames(train)[1:(ncol(train) - 1)]
  names(grid.cdf.datos)[1:(d - 1)] <- var.indep
  names(grid.cdf.datos)[d] <- 'y'
  
  grid.datos.orig <- data.table(grid.datos.orig)
  names(grid.datos.orig) <- paste(colnames(train), '.orig', sep = '')
  
  names(grid.datos.orig)[1:(d-1)] <- paste(var.indep, '.orig', sep = '')
  names(grid.datos.orig)[d] <- 'y.orig'
  
  #uniq.cdf.datos <- unique(grid.cdf.datos)

  distr.cop <- copula::dCopula(u = as.matrix(grid.cdf.datos), copula = copulaoptima)
  #print(Sys.time() - ini)
  #uniq.cdf.datos <- cbind(uniq.cdf.datos, distr.cop)
  #grid.cdf.datos <- grid.cdf.datos %>% left_join(uniq.cdf.datos)
  grid.cdf.datos <- cbind(grid.cdf.datos, distr.cop)
  grid.cdf.datos <- data.table(grid.cdf.datos)
  
  var.agrup <- paste(var.indep, collapse = ',')
  grid.cdf.datos <- grid.cdf.datos[, prob.margin := estim.area(y, distr.cop), by = var.agrup]
  
  grid.cdf.datos$distr.condic <- grid.cdf.datos$distr.cop/grid.cdf.datos$prob.margin
  
  grid.cdf.datos$esperanza <- grid.cdf.datos$distr.condic*grid.cdf.datos[, y]
  
  grid.cdf.datos <- grid.cdf.datos[, esp.condic := estim.area(y, esperanza), by = var.agrup]
  
  # grid.cdf.datos$esperanza.ord2 <- grid.cdf.datos$distr.condic*(grid.cdf.datos[, y]^2)
  
  # grid.cdf.datos <- grid.cdf.datos[, momen.ord2 := estim.area(y, esperanza.ord2), by = var.agrup]
  
  # grid.cdf.datos$std.condic <- sqrt(grid.cdf.datos$momen.ord2-(grid.cdf.datos$esp.condic)^2)
  
  # grid.esp.condic <- grid.cdf.datos[, c(var.indep, 'esp.condic','std.condic'), with = FALSE]
  
  grid.esp.condic <- grid.cdf.datos[, c(var.indep, 'esp.condic'), with = FALSE]
  
  #grid.esp.condic <- unique(grid.cdf.datos[, c(var.indep, 'esp.condic'), with = FALSE])
  
  grid.esp.condic <- cbind(grid.esp.condic, grid.datos.orig) 
  
  grid.esp.condic <- unique(grid.esp.condic[,1:(ncol(grid.esp.condic)-1)])
  
  y2 <- unique(ptos.sim[,d])
  dato_real <- unique(grid.cdf.datos[, d, with = FALSE])
  
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
  
  #var.orig <- paste(var.indep, '.orig', sep = '')
  nombres <- c(paste(names(train)[-d], '.orig', sep = ''), 'estim.copula', 'int.inferior', 'int.superior')
  
  final <- grid.esp.condic[, nombres, with = FALSE]
  names(final) <- c(names(train),'int.inferior', 'int.superior')
  #print(Sys.time() - ini)
  return(final)
}

