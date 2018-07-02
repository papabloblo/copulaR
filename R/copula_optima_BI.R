set.seed(1234)

ker.cdf <- function(x){
  predict(ks::kcde(x), x = as.matrix(x))
}

vector.datos <- function(x, n.ptos = n2){
  seq(min(x), max(x), length.out = n.ptos)
}

estim.area <- function(dis, par2){
  sum(abs(diff(dis)*caTools::runmean(par2, 3)[1:(length(dis)-1)]))
}

copula.optima <- function(puntos, num_obs_fit){
      
    d <- ncol(puntos)
    cdf.ptos <- apply(puntos, 2, ker.cdf)
    if (!is.null(num_obs_fit)){
      ind_muestra <- sample(1:nrow(cdf.ptos), num_obs_fit, replace = TRUE)
      cdf.ptos.aux <- cdf.ptos[ind_muestra,]
      if (sum(apply(as.matrix(cdf.ptos.aux[,-d]),
                    2,
                    function(x){length(unique(x))}) == 1) == 0){
        cdf.ptos <- cdf.ptos.aux
      }
    }
    
    familias <- c(0,
                  1,
                  2,
                  3,
                  4,
                  5,
                  6,
                  7,
                  8,
                  9,
                  10,
                  13,
                  14,
                  16,
                  17,
                  18,
                  19,
                  20,
                  23,
                  24,
                  26,
                  27,
                  28,
                  29,
                  30,
                  33,
                  34,
                  36,
                  37,
                  38,
                  39,
                  40,
                  104,
                  114,
                  124,
                  134,
                  204,
                  214,
                  224,
                  234)
    familias_itau <- c(1,2,3,4,5,6,13,14,16,23,24,26,33,34,36)
    fit <- list()
    for (i in 1:length(familias)){
      if (familias[i] %in% familias_itau){
        capture.output(try(fit[[i]] <- BiCopEst(cdf.ptos[,1], cdf.ptos[,2], familias[i], method = "itau", 
                                                max.df = 3000,
                                                max.BB = list(BB1=c(500,600),BB6=c(600,600),BB7=c(500,600),BB8=c(600,1))),
                           silent = TRUE))
      } else {
        capture.output(try(fit[[i]] <- BiCopEst(cdf.ptos[,1], cdf.ptos[,2], familias[i], method = "mle", 
                                                max.df = 3000,
                                                max.BB = list(BB1=c(500,600),BB6=c(600,600),BB7=c(500,600),BB8=c(600,1))),
                           silent = TRUE))
      }
      
      if (length(fit)<i){
        fit[[i]] <- BiCopEst(cdf.ptos[,1], cdf.ptos[,2], familias[1], method = "itau")
      }
    }
    
    #Selecci?n del ajuste con menor AIC
    aic <- unlist(sapply(fit,function(x){x$AIC}))
    mejorcopula <- which.min(aic)
    copulas <- unlist(sapply(fit,function(x){x$familyname}))
    
    copulaoptima <- fit[[mejorcopula]]
    
    copulas <- data.frame(copula = copulas,
                          aic = aic)
    
    copulas <- copulas[order(copulas$aic),]
    
    copula_final <- list(aic = copulas,
                         copulaoptima = copulaoptima,
                         indep = min(copulas$aic)>=0)

  return(copula_final)
}

