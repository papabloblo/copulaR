set.seed(1234)

ker.cdf <- function(ptos){sapply(ptos,
                                 KernelSmoothing.cdf,
                                 xx=ptos, 
                                 bw=0.1)}

vector.datos <- function(x, n.ptos = n2){
  seq(min(x), max(x), length.out = n.ptos)
}

estim.area <- function(dis, par2){
  sum(abs(diff(dis)*rollmean(par2, 2)))
}

copula.optima <- function(puntos){
      
    d <- ncol(puntos)
    cdf.ptos <- apply(puntos, 2, ker.cdf) 
    test_independencia <- chisq.test(cdf.ptos[,1],
                                     cdf.ptos[,2],
                                     correct =FALSE)$p.value
    if (test_independencia<=0.05){
      #print(head(cdf.ptos))
      copulas <- c('Clayton',
                   'Gumbel',
                   'Frank',
                   'Normal',
                   'Joe',
                   't')
      #Definición de las cópulas a comprobar
      
      cop1 <- claytonCopula( ,dim = d)
      cop2 <- gumbelCopula(  ,dim = d)
      cop3 <- frankCopula(   ,dim = d)
      cop4 <- normalCopula(  ,dim = d)
      cop5 <- joeCopula( ,dim = d)
      
      cop6 <- tCopula(       ,dim = d, dispstr='un')
      
      #Estimación de las cópulas
      fit1 <- fitCopula(cop1, cdf.ptos, method = "itau")
      fit2 <- fitCopula(cop2, cdf.ptos, method = "itau")
      fit3 <- fitCopula(cop3, cdf.ptos, method = "itau")
      fit4 <- fitCopula(cop4, cdf.ptos, method = "ml")
      
      # print(5)
      fit5 <- fitCopula(cop5, cdf.ptos, method = "mpl",
                        start = 1.001,
                        lower = 1.001,
                        upper = Inf)
      
      fit6 <- fitCopula(cop6, cdf.ptos, method = "itau.mpl")
      
      fit <- c(fit1, fit2, fit3, fit4,
               fit5, 
               fit6
      )
      # print(fit1)
      # print(fit2)
      # print(fit3)
      # print(fit4)
      # print(fit5)
      # print(fit6)
      #Cálculo del logaritmo de la m?xima verosimilitud (log-likelihood)
      log.lik1 <- loglikCopula(fit1@copula@parameters,
                               cdf.ptos,
                               cop1)
      
      log.lik2 <- loglikCopula(fit2@copula@parameters,
                               cdf.ptos,
                               cop2)
      
      log.lik3 <- loglikCopula(fit3@copula@parameters,
                               cdf.ptos,
                               cop3)
      
      log.lik4 <- loglikCopula(fit4@copula@parameters,
                               cdf.ptos,
                               cop4)
      
      log.lik5 <- loglikCopula(fit5@copula@parameters,
                               cdf.ptos,
                               cop5)
      
      log.lik6 <- loglikCopula(fit6@copula@parameters,
                               cdf.ptos,
                               cop6)
      
      # log.lik6 <- c()
      # for (i in 1:100){
      #   log.lik6[i] <- loglikCopula(c(fit6@copula@parameters, i),
      #                               cdf.ptos,
      #                               cop6)
      # }
      
      # cop6 <- tCopula(fit6@copula@parameters[1],
      #                 df = which.max(log.lik6))
      # log.lik6 <- log.lik6[which.max(log.lik6)]
      
      log.lik <- c(log.lik1, log.lik2, log.lik3, log.lik4,
                   log.lik5, 
                   log.lik6
      )
      
      if (d == 2) {
        copulas <- c(copulas,
                     'husler', 
                     'tawn')
        ## Estas c?pulas solo est?n implementadas para dimensi?n 2 ##
        # cop7 <- amhCopula(     ,dim = d)
        
        ev <- c('huslerReiss', 'tawn')
        
        cop8 <- evCopula(ev[1],  dim = d)
        cop9 <- evCopula(ev[2],  dim = d)
        ##
        
        # fit7 <- fitCopula(cop7, cdf.ptos, method = "itau")
        
        fit8 <- fitCopula(cop8, cdf.ptos, method = "ml",
                          start = 1.001,
                          lower = 0,
                          upper = 1,
                          optim.method="Brent")
        
        fit9 <- fitCopula(cop9, cdf.ptos, method = "mpl",
                          start = 1.001,
                          lower = 0,
                          upper = 1,
                          optim.method="Brent")
        
        # print(fit8)
        # print(fit9)
        fit <- c(fit, #fit7,
                 fit8, fit9)
        
        
        # log.lik7 <- loglikCopula(fit7@copula@parameters,
        #                          cdf.ptos,
        #                          cop7)
        
        log.lik8 <- loglikCopula(fit8@copula@parameters,
                                 cdf.ptos,
                                 cop8)
        
        log.lik9 <- loglikCopula(fit9@copula@parameters,
                                 cdf.ptos,
                                 cop9)
        
        log.lik <- c(log.lik,
                     #log.lik7,
                     log.lik8, log.lik9)
      }
      
      log.lik[(log.lik =="NaN") | (log.lik==(Inf))] <- -Inf
      
      #C?lculo AIC
      aic <- 2 - (2*log.lik)
      
      
      
      #Selecci?n del ajuste con menor AIC
      mejorcopula <- which.min(aic)
      
      copulaoptima <- fit[[mejorcopula]]
      
      copulas <- data.frame(copula = copulas,
                            aic = aic)
      
      copulas <- copulas[order(copulas$aic),]
      
      copula_final <- list(aic = copulas,
                           copulaoptima = copulaoptima,
                           indep = test_independencia)
    } else {
      copula_final <- list(aic = '',
                           copulaoptima = '',
                           indep = test_independencia)
    }

  
  return(copula_final)
}

