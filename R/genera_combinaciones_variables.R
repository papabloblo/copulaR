#' Title
#'
#' @param num_variables 
#' @param dim_copulas 
#'
#' @return
#' @export
#'
#' @examples
genera_combinaciones_variables <- function(num_variables, dim_copulas){
  
  n <- num_variables
  m <- dim_copulas - 1
  
  vector <- 1:n
  for (i in 1:m){
    aux <- list()
    if (i == 1){
      for (j in 1:n){
        aux[[j]] <- vector[j]
      }
      assign(paste('comb_', i, sep = ''), aux)
    } else {
      datos <- get(paste('comb_', (i-1), sep = ''))
      for (j in 1:length(datos)){
        if (i == 2){
          a <- datos[[j]]
        } else {
          a <- datos[[j]][length(datos[[j]])]
        }
        if (a < n){
          for (k in (a+1):n){
            aux[[length(aux) + 1]] <- c(datos[[j]], k)
          }
        }
      }
      assign(paste('comb_', i, sep = ''), aux)
    }
  }
  comb <- list()
  for (i in 1:m){
    a <- get(paste('comb_', i, sep = ''))
    comb <- c(comb, a)
  }
  return(comb)
}

