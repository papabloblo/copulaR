comprobacion_escalar <- function(x){
  if (length(x) > 1) {
    return('must be of length 1')
  } else if (!is.numeric(x) || x <= 0) {
    return('must be a positive number')
  } else {
    return(floor(x))
  }
}
