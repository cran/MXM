vara <- function( x ) {
  n <- length(x)
  sum(x^2) - sum(x)^2 / n 
}
