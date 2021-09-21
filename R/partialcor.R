partialcor <- function(R, indx, indy, indz, n) {
  ## R is a correlation matrix
  ## i and j denote the two variables whose conditional correlation is to be estimated 
  ## k denotes the set of conditioning variables
  if ( indz == 0 ) { ## classical correlation
    r <- R[indx, indy]   
  } else {
    rho <- try( solve( R[c(indx, indy, indz), c(indx, indy, indz)] ), silent = TRUE)
    if ( !identical( class(rho), "try-error" ) ) {
	    r <-  - rho[1, 2] / sqrt(rho[1, 1] * rho[2, 2])
  	} else r <- 0.99999
  }
  if ( abs(r) >= 1 | is.na(r) )  r <- 0.99999
  z <- 0.5 * log( (1 + r) / (1 - r) ) * sqrt( n - sum(indz > 0) - 3 )
  pval <- log(2) + pt( abs(z), n - sum(indz > 0) - 3, lower.tail = FALSE, log.p = TRUE )
  res <- c(r, pval)
  names(res) <- c( "partial correlation", "logged p-value")
  res
}
