big.model <- function(y, x, test) {
  
  if ( test == "testIndLogistic" ) {
    mod <- Rfast::glm_logistic(x, y, maxiters = 5000)
    res <- list(be = mod$be, devi = mod$devi, phi = NA)
    
  } else if ( test == "testIndMultinom" ) {
    mod <- try( Rfast::multinom.reg(x, y, maxiters = 10000), silent = TRUE )
    if ( identical(class, "try-error") ) {
      res <- list(be = NULL, devi = NULL, phi = NA)
    } else   res <- list(be = mod$be, devi = mod$devi, phi = NULL)
    
  } else if ( test == "testIndPois" ) {
    mod <- Rfast::glm_poisson(x, y)
    res <- list(be = mod$be, devi = mod$devi, phi = NA)
    
  } else if ( test == "testIndQPois" ) {
    mod <- Rfast::qpois.reg(x, y, maxiters = 5000)
    res <- list(be = mod$be, devi = mod$devi, phi = mod$phi)
    
  } else if ( test == "censIndWR" ) {
    mod <- survival::survreg(y ~ x, control = list(iter.max = 10000) )
    res <- list(be = mod$coefficients, devi = 2 * mod$loglik[2], phi = NA)
    
  } else if ( test == "testIndFisher" ) {
    mod <- try( Rfast::lmfit( cbind(1, x), y ), silent = TRUE )
    if ( identical(class, "try-error") ) {
      res <- list(be = NULL, devi = NULL, phi = NA)
    } else {
      dm <- dim(x)
      devi <- sum( (mod$residuals)^2 )
      res <- list(be = mod$be, devi = devi, phi = devi/(dm[1] - dm[2] - 1) )
    }  
  }
  
  res
  
}
