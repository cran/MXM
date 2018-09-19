################################
#### Projected multivariate normal for circular or angular regression
#### Tsagris Michail 1/2014
#### mtsagris@yahoo.gr
#### Presnell, Morrison and Littell (1998), JASA
################################
spml.reg2 <- function(y, x, tol = 1e-07, seb = FALSE, maxiters = 100) {
  ## y is the angular dependent variable
  ## x contains the independent variable(s)
  x <- model.matrix(~., data.frame(x) )
  if (is.matrix(y)) {
    u <- y
    ci <- u[, 1]
    si <- u[, 2]
  } else {
    ci <- cos(y)
    si <- sin(y)
    u <- cbind(ci, si)
  }
  n <- dim(u)[1]
  XX <- solve( crossprod(x), t(x) )
  be <-  XX %*% u 
  mu <- x %*% be
  f <-  - 0.5   ;   con <- sqrt(2 * pi) 
  tau <- rowsums(u * mu)
  ptau <- pnorm(tau)
  lik1 <-  - 0.5 * sum( mu^2 ) + sum( log1p( tau * ptau * con / exp(f * tau^2) ) )
  psit <- tau + ptau / ( exp(f * tau^2)/con + tau * ptau )
  be <- XX %*% (psit * u)
  mu <- mu * psit        
  tau <- rowsums(u * mu)
  ptau <- pnorm(tau)
  lik2 <-  - 0.5 * sum( mu^2 ) + sum( log1p( tau * ptau * con / exp(f * tau^2) ) )
  i <- 2
  ## mono th while
  while ( lik2 - lik1 > tol  & i < maxiters ) {
    lik1 <- lik2
    i <- i + 1
    rat <- ptau / ( exp(f * tau^2)/con + tau * ptau )
    psit <- tau + rat    
    psit2 <- 2 - tau * rat - rat^2
    der <- as.vector( crossprod(x, - mu + psit * u) )
    a11 <- crossprod(x, x * (psit2 * ci^2 - 1) )
    a12 <- crossprod(x, x * (psit2 * ci * si ) )
    a22 <- crossprod(x, x * (psit2 * si^2 - 1 ) )
    der2 <- cbind( rbind(a11, a12), rbind(a12, a22) )
    be <- be - solve(der2, der)
    mu <- x %*% be
    tau <- rowsums(u * mu)
    ptau <- pnorm(tau)
    lik2 <-  - 0.5 * sum( mu^2 ) + sum( log1p( tau * ptau * con / exp(f * tau^2) ) )
  }
  ###
  if ( seb ) {
    seb <-  sqrt( diag( solve( -der2 ) ) )
    seb <- matrix(seb, ncol = 2)
    colnames(seb) <- c("Cosinus of y", "Sinus of y")
    rownames(seb) <- colnames(x)    
  } else   seb <- NULL
  loglik <- lik2 - n * log(2 * pi)
  colnames(be) <- c("Cosinus of y", "Sinus of y")
  rownames(be) <- colnames(x)
  list(iters = i, loglik = loglik, be = be, seb = seb)  
}