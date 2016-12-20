####################
#### Univariate ridge regression
####################

### usage: ridge.reg(target, dataset, lambda, B = 1, newdata = NULL) 


ridge.reg <- function(target, dataset, lambda, B = 1, newdata = NULL) {
  ## target is the dependent variable and can be a matrix as well
  ## However we only use it with a univariate target variable
  ## dataset contains the independent, CONTINUOUS ONLY, variables
  ## lambda is the ridge regularization parameter
  ## if lambda=0, the classical multivariate regression is implemented
  ## B is for bootstrap estimation of the standard errors of the betas
  ## newdata is the new independent variables values 
  ## whose values of y you want to estimate
  ## by default newdata is NULL
  
  target <- as.vector(target)
  dataset <- as.matrix(dataset)
  n <- length(target)  ## sample size
  p <- ncol(dataset)  ## dimensionality of dataset
  my <- sum(target) / n
  yy <- target - my  ## center the dependent variables
  s <- Rfast::colVars(dataset, std = TRUE) 
  mx <- as.vector( Rfast::colmeans(dataset) )
  xx <- t( ( t(dataset) - mx ) / s )

  xtx <- crossprod(xx)
  lamip <- lambda * diag(p)
  
  W <- solve( xtx + lamip )
  betas <- W %*% crossprod(xx, yy) 
  est <- xx %*% betas + my
  va <- var(target - est) * (n - 1) / (n - p - 1)
  # vab <- kronecker(va, W %*% xtx %*% W  ) 
  # seb <- as.vector( sqrt( diag(vab) ) )
  vab <- va * mahalanobis(W, numeric(p), xtx, inverted = TRUE)
  seb <- sqrt( vab )

  if (B > 1) { ## bootstrap estimation of the standard errors
    be <- matrix(nrow = B, ncol = p )
    for ( i in 1:B) {
       id <- sample(1:n, n, replace = TRUE)
       yb <- yy[id]     ;     xb <- xx[id, ]
       be[i, ] <- solve( crossprod(xb) + lamip, crossprod(xb, yb) )
    }
    seb <- Rfast::colVars(be, std = TRUE) ## bootstrap standard errors of betas
  } 
  
  be <- as.vector(betas)
  ## seb contains the standard errors of the coefficients
  
  if ( is.null( colnames(dataset) ) ) {
    names(seb) <- paste("X", 1:p, sep = "")
    names(be) <- paste("X", 1:p, sep = "")
    
  } else  names(seb) <- names(be) <- colnames(dataset)
  if ( !is.null(newdata) ) {
    
    newdata <- as.matrix(newdata)
    newdata <- matrix(newdata, ncol = p)
    newdata <- t( ( t(newdata) - mx ) / s ) ## standardize the independent variables of the new data
    est <- newdata %*% beta + my 
  } else est <- est
  
  est <- as.vector(est) 
  list(beta = be, seb = seb, est = est)
  
}



### References
### Hoerl, A.E. and R.W. Kennard (1970). Ridge regression: Biased estimation for nonorthogonal problems". Technometrics, 12(1):55?67
### Brown, P. J. (1994). Measurement, Regression and Calibration. Oxford Science Publications.

