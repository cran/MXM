#######################
### K-fold cross validation for selecting the optimal lambda in ridge regression
#######################
 
### usage:  ridgereg.cv( target, dataset, K = 10, lambda = seq(0, 1, by = 0.01), 
### auto = FALSE, seed = FALSE, ncores = 1, mat = NULL )

ridgereg.cv <- function( target, dataset, K = 10, lambda = seq(0, 2, by = 0.1), 
                        auto = FALSE, seed = FALSE, ncores = 1, mat = NULL ) {
  ## target is a dependent continuous variable or a matrix with continuous variables
  ## dataset is a matrix with continuous variables only
  ## K is the number of folds for the K-fold cross validation, set to 10 by default
  ## lambda is a vector containing a grid of values
  ## auto is a boolean variable. If TRUE, the GCV criterion is used to return the best lambda automatically.
  ## Otherwise, a K-fold cross validation is performed
  ## seed is boolean. Should the same K-folds be used always or not?
  ## ncores specifies how many cores to be used

  target <- as.vector(target)
  dataset <- as.matrix(dataset)
  n <- length(target)
  p <- ncol(dataset)
  mspe <- NULL
  performance <- NULL

  if ( auto == TRUE ) {
    runtime <- proc.time()
    mod <- MASS::lm.ridge( target ~ dataset, lambda = lambda ) 
    gcv <- min(mod$GCV)
    plot(lambda, mod$GCV, type = "b", xlab = expression(paste(lambda, " values")), ylab = "GCV")
    lam <- lambda[ which.min(mod$GCV) ]
    runtime <- proc.time() - runtime

  } else {
    if ( is.null(mat) ) { ## no folds were given by the user
      if (seed == TRUE)  set.seed(1234567) ## the folds will always be the same
      nu <- sample(1:n, min( n, round(n / K) * K ) )
      ## It may be the case this new nu is not exactly the same
      ## as the one specified by the user
      options(warn = -1)  # if the length of nu does not fit 
      ## to a matrix a warning message should appear 
      mat <- matrix( nu, ncol = K ) 
    } else mat <- mat
      rmat <- nrow(mat)
      
    if ( ncores == 1 ) {
      
      runtime <- proc.time()
      mi <- length(lambda)
      per <- matrix( nrow = K, ncol = mi )
      
      for (vim in 1:K) {
        ytest <- as.vector( target[ mat[, vim] ] )  ## test set dependent vars
        ytrain <- as.vector( target[ -mat[, vim] ] )  ## train set dependent vars
        my <- sum(ytrain) / (n - rmat)
        xtrain <- as.matrix( dataset[ -mat[, vim], ] )  ## train set independent vars
        xtest <- as.matrix( dataset[ mat[, vim], ] )  ## test set independent vars
        mx <- as.vector( Rfast::colmeans(xtrain) )
        yy <- ytrain - my  ## center the dependent variables
        s <- Rfast::colVars(xtrain, std = TRUE) 
        xtest <- ( t(xtest) - mx ) / s ## standardize the newdata values 
        xtest <- t(xtest)
        
        xx <- ( t(xtrain) - mx ) / s
        xx <- t(xx)
        sa <- svd(xx)
        tu <- t(sa$u)    ;    d <- sa$d    ;    v <- sa$v
        
        for (i in 1:mi) {
          betas <- ( v %*% (tu * ( d / ( d^2 + lambda[i] ) ) ) ) %*% yy 
          est <- xtest %*% betas + my 
            per[vim, i] <- sum( (ytest - est)^2 ) / rmat
        }
        
      }
      
      runtime <- proc.time() - runtime

    } else {
      
      runtime <- proc.time()
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      mi <- length(lambda)
      pe <- numeric(mi)
	  
      per <- foreach(vim = 1:K, .combine = rbind) %dopar% {
        ytest <- as.vector( target[mat[, vim] ] )  ## test set dependent vars
        ytrain <- as.vector( target[-mat[, vim] ] )  ## train set dependent vars
        my <- sum(ytrain) / (n - rmat)
        xtrain <- as.matrix( dataset[-mat[, vim], ] )  ## train set independent vars
        xtest <- as.matrix( dataset[mat[, vim], ] )  ## test set independent vars
        mx <- as.vector( Rfast::colmeans(xtrain) )
        yy <- ytrain - my  ## center the dependent variables
        s <- Rfast::colVars(xtrain, std = TRUE) 
        
        xx <- ( t(xtrain) - mx ) / s
        xx <- t(xx)
        xtest <- ( t(xtest) - mx ) / s ## standardize the newdata values   
        xtest <- t(xtest)

        sa <- svd(xx)
        tu <- t(sa$u)    ;    d <- sa$d    ;    v <- sa$v
        for ( i in 1:mi ) {
          betas <- ( v %*% ( tu *( d / ( d^2 + lambda[i] ) ) ) ) %*% yy 
          est <- xtest %*% betas + my 
          pe[i] <- sum( (ytest - est)^2 ) / rmat
        }
        return(pe)
      }
	  
      stopCluster(cl)
	  
      runtime <- proc.time() - runtime
    }

    per <- as.matrix(per)
    mspe <- as.vector( Rfast::colmeans(per) )
    bias <- per[ , which.min(mspe)] - apply(per, 1, min)  ## TT estimate of bias
    estb <- mean( bias )  ## TT estimate of bias
    names(mspe) <- lambda
    lam <- lambda[ which.min(mspe) ]
    
    plot(lambda, mspe, xlab = expression(paste(lambda, " values")), ylab = "MSPE", type = "b")
    names(mspe) <- lambda
    
    performance <- c( min(mspe) + estb, estb)
    names(performance) <- c("Estimated MSPE", "Estimated bias")
  }
  
  list(mspe = mspe, lambda = lam, performance = performance, runtime = runtime)
}


### References
### Hoerl, A.E. and R.W. Kennard (1970). Ridge regression: Biased estimation for nonorthogonal problems". Technometrics, 12(1):55?67.

### Brown, P. J. (1994). Measurement, Regression and Calibration. Oxford Science Publications.

### Tibshirani Ryan J., and Tibshirani Robert (2009). A bias correction for the minimum error rate in cross-validation.
#The Annals of Applied Statistics 3(2): 822-829.
