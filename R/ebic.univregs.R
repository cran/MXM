ebic.univregs <- function(target, dataset, targetID = -1, test = NULL, user_test = NULL, wei = NULL, ncores = 1, gam = NULL) {
  
  dm <- dim(dataset)
  rows <- dm[1]
  cols <- dm[2]
  if (targetID != -1 ) {
    target <- dataset[, targetID]
    dataset[, targetID] <- rbinom(rows, 1, 0.5)
  }   
  id <- NULL
  ina <- 1:cols
  id <- Rfast::check_data(dataset)
  if ( sum(id > 0) )  dataset[, id] <- rnorm(rows * length(id) )
  logn <- log(rows)
  
  if ( is.null(gam) ) {
    con <- 2 - log(cols) / log(rows)
  } else con <- 2 * gam
  if ( (con) < 0 )  con <- 0
  
  if ( !is.null(user_test) ) {
    ebic <- ebicScore(target, dataset, test = user_test, wei, targetID)

  } else if ( identical(test, testIndBeta) ) {  ## Beta regression
    ebic <- beta.regs(target, dataset, wei, logged = TRUE, ncores = ncores)[, 3]

  } else if ( identical(test, testIndMMReg) ) {  ## M (Robust) linear regression

    if ( ncores <= 1 | is.null(ncores) ) {
      ebic <- numeric(cols)
      for ( i in 1:cols ) {
        fit2 <- MASS::rlm(target ~ dataset[, i], maxit = 2000, method = "MM" )
        ebic[i] <- BIC(fit2) 
      } 

    } else {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      ebic <- foreach(i = 1:cols, .combine = rbind, .packages = "MASS") %dopar% {
        fit2 <- MASS::rlm(target ~ dataset[, i], maxit = 2000, method = "MM" )
        return( BIC(fit2) ) 
      }  
      stopCluster(cl)
    }   
    
  } else if ( identical(test, testIndReg) ) {  ## linear regression

    if ( ncores <= 1 | is.null(ncores) ) {
      ebic <- numeric(cols)
      for ( i in 1:cols ) {
        fit2 <- lm( target ~ dataset[, i], weights = wei, y = FALSE, model = FALSE )
        ebic[i] <- BIC(fit2)
      }

    } else {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      ebic <- foreach(i = 1:cols, .combine = rbind) %dopar% {
        ww <- lm( target ~ dataset[, i], weights = wei, y = FALSE, model = FALSE )
        return( BIC( ww ) )
      }
      stopCluster(cl)
    }   
    
  } else if ( identical(test, testIndOrdinal) ) {  ## ordinal regression

    if ( ncores <= 1 | is.null(ncores) ) {
      ebic <- numeric(cols)
      for ( i in 1:cols ) {
        fit2 <- ordinal::clm(target ~ dataset[, i], weights = wei)
        ebic[i] <- BIC(fit2)
      }

    } else {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      ebic <- foreach(i = 1:cols, .combine = rbind, .packages = "ordinal") %dopar% {
        fit2 <- ordinal::clm(target ~ dataset[, i], weights = wei)
        return( BIC(fit2) )
      } 
      stopCluster(cl)
    }
    
  } else if ( identical(test, testIndMultinom) ) {  ## multinomial regression

    if ( ncores <= 1 | is.null(ncores) ) {
      ebic <- numeric(cols)
      for ( i in 1:cols ) {
        fit2 <- nnet::multinom(target ~ dataset[, i], trace = FALSE, weights = wei )
        ebic[i] <- BIC(fit2)
      }
      
    } else {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      ebic <- foreach(i = 1:cols, .combine = rbind, .packages = "nnet") %dopar% {
        fit2 <- nnet::multinom(target ~ dataset[, i], weights = wei)
        return( BIC(fit2) )
      }
      stopCluster(cl)
    }
    
  } else if ( identical(test, testIndLogistic)  &  is.matrix(dataset)  &  is.null(wei) ) {  ## logistic regression
    if ( is.factor(target) )   target <- as.numeric(target) - 1
    ebic <- Rfast::logistic_only(dataset, target) + 2 * logn
    
  } else if ( identical(test, testIndLogistic) ) {  ## Logistic regression

    if ( ncores <= 1 | is.null(ncores) ) {
      ebic <- numeric(cols)
      for ( i in 1:cols ) {
        fit2 <- glm( target ~ dataset[, i], binomial, weights = wei )
        ebic[i] <- BIC(fit2)
      }
      
    } else {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      ebic <- foreach(i = 1:cols, .combine = rbind) %dopar% {
        fit2 <- glm( target ~ dataset[, i], binomial, weights = wei )
        return( BIC(fit2) )
      }
      stopCluster(cl)
    }
    
  } else if ( identical(test, testIndBinom) ) {  ## Binomial regression
    wei <- target[, 2] 
    y <- target[, 1] / wei

    if ( ncores <= 1 | is.null(ncores) ) {
      ebic <- numeric(cols)
      for ( i in 1:cols ) {
        fit2 = glm( y ~ dataset[, i], binomial, weights = wei )
        ebic[i] = BIC(fit2)       
      }
      
    } else {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      wei <- target[, 2] 
      y <- target[, 1] / wei
      ebic <- foreach(i = ina, .combine = rbind) %dopar% {
        fit2 <- glm( y ~ dataset[, i], binomial, weights = wei )
        return( BIC(fit2) )
      }
      stopCluster(cl)
    }
    
  } else if ( identical(test, testIndPois)  &  is.matrix(dataset)  &  is.null(wei) ) {  ## logistic regression
    if ( is.factor(target) )   target <- as.numeric(target) - 1
    ebic <- Rfast::poisson_only(dataset, target) + 2 * logn

  } else if ( identical(test, testIndPois) ) {  ## Poisson regression

    if ( ncores <= 1 | is.null(ncores) ) {
      ebic <- numeric(cols)
      for ( i in ina ) {
        fit2 <- glm( target ~ dataset[, i], poisson, weights = wei )
        ebic[i] <- BIC(fit2)
      }

    } else {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      ebic <- foreach(i = ina, .combine = rbind) %dopar% {
        fit2 <- glm( target ~ dataset[, i], poisson, weights = wei )
        return( BIC(fit2) )
      }
      stopCluster(cl)
    }
    
  } else if ( identical(test, testIndNB) ) {  ## Negative binomial regression

    if ( ncores <= 1 | is.null(ncores) ) {
      ebic <- numeric(cols)
      for ( i in 1:cols ) {
        fit2 <- MASS::glm.nb( target ~ dataset[, i], weights = wei )
        ebic[i] <- BIC(fit2)
      }

    } else {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      ebic <- foreach(i = ina, .combine = rbind, .packages = "MASS") %dopar% {
        fit2 <- MASS::glm.nb( target ~ dataset[, i], weights = wei )
        return( BIC(fit2) ) 
      }
      stopCluster(cl)
    }
    
  } else if ( identical(test, testIndNormLog)   ) {  ## Normal log link regression

    if ( ncores <= 1 | is.null(ncores) ) {
      ebic <- numeric(cols)
      for ( i in ina ) {
        fit2 <- glm( target ~ dataset[, i], family = gaussian(link = log), weights = wei )
        ebic[i] <- BIC(fit2)
      }

    } else {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      ebic <- foreach(i = ina, .combine = rbind) %dopar% {
        fit2 <- glm( target ~ dataset[, i], family = gaussian(link = log), weights = wei )
        return( BIC(fit2) )
      }
      stopCluster(cl)
    }
    
  } else if ( identical(test, testIndGamma)   ) {  ## Gamma regression

    if ( ncores <= 1 | is.null(ncores) ) {
      ebic <- numeric(cols)
      for ( i in ina ) {
        fit2 <- glm( target ~ dataset[, i], family = Gamma(link = log), weights = wei )
        ebic[i] <- BIC(fit2)
      }

    } else {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      ebic <- foreach(i = ina, .combine = rbind) %dopar% {
        fit2 <- glm( target ~ dataset[, i], family = Gamma(link = log), weights = wei )
        return( BIC(fit2) )
      }
      stopCluster(cl)
    } 
    
  } else if ( identical(test, testIndZIP) ) {  ## Zero-inflated Poisson regression
    ebic <- zip.regs(target, dataset, wei, logged = TRUE, ncores = ncores)[, 3] 

  } else if ( identical(test, testIndIGreg) ) {  ## Inverse Gaussian regression
    
    if ( ncores <= 1 | is.null(ncores) ) {
      ebic <- numeric(cols)
      for ( i in 1:cols ) {
        fit2 <- glm( target ~ dataset[, i], family = inverse.gaussian(link = log), weights = wei )
        ebic[i] <- BIC(fit2)
      }

    } else {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      ebic <- foreach(i = 1:cols, .combine = rbind) %dopar% {
        fit2 <- glm( target ~ dataset[, i], family = inverse.gaussian(link = log), weights = wei )
        return( BIC(fit2) )
      }
      stopCluster(cl)
    }
    
  } else if ( identical(test, censIndCR) ) {  ## Cox regression
    
    if ( ncores <= 1 | is.null(ncores) ) {
      ebic <- numeric(cols)
      for ( i in 1:cols ) {
        fit2 <- survival::coxph( target ~ dataset[, i], weights = wei)
        ebic[i] <- BIC(fit2)
      }

    } else {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      ebic <- foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
        fit2 <- survival::coxph( target ~ dataset[, i], weights = wei )
        return( BIC(fit2) )
      }
      stopCluster(cl)
    }
    
  } else if ( identical(test, censIndWR) ) {  ## Weibull regression

    if ( ncores <= 1 | is.null(ncores) ) {
      ebic <- numeric(cols)
      for ( i in 1:cols ) {
        fit2 <- survival::survreg( target ~ dataset[, i], weights = wei )
        ebic[i] <-  - 2 * logLik(fit2) + (length(fit2$coefficients) + 1) * logn 
      }

    } else {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      ebic <- foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
        fit2 <- survival::survreg( target ~ dataset[, i], weights = wei )
        return( - 2 * logLik(fit2) + (length(fit2$coefficients) + 1) * logn )
      }
      stopCluster(cl)
    }
    
  } else if ( identical(test, testIndTobit) ) {  ## Tobit regression

    if ( ncores <= 1 | is.null(ncores) ) {
      ebic <- numeric(cols)
      for ( i in 1:cols ) {
        fit2 <- survival::survreg( target ~ dataset[, i], weights = wei, dist = "gaussian" )
        ebic[i] <-  - 2 * logLik(fit2) + (length(fit2$coefficients) + 1) * logn
      }

    } else {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      ebic <- foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
        fit2 <- survival::survreg( target ~ dataset[, i], weights = wei, dist = "gaussian" )
        return( - 2 * logLik(fit2) + (length(fit2$coefficients) + 1) * logn )
      }
      stopCluster(cl)
    }
    
  } else if ( identical(test, testIndClogit) ) {  ## Conditional logistic regression
    subject <- target[, 2] #the patient id
    case <- as.logical(target[, 1]);  ## case 

    if ( ncores <= 1 | is.null(ncores) ) {
      ebic <- numeric(cols)
      for ( i in 1:cols ) {
        fit2 <- survival::clogit( case ~ dataset[, i] + strata(subject) ) 
        ebic[i] <- BIC(fit2)
      }

    } else {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      ebic <- foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
        fit2 <- survival::clogit(case ~ dataset[, i] + strata(subject) ) 
        return( BIC(fit2) )
      }
      stopCluster(cl)
    }
    
  } else if ( identical(test, censIndER) ) {  ## Exponential regression

    if ( ncores <= 1 | is.null(ncores) ) {
      ebic <- numeric(cols)
      for ( i in 1:cols ) {
        fit2 <- survival::survreg( target ~ dataset[, i], dist = "exponential", weights = wei )
        ebic[i] <-  - 2 * logLik(fit2) + (length(fit2$coefficients) + 1) * logn 
      }

    } else {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      ebic <- foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
        fit2 <- survival::survreg( target ~ dataset[, i], dist = "exponential", weights = wei )
        return( - 2 * logLik(fit2) + (length(fit2$coefficients) + 1) * logn )
      }
      stopCluster(cl)
    }
    
  }  else  ebic <- NULL
  
  if ( !is.null(ebic) )  {
    ebic <- as.numeric(ebic)
    if ( gam != 0 ) {
      ebic <- ebic + con * log(cols)
    } else ebic <- ebic
    if ( targetID != - 1 )   ebic[targetID] <- Inf
    if ( sum(id>0) > 0 )   ebic[id] <- Inf
  }
  
  list(ebic = ebic)
  
}