gomp <- function(target, dataset, tol = qchisq(0.95, 1) + log( dim(dataset)[1] ), test = "testIndFisher" ) {

  if ( test == "testIndReg" | test == "testIndFisher" ) {
    tic <- proc.time()
    res <- Rfast::ompr(target, dataset, method = "BIC", tol)
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, res = res)
  } else if ( test == "testIndLogistic" ) {
    tic <- proc.time()
    res <- Rfast::omp(target, dataset, tol, type = "logistic")$info
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, res = res)
  } else if ( test == "testIndPois" ) {
    tic <- proc.time()
    res <- Rfast::omp(target, dataset, tol, type = "poisson")$info
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, res = res)
    
  } else if ( test == "testIndQPois" ) {
    result <- gomp2(target, dataset, tol, type = "quasipoisson")
    
  } else if (test == "testIndMVreg") {	
    tic <- proc.time()
    res <- Rfast::omp(target, dataset, tol, type = "mv")$info
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, res = res)
  } else if ( test == "testIndQBinom" ) {
    tic <- proc.time()
    res <- Rfast::omp(target, dataset, tol, type = "quasibinomial")$info
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, res = res)
  } else if ( test == "testIndNormLog" ) {
    tic <- proc.time()
    res <- Rfast::omp(target, dataset, tol, type = "normlog")$info
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, res = res)
    
  } else if ( test == "testIndGamma" ) {
    result <- gomp2(target, dataset, tol, type = "gamma")
    
  } else {
    d <- dim(dataset)[2]
    ind <- 1:d
    dataset <- Rfast::standardise(dataset)
    can <- which( is.na( Rfast::colsums(dataset) ) )
	  ind[can] <- 0
	
    if (test == "testIndNB") {
      tic <- proc.time()
     	mod <- MASS::glm.nb(target ~ 1)
    	rho <-  - 2 * as.numeric( logLik(mod) ) 
	    res <-  target - fitted(mod)
	    ela <- as.vector( cov(res, dataset) )
	    sel <- which.max( abs(ela) )  
      sela <- sel
      names(sela) <- NULL
      mod <- MASS::glm.nb(target ~ dataset[, sela])
      res <-  mod$residuals
      rho[2] <-  - 2 * as.numeric( logLik(mod) )
      ind[sel] <- 0
      i <- 2
      r <- numeric(d)
      while ( (rho[i - 1] - rho[i]) > tol ) {
        i <- i + 1
        r[ind] <- Rfast::colsums( dataset[, ind, drop = FALSE] * res )
        sel <- which.max( abs(r) )
        sela <- c(sela, sel)
        mod <- MASS::glm.nb(target ~ dataset[, sela])        
        res <- target - fitted(mod)
        rho[i] <-  - 2 * as.numeric( logLik(mod) )
        ind[sela] <- 0
        r[sela] <- 0
      } ## end while ( (rho[i - 1] - rho[i]) > tol )
      runtime <- proc.time() - tic
      len <- length(sela)
      res <- cbind(c(0, sela[-len]), rho[1:len])
      colnames(res) <- c("Selected Vars", "Deviance") 
      result <- list(runtime = runtime, res = res)
       
    } else if (test == "testIndBeta") {
      tic <- proc.time()
	    mod <- Rfast::beta.mle(target)
      rho <-  - 2 * mod$loglik
	    res <- target - mod$param[1]/sum(mod$param)
      ela <- as.vector( cov(res, dataset) )
      sel <- which.max( abs(ela) )     
      sela <- sel
      names(sela) <- NULL
      mod <- try( beta.reg(target, dataset[, sela]), silent = TRUE )
      if ( identical( class(mod), "try-error" ) ) {
        rho[2] <- rho[1]  
      } else {  
        est <- exp( mod$be[1] + dataset[, sela] * mod$be[2] )
        res <- target - est / (1 + est)
        rho[2] <-  - 2 * mod$loglik
        ind[sel] <- 0
      }  
      i <- 2
      r <- numeric(d)
      while ( (rho[i - 1] - rho[i]) > tol ) {
        i <- i + 1
        r[ind] <- Rfast::colsums( dataset[, ind, drop = FALSE] * res )
        sel <- which.max( abs(r) )
        sela <- c(sela, sel)
        mod <- try( beta.reg(target, dataset[, sela]), silent = TRUE )        
        if ( identical( class(mod), "try-error" ) ) {
          rho[i] <- rho[i - 1]  
        } else {  
          est <- exp( mod$be[1] + dataset[, sela] %*% mod$be[-1] )
          res <- target - est / (1 + est)
          rho[i] <-  - 2 * mod$loglik
          ind[sela] <- 0
          r[sela] <- 0
        }  
      } ## end while ( (rho[i - 1] - rho[i]) > tol )
      runtime <- proc.time() - tic
      len <- length(sela)
      res <- cbind(c(0, sela[-len]), rho[1:len])
      colnames(res) <- c("Selected Vars", "Deviance") 
      result <- list(runtime = runtime, res = res)
      
    } else if ( test == "testIndMMReg") {
      tic <- proc.time()
      mod <- MASS::rlm(target ~ 1, method = "MM", maxit = 2000)
      rho <-  - 2 * as.numeric( logLik(mod) ) 
      res <- mod$residuals
      ela <- as.vector( cov(res, dataset) )
      sel <- which.max( abs(ela) )  
      sela <- sel
      names(sela) <- NULL
      mod <- MASS::rlm(target ~ dataset[, sela], method = "MM", maxit = 2000 )
      res <- mod$residuals
      rho[2] <-  - 2 * as.numeric( logLik(mod) ) 
      ind[sel] <- 0
      i <- 2
      r <- numeric(d)
      while ( (rho[i - 1] - rho[i]) > tol ) {
        i <- i + 1
        r[ind] <- Rfast::colsums( dataset[, ind, drop = FALSE] * res )
        sel <- which.max( abs(r) )
        sela <- c(sela, sel)
        mod <- try( MASS::rlm(target ~ dataset[, sela], method = "MM", maxit = 2000 ), silent = TRUE )
        if ( identical( class(mod), "try-error" ) ) {
          rho[i] <- rho[i - 1]
        } else {
          res <- mod$residuals
          rho[i] <-  - 2 * as.numeric( logLik(mod) ) 
          ind[sela] <- 0
          r[sela] <- 0
        }  
      } ## end while ( (rho[i - 1] - rho[i]) > tol )
      runtime <- proc.time() - tic
      len <- length(sela)
      res <- cbind(c(0, sela[-len]), rho[1:len])
      colnames(res) <- c("Selected Vars", "Deviance") 
      result <- list(runtime = runtime, res = res)
      
    } else if ( test == "testIndRQ") {
      tic <- proc.time()
      mod <- quantreg::rq(target ~ 1)
      rho <-  - 2 * as.numeric( logLik(mod) ) 
      res <- mod$residuals
      ela <- as.vector( cov(res, dataset) )
      sel <- which.max( abs(ela) )  
      sela <- sel
      names(sela) <- NULL
      mod <- quantreg::rq(target ~ dataset[, sela])
      res <- mod$residuals
      rho[2] <-  - 2 * as.numeric( logLik(mod) ) 
      ind[sel] <- 0
      i <- 2
      r <- numeric(d)
      while ( (rho[i - 1] - rho[i]) > tol ) {
        i <- i + 1
        r[ind] <- Rfast::colsums( dataset[, ind, drop = FALSE] * res )
        sel <- which.max( abs(r) )
        sela <- c(sela, sel)
        mod <- try( quantreg::rq(target ~ dataset[, sela]), silent = TRUE )
        if ( identical( class(mod), "try-error" ) ) {
          rho[i] <- rho[i - 1]
        } else {
          res <- mod$residuals
          rho[i] <-  - 2 * as.numeric( logLik(mod) ) 
          ind[sela] <- 0
          r[sela] <- 0
        }  
      } ## end while ( (rho[i - 1] - rho[i]) > tol )
      runtime <- proc.time() - tic
      len <- length(sela)
      res <- cbind(c(0, sela[-len]), rho[1:len])
      colnames(res) <- c("Selected Vars", "Deviance") 
      result <- list(runtime = runtime, res = res)
      
    } else if ( test == "testIndTobit") {
      tic <- proc.time()
      mod <- survival::survreg(target ~ 1, dist = "gaussian")
      rho <-  - 2 * as.numeric( logLik(mod) ) 
	    res <- resid(mod)
      ela <- as.vector( cov(res, dataset) )
      sel <- which.max( abs(ela) )  
      sela <- sel
      names(sela) <- NULL
      mod <- survival::survreg(target ~ dataset[, sela], dist = "gaussian" )
	    res <- resid(mod)
      rho[2] <-  - 2 * as.numeric( logLik(mod) ) 
      ind[sel] <- 0
      i <- 2
      r <- numeric(d)
      while ( (rho[i - 1] - rho[i]) > tol ) {
        i <- i + 1
        r[ind] <- Rfast::colsums( dataset[, ind, drop = FALSE] * res )
        sel <- which.max( abs(r) )
        sela <- c(sela, sel)
        mod <- try( survival::survreg(target ~ dataset[, sela], dist = "gaussian" ), silent = TRUE )
        if ( identical( class(mod), "try-error" ) ) {
          rho[i] <- rho[i - 1]
        } else {
          res <- resid(mod)
          rho[i] <-  - 2 * as.numeric( logLik(mod) ) 
          ind[sela] <- 0
          r[sela] <- 0
        }  
      } ## end while ( (rho[i - 1] - rho[i]) > tol )
      runtime <- proc.time() - tic
      len <- length(sela)
      res <- cbind(c(0, sela[-len]), rho[1:len])
      colnames(res) <- c("Selected Vars", "Deviance") 
      result <- list(runtime = runtime, res = res)

    } else if ( test == "censIndCR") {
      tic <- proc.time()
      mod <- survival::coxph(target ~ 1)
      rho <-  - 2 * summary( mod) [[1]]
      res <- mod$residuals   ## martingale residuals
      ela <- as.vector( cov(res, dataset) )
      sel <- which.max( abs(ela) )  
      sela <- sel
      names(sela) <- NULL
      mod <- survival::coxph(target ~ dataset[, sela] )
      res <-  mod$residuals
      rho[2] <-  - 2 * as.numeric( logLik(mod) ) 
      ind[sel] <- 0
      i <- 2
      r <- numeric(d)
      while ( (rho[i - 1] - rho[i]) > tol ) {
        i <- i + 1
        r[ind] <- Rfast::colsums( dataset[, ind, drop = FALSE] * res )
        sel <- which.max( abs(r) )
        sela <- c(sela, sel)
        mod <- try( survival::coxph(target ~ dataset[, sela] ), silent = TRUE )
        if ( identical( class(mod), "try-error" ) ) {
          rho[i] <- rho[i - 1]
        } else {
          res <- mod$residuals   ## martingale residuals
          rho[i] <-  - 2 * as.numeric( logLik(mod) ) 
          ind[sela] <- 0
          r[sela] <- 0
        }  
      } ## end while ( (rho[i - 1] - rho[i]) > tol )
      runtime <- proc.time() - tic
      len <- length(sela)
      res <- cbind(c(0, sela[-len]), rho[1:len])
      colnames(res) <- c("Selected Vars", "Deviance") 
      result <- list(runtime = runtime, res = res)
    
    } else if ( test == "censIndWR") {
      tic <- proc.time()
      mod <- survival::survreg(target ~ 1)
      rho <-  - 2 * as.numeric( logLik(mod) ) 
      res <- resid(mod)
      ela <- as.vector( cov(res, dataset) )
      sel <- which.max( abs(ela) )  
      sela <- sel
      names(sela) <- NULL
      mod <- survival::survreg(target ~ dataset[, sela] )
      res <- resid(mod)
      rho[2] <-  - 2 * as.numeric( logLik(mod) ) 
      ind[sel] <- 0
      i <- 2
      r <- numeric(d)
      while ( (rho[i - 1] - rho[i]) > tol ) {
        i <- i + 1
        r[ind] <- Rfast::colsums( dataset[, ind, drop = FALSE] * res )
        sel <- which.max( abs(r) )
        sela <- c(sela, sel)
        mod <- try( survival::survreg(target ~ dataset[, sela] ), silent = TRUE )
        if ( identical( class(mod), "try-error" ) ) {
          rho[i] <- rho[i - 1]
        } else {  
          res <- resid(mod)
          rho[i] <-  - 2 * as.numeric( logLik(mod) ) 
          ind[sela] <- 0  
          r[sela] <- 0
        }  
      } ## end while ( (rho[i - 1] - rho[i]) > tol )
      runtime <- proc.time() - tic
      len <- length(sela)
      res <- cbind(c(0, sela[-len]), rho[1:len])
      colnames(res) <- c("Selected Vars", "Deviance") 
      result <- list(runtime = runtime, res = res)
    } ##  end if (test == "testIndNB")
  }  ##  end if ( test == "testIndReg" | test == "testIndfisher" ) 
  
  result  
}
