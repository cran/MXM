gomp <- function(target, dataset, tol = qchisq(0.95, 1), test = "testIndLogistic", method = "ar2") {

  if ( test == "testIndReg" | test == "testIndFisher" ) {
    tic <- proc.time()
    res <- Rfast::ompr(target, dataset, method = method, tol)
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = NULL, res = res)
  } else if ( test == "testIndLogistic" ) {
    tic <- proc.time()
    res <- Rfast::omp(target, dataset, tol, type = "logistic")
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = NULL, res = res$info)
  } else if ( test == "testIndPois" ) {
    tic <- proc.time()
    res <- Rfast::omp(target, dataset, tol, type = "poisson")
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = NULL, res = res$info)

  } else if ( test == "testIndQPois" ) {
    tic <- proc.time()
    res <- Rfast::omp(target, dataset, tol, type = "quasipoisson")
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = res$phi, res = res$info)
    
  } else if (test == "testIndMVreg") {	
    tic <- proc.time()
    res <- Rfast::omp(target, dataset, tol, type = "mv")
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = NULL, res = res$info)
  } else if ( test == "testIndQBinom" ) {
    tic <- proc.time()
    res <- Rfast::omp(target, dataset, tol, type = "quasibinomial")
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = res$phi, res = res$info)
  } else if ( test == "testIndNormLog" ) {
    tic <- proc.time()
    res <- Rfast::omp(target, dataset, tol, type = "normlog")
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, phi = res$phi, res = res$info)
    
  } else if ( test == "testIndGamma" ) {
    result <- gomp2(target, dataset, tol, test = test)
    
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
      mod <- MASS::glm.nb( target ~ dataset[, sela], control = list(epsilon = 1e-08, maxit = 50, trace = FALSE) )
      res <-  mod$residuals
      rho[2] <-  - 2 * as.numeric( logLik(mod) )
      ind[sel] <- 0
      i <- 2
      r <- rep(NA, d)
      while ( (rho[i - 1] - rho[i]) > tol ) {
        i <- i + 1
        ## r[ind] <- Rfast::colsums( dataset[, ind, drop = FALSE] * res )
        r[ind] <- Rfast::eachcol.apply(dataset, res, indices = ind[ind > 0 ], oper = "*", apply = "sum") 
        sel <- which.max( abs(r) )
        sela <- c(sela, sel)
        mod <- MASS::glm.nb( target ~ dataset[, sela], control = list(epsilon = 1e-08, maxit = 50, trace = FALSE) )        
        res <- target - fitted(mod)
        rho[i] <-  - 2 * as.numeric( logLik(mod) )
        ind[sela] <- 0
        r[sela] <- NA
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
      r <- rep(NA, d)
      while ( (rho[i - 1] - rho[i]) > tol ) {
        i <- i + 1
        ##r[ind] <- Rfast::colsums( dataset[, ind, drop = FALSE] * res )
        r[ind] <- Rfast::eachcol.apply(dataset, res, indices = ind[ind > 0 ], oper = "*", apply = "sum") 
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
          r[sela] <- NA
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
      r <- rep(NA, d)
      while ( (rho[i - 1] - rho[i]) > tol ) {
        i <- i + 1
        ##r[ind] <- Rfast::colsums( dataset[, ind, drop = FALSE] * res )
        r[ind] <- Rfast::eachcol.apply(dataset, res, indices = ind[ind > 0 ], oper = "*", apply = "sum") 
        sel <- which.max( abs(r) )
        sela <- c(sela, sel)
        mod <- try( MASS::rlm(target ~ dataset[, sela], method = "MM", maxit = 2000 ), silent = TRUE )
        if ( identical( class(mod), "try-error" ) ) {
          rho[i] <- rho[i - 1]
        } else {
          res <- mod$residuals
          rho[i] <-  - 2 * as.numeric( logLik(mod) ) 
          ind[sela] <- 0
          r[sela] <- NA
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
      r <- rep(NA, d)
      while ( (rho[i - 1] - rho[i]) > tol ) {
        i <- i + 1
        ##r[ind] <- Rfast::colsums( dataset[, ind, drop = FALSE] * res )
        r[ind] <- Rfast::eachcol.apply(dataset, res, indices = ind[ind > 0 ], oper = "*", apply = "sum") 
        sel <- which.max( abs(r) )
        sela <- c(sela, sel)
        mod <- try( quantreg::rq(target ~ dataset[, sela]), silent = TRUE )
        if ( identical( class(mod), "try-error" ) ) {
          rho[i] <- rho[i - 1]
        } else {
          res <- mod$residuals
          rho[i] <-  - 2 * as.numeric( logLik(mod) ) 
          ind[sela] <- 0
          r[sela] <- NA
        }  
      } ## end while ( (rho[i - 1] - rho[i]) > tol )
      runtime <- proc.time() - tic
      len <- length(sela)
      res <- cbind(c(0, sela[-len]), rho[1:len])
      colnames(res) <- c("Selected Vars", "Deviance") 
      result <- list(runtime = runtime, res = res)
      
    } else if (test == "testIndOrdinal") {
      tic <- proc.time()
      mod <- MASS::polr(target ~ 1)
      rho <-  - 2 * as.numeric( logLik(mod) )
      res <- ord.resid(target, mod$fitted.values)
      ela <- as.vector( cov(res, dataset) )
      sel <- which.max( abs(ela) )     
      sela <- sel
      names(sela) <- NULL
      mod <- MASS::polr(target ~ dataset[, sela])
      res <- ord.resid(target, mod$fitted.values) 
      rho[2] <-  - 2 * as.numeric( logLik(mod) )
      ind[sel] <- 0
      i <- 2
      r <- rep(NA, d)
      while ( (rho[i - 1] - rho[i]) > tol ) {
        i <- i + 1
        ##r[ind] <- Rfast::colsums( dataset[, ind, drop = FALSE] * res )
        r[ind] <- Rfast::eachcol.apply(dataset, res, indices = ind[ind > 0 ], oper = "*", apply = "sum") 
        sel <- which.max( abs(r) )
        sela <- c(sela, sel)
        mod <- MASS::polr(target ~ dataset[, sela])        
        res <- ord.resid(target, mod$fitted.values)
        rho[i] <-  - 2 * as.numeric( logLik(mod) )
        ind[sela] <- 0
        r[sela] <- NA
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
      r <- rep(NA, d)
      while ( (rho[i - 1] - rho[i]) > tol ) {
        i <- i + 1
        ##r[ind] <- Rfast::colsums( dataset[, ind, drop = FALSE] * res )
        r[ind] <- Rfast::eachcol.apply(dataset, res, indices = ind[ind > 0 ], oper = "*", apply = "sum") 
        sel <- which.max( abs(r) )
        sela <- c(sela, sel)
        mod <- try( survival::survreg(target ~ dataset[, sela], dist = "gaussian" ), silent = TRUE )
        if ( identical( class(mod), "try-error" ) ) {
          rho[i] <- rho[i - 1]
        } else {
          res <- resid(mod)
          rho[i] <-  - 2 * as.numeric( logLik(mod) ) 
          ind[sela] <- 0
          r[sela] <- NA
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
      r <- rep(NA, d)
      while ( (rho[i - 1] - rho[i]) > tol ) {
        i <- i + 1
        ##r[ind] <- Rfast::colsums( dataset[, ind, drop = FALSE] * res )
        r[ind] <- Rfast::eachcol.apply(dataset, res, indices = ind[ind > 0 ], oper = "*", apply = "sum") 
        sel <- which.max( abs(r) )
        sela <- c(sela, sel)
        mod <- try( survival::coxph(target ~ dataset[, sela] ), silent = TRUE )
        if ( identical( class(mod), "try-error" ) ) {
          rho[i] <- rho[i - 1]
        } else {
          res <- mod$residuals   ## martingale residuals
          rho[i] <-  - 2 * as.numeric( logLik(mod) ) 
          ind[sela] <- 0
          r[sela] <- NA
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
      mod <- survival::survreg(target ~ dataset[, sela], control = list(iter.max = 5000) )
      res <- resid(mod)
      rho[2] <-  - 2 * as.numeric( logLik(mod) ) 
      if ( is.na(rho[2]) )  rho[2] <- rho[1]
      ind[sel] <- 0
      i <- 2
      r <- rep(NA, d)
      while ( (rho[i - 1] - rho[i]) > tol ) {
        i <- i + 1
        ##r[ind] <- Rfast::colsums( dataset[, ind, drop = FALSE] * res )
        r[ind] <- Rfast::eachcol.apply(dataset, res, indices = ind[ind > 0 ], oper = "*", apply = "sum") 
        sel <- which.max( abs(r) )
        sela <- c(sela, sel)
        mod <- try( survival::survreg(target ~ dataset[, sela], control = list(iter.max = 5000) ), silent = TRUE )
        if ( identical( class(mod), "try-error" ) ) {
          rho[i] <- rho[i - 1]
        } else {  
          res <- resid(mod)
          rho[i] <-  - 2 * as.numeric( logLik(mod) ) 
          ind[sela] <- 0  
          r[sela] <- NA
        }  
      } ## end while ( (rho[i - 1] - rho[i]) > tol )
      runtime <- proc.time() - tic
      len <- length(sela)
      res <- cbind(c(0, sela[-len]), rho[1:len])
      colnames(res) <- c("Selected Vars", "Deviance") 
      result <- list(runtime = runtime, phi = NULL, res = res)
    } ##  end if (test == "testIndNB")
    
  }  ##  end if ( test == "testIndReg" | test == "testIndfisher" ) 
  
  result  
}
