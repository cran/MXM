omp <- function(target, dataset, tol = qchisq(0.95, 1) + log( dim(dataset)[1] ), test = "testIndFisher" ) {

  if ( test == "testIndReg" | test == "testIndFisher" ) {
    tic <- proc.time()
    res <- Rfast::ompr(target, dataset, method = "BIC", tol)
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, res = res)
  } else if ( test == "testIndLogistic" ) {
    tic <- proc.time()
    res <- omp2(target, dataset, tol, type = "logistic")
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, res = res)
  } else if ( test == "testIndPois" ) {
    tic <- proc.time()
    res <- omp2(target, dataset, tol, type = "poisson")
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, res = res)
  } else if (test == "testIndMVreg") {	
    tic <- proc.time()
    res <- omp2(target, dataset, tol, type = "mv")
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, res = res)
  } else if ( test == "testIndQBinom" ) {
    tic <- proc.time()
    res <- omp2(target, dataset, tol, type = "quasibinomial")	
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, res = res)
  } else if ( test == "testIndNormLog" ) {
    tic <- proc.time()
    res <- omp2(target, dataset, tol, type = "normlog")
    runtime <- proc.time() - tic
    result <- list(runtime = runtime, res = res)
    
	} else {
    d <- dim(dataset)[2]
    ind <- 1:d
    dataset <- Rfast::standardise(dataset)
    
    if (test == "testIndNB") {
      tic <- proc.time()
     	mod <- MASS::glm.nb(target ~ 1)
	    rho <-  - 2 * as.numeric( logLik(mod) ) 
	    res <-  target - fitted(mod)
	    ela <- as.vector( cor(res, dataset) )
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
       res <- cbind(c(0, sela), rho[1:(len+1)])
       colnames(res) <- c("Selected Vars", "Deviance") 
       result <- list(runtime = runtime, res = res)
       
    } else if (test == "testIndBeta") {
      tic <- proc.time()
	    mod <- Rfast::beta.mle(target)
      rho <-  - 2 * mod$loglik
	    res <- target - mod$param[1]/sum(mod$param)
      ela <- as.vector( cor(res, dataset) )
      sel <- which.max( abs(ela) )     
      sela <- sel
      names(sela) <- NULL
      mod <- beta.reg(target, dataset[, sela])
      est <- exp( mod$be[1] + dataset[, sela] * mod$be[2] )
      res <- target - est / (1 + est)
      rho[2] <-  - 2 * mod$loglik
      ind[sel] <- 0
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
      res <- cbind(c(0, sela), rho[1:(len+1)])
      colnames(res) <- c("Selected Vars", "Deviance") 
      result <- list(runtime = runtime, res = res)
      
    } else if ( test == "testIndGamma") {
      tic <- proc.time()
      mod <- glm(target ~ 1, family = Gamma(log))
	    rho <-  - 2 * as.numeric( logLik(mod) ) 
	    res <- target - fitted(mod)
	    ela <- as.vector( cor(res, dataset) )
      sel <- which.max( abs(ela) )  
      sela <- sel
      names(sela) <- NULL
      mod <- glm(target ~ dataset[, sela], family = Gamma(log) )
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
        mod <- glm(target ~ dataset[, sela], family = Gamma(log) )
        res <- target - fitted(mod)
        rho[i] <-  - 2 * as.numeric( logLik(mod) ) 
        ind[sela] <- 0
        r[sela] <- 0
      } ## end while ( (rho[i - 1] - rho[i]) > tol )
      runtime <- proc.time() - tic
      len <- length(sela)
      res <- cbind(c(0, sela), rho[1:(len+1)])
      colnames(res) <- c("Selected Vars", "Deviance") 
      result <- list(runtime = runtime, res = res)  
      
    } else if ( test == "testIndMMReg") {
      tic <- proc.time()
      mod <- MASS::rlm(target ~ 1, method = "MM", maxit = 2000)
      rho <-  - 2 * as.numeric( logLik(mod) ) 
      res <- mod$residuals
      ela <- as.vector( cor(res, dataset) )
      sel <- which.max( abs(ela) )  
      sela <- sel
      names(sela) <- NULL
      mod <- MASS::rlm(target ~ dataset[, sela], method = "MM", maxit = 2000 )
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
        mod <- MASS::rlm(target ~ dataset[, sela], method = "MM", maxit = 2000 )
        res <-  mod$residuals
        rho[i] <-  - 2 * as.numeric( logLik(mod) ) 
        ind[sela] <- 0
        r[sela] <- 0
      } ## end while ( (rho[i - 1] - rho[i]) > tol )
      runtime <- proc.time() - tic
      len <- length(sela)
      res <- cbind(c(0, sela), rho[1:(len+1)])
      colnames(res) <- c("Selected Vars", "Deviance") 
      result <- list(runtime = runtime, res = res)  
      
    } else if ( test == "testIndTobit") {
      tic <- proc.time()
      mod <- survival::survreg(target ~ 1, dist = "gaussian")
      rho <-  - 2 * as.numeric( logLik(mod) ) 
	    res <- target - coef(mod)
      ela <- as.vector( cor(res, dataset) )
      sel <- which.max( abs(ela) )  
      sela <- sel
      names(sela) <- NULL
      mod <- survival::survreg(target ~ dataset[, sela], dist = "gaussian" )
      res <- target - predict( mod, data.frame(dataset[, sela]) ) 
      rho[2] <-  - 2 * as.numeric( logLik(mod) ) 
      ind[sel] <- 0
      i <- 2
      r <- numeric(d)
      while ( (rho[i - 1] - rho[i]) > tol ) {
        i <- i + 1
        r[ind] <- Rfast::colsums( dataset[, ind, drop = FALSE] * res )
        sel <- which.max( abs(r) )
        sela <- c(sela, sel)
        mod <- survival::survreg(target ~ dataset[, sela], dist = "gaussian" )
        res <- target - predict( mod, data.frame(dataset[, sela]) ) 
        rho[i] <-  - 2 * as.numeric( logLik(mod) ) 
        ind[sela] <- 0
        r[sela] <- 0
      } ## end while ( (rho[i - 1] - rho[i]) > tol )
      runtime <- proc.time() - tic
      len <- length(sela)
      res <- cbind(c(0, sela), rho[1:(len+1)])
      colnames(res) <- c("Selected Vars", "Deviance") 
      result <- list(runtime = runtime, res = res)

    } else if ( test == "censIndCR") {
      tic <- proc.time()
      mod <- survival::coxph(target ~ 1)
      rho <-  - 2 * summary( mod) [[1]]
      res <- mod$residuals   ## martingale residuals
      ela <- as.vector( cor(res, dataset) )
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
      res <- cbind(c(0, sela), rho[1:(len+1)])
      colnames(res) <- c("Selected Vars", "Deviance") 
      result <- list(runtime = runtime, res = res)
    
    } else if ( test == "censIndWR") {
      tic <- proc.time()
      mod <- survival::survreg(target ~ 1)
      rho <-  - 2 * as.numeric( logLik(mod) ) 
      res <- target - coef(mod)
      ela <- as.vector( cor(res, dataset) )
      sel <- which.max( abs(ela) )  
      sela <- sel
      names(sela) <- NULL
      mod <- survival::survreg(target ~ dataset[, sela] )
      res <- target - predict( mod, data.frame(dataset[, sela]) ) 
      rho[2] <-  - 2 * as.numeric( logLik(mod) ) 
      ind[sel] <- 0
      i <- 2
      r <- numeric(d)
      while ( (rho[i - 1] - rho[i]) > tol ) {
        i <- i + 1
        r[ind] <- Rfast::colsums( dataset[, ind, drop = FALSE] * res )
        sel <- which.max( abs(r) )
        sela <- c(sela, sel)
        mod <- survival::survreg(target ~ dataset[, sela] )
        res <- target - predict( mod, data.frame(dataset[, sela]) ) 
        rho[i] <-  - 2 * as.numeric( logLik(mod) ) 
        ind[sela] <- 0  
        r[sela] <- 0
      } ## end while ( (rho[i - 1] - rho[i]) > tol )
      runtime <- proc.time() - tic
      len <- length(sela)
      res <- cbind(c(0, sela), rho[1:(len+1)])
      colnames(res) <- c("Selected Vars", "Deviance") 
      result <- list(runtime = runtime, res = res)
    } ##  end if (test == "testIndNB")
  }  ##  end if ( test == "testIndReg" | test == "testIndfisher" ) 
  
  result  
}
