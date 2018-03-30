gomp2 <- function (y, x, tol = qchisq(0.95, 1) + log( length(y) ), type = "gamma" ) {
    tic <- proc.time()
    dm <- dim(x)
    d <- dm[2]
    n <- dm[1]
    ind <- 1:d
    x <- Rfast::standardise(x)
    can <- which( is.na(Rfast::colsums(x)) ) 
	  ind[can] <- 0
	
  if (type == "gamma") {
    
      tic <- proc.time()
      mod <- glm(y~ 1, family = Gamma(log))
      rho <- mod$deviance
      res <- y - fitted(mod)
      ela <- as.vector( cov(res, x) )
      sel <- which.max( abs(ela) )  
      sela <- sel
      names(sela) <- NULL
      mod <- glm(y ~ x[, sela], family = Gamma(log) )
      res <-  mod$residuals
      rho[2] <- mod$deviance
      phi <- summary(mod)[[ 14 ]]
      ind[sel] <- 0
      i <- 2
      r <- numeric(d)
      while ( (rho[i - 1] - rho[i]) / phi > tol ) {
        i <- i + 1
        r[ind] <- Rfast::colsums( x[, ind, drop = FALSE] * res )
        sel <- which.max( abs(r) )
        sela <- c(sela, sel)
        mod <- glm(y ~ x[, sela], family = Gamma(log) )
        res <- y - fitted(mod)
        rho[i] <-  - 2 * as.numeric( logLik(mod) ) 
        phi <- summary(mod)[[ 14 ]]
        ind[sela] <- 0
        r[sela] <- 0
      } ## end while ( (rho[i - 1] - rho[i]) > tol )

  } else if ( type == "quasipoisson" ) {
    m <- sum(y)/n
    rho <- 2 * sum(y * log(y), na.rm = TRUE) - 2 * n * m * log(m)
    ela <- as.vector( cor(y - m, x) )
    sel <- which.max( abs(ela) )
    sela <- sel
    names(sela) <- NULL
    options(warn = -1)
    mod <- Rfast::qpois.reg(x[, sel], y)
    phi <- mod$phi
    res <- y - exp( mod$be[1] + x[, sel] * mod$be[2] )
    rho[2] <- mod$devi
    ind[sel] <- 0
    i <- 2
    r <- numeric(d)
    while ( (rho[i - 1] - rho[i]) / phi > tol ) {
      i <- i + 1
      r[ind] <- colsums(res * x[, ind] )
      sel <- which.max( abs(r) )
      sela <- c(sela, sel)
      options(warn = -1)
      mod <- Rfast::qpois.reg(x[, sela], y)        
      res <- y - as.vector( exp( mod$be[1] + x[, sela] %*% mod$be[-1] ) ) 
      rho[i] <- mod$devi
      phi <- mod$phi 
      ind[sela] <- 0
      r[sela] <- 0
    }
	
  }
  runtime <- proc.time() - tic
  len <- length(sela)
  res <- cbind(c(0, sela[-len]), rho[1:len])
  colnames(res) <- c("Selected Vars", "Deviance") 
  list(runtime = runtime, res = res)
}

