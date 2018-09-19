gomp2 <- function (y, x, tol = qchisq(0.95, 1) + log( length(y) ), test = "testIndGamma" ) {
    tic <- proc.time()
    dm <- dim(x)
    d <- dm[2]
    ind <- 1:d
    x <- Rfast::standardise(x)
    can <- which( is.na(Rfast::colsums(x)) ) 
    ind[can] <- 0
    tic <- proc.time()
    mod <- glm(y ~ 1, family = Gamma(log))
    rho <- mod$deviance
    phi[1] <- summary(mod)[[ 14 ]]
    res <- y - fitted(mod)
    ela <- as.vector( cov(res, x) )
    sel <- which.max( abs(ela) )  
    sela <- sel
    names(sela) <- NULL
    mod <- glm(y ~ x[, sela], family = Gamma(log) )
    res <-  mod$residuals
    rho[2] <- mod$deviance
    phi[2] <- summary(mod)[[ 14 ]]
    ind[sel] <- 0
    i <- 2
    while ( (rho[i - 1] - rho[i]) / phi[i] > tol ) {
      r <- rep(NA, d)
      i <- i + 1
      ##r[ind] <- Rfast::colsums( x[, ind, drop = FALSE] * res )
      r[ind] <- Rfast::eachcol.apply(x, res, indices = ind[ind > 0 ], oper = "*", apply = "sum") 
      sel <- which.max( abs(r) )
      sela <- c(sela, sel)
      mod <- glm(y ~ x[, sela], family = Gamma(log) )
      res <- y - fitted(mod)
      rho[i] <- mod$deviance
      phi[i] <- summary(mod)[[ 14 ]]
      ind[sela] <- 0
    } ## end while ( (rho[i - 1] - rho[i]) > tol )

  runtime <- proc.time() - tic
  len <- length(sela)
  res <- cbind(c(0, sela[-len]), rho[1:len])
  colnames(res) <- c("Selected Vars", "Deviance") 
  list(runtime = runtime, phi = phi[1:len], res = res)
}

