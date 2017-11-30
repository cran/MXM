rint.regs <- function(target, dataset, targetID = -1, id, reps = NULL, tol = 1e-08) {
  
  univariateModels <- list();
  dm <- dim(dataset)
  n <- dm[1]
  D <- dm[2]
  
  if (targetID != -1 ) {
    target <- dataset[, targetID]
    dataset[, targetID] <- rnorm(n)
  }   
  poia <- NULL
  poia <- Rfast::check_data(dataset)
  if ( sum(poia > 0) )  dataset[, poia] <- rnorm(n * length(poia) )

  if ( is.null(reps) ) {
    y <- target
    dataset <- as.matrix(dataset)
    mod <- Rfast::rint.regs(target, dataset, id, logged = TRUE)
    univariateModels$stat <- mod[, 1]
    univariateModels$pvalue <- mod[, 2]
  } else {
    
    sy <- as.vector( Rfast::group.sum(target, id) )
    ni <- tabulate(id)
    my <- sy / ni
    ni2 <- ni^2
    Xi <- cbind(1, reps, dataset[, 1])

    funa <- function(d, n, ni, ni2, S, hi2)   sum( log1p(ni * d) ) + n * log(S - d * sum(ni2 * hi2/ (1 + ni * d) ) )    
    stat <- numeric(D)
    for (i in 1:D) {
      Xi[, 3] <- dataset[, i]
      xx <- crossprod(Xi)
      sx <- rowsum(Xi, id)
      sxy <- crossprod(Xi, target)
      mx <- sx / ni
      #############
      mod <- .lm.fit(Xi, target)
      b1 <- mod$coefficients
      S <- sum( mod$residuals^2 )
      hi2 <- ( my - mx %*% b1 )^2
      mod <- optimise(funa, c(0, 70), n = n, ni = ni, ni2 = ni2, S = S, hi2 = hi2, tol = tol)
      d <- mod$minimum 
      b2 <- solve( xx - d * crossprod(sx/(1+ ni * d), sx), sxy - d * crossprod( sx, sy/(1 + ni * d) ) )   
      k <- 2
      while ( sum( abs(b2 - b1) ) > tol  &  k < 100) {
        k <- k + 1
        b1 <- b2
        S <- sum( (target - Xi %*% b1)^2 )
        hi2 <- ( my - mx %*% b1 )^2
        mod <- optimise(funa, c(0, 70), n = n, ni = ni, ni2 = ni2, S = S, hi2 = hi2, tol = tol)
        d <- mod$minimum 
        b2 <- solve( xx - d * crossprod(sx/(1+ ni * d), sx), sxy - d * crossprod( sx, sy/(1 + ni * d) ) ) 
      }
      se <- (S - d * sum(ni^2 * hi2/ (1 + ni * d) ) )/n
      seb <- solve( xx - d * crossprod(sx/(1+ ni * d), sx) )[3, 3] * se 
      stat[i] <- b2[3]^2 / seb
    }  
    univariateModels$stat = stat
    univariateModels$pvalue = pf(stat, 1, n - 5, lower.tail = FALSE, log.p = TRUE)
  
  }  ## else if ( is.null(reps) )
  if ( sum(poia>0) > 0 ) {
    univariateModels$stat[poia] = 0
    univariateModels$pvalue[poia] = log(1)
  }

  univariateModels
}   



