zip.regs <- function(target, dataset, wei = NULL, logged = FALSE, ncores = 1) {
  
  if ( ncores <= 1 ) {
    D <- ncol(dataset)
    n <- length(target)
    dof <- loglik <- numeric(D)   
    poia <- which( target == 0 ) 
    n0 <- length(poia)    ;     n1 <-  n - n0
    target1 <- target[ -poia ]   

    if ( is.null(wei) ) {
      lgy <- sum( lgamma(target1 + 1) )  
      ini <- 2 * Rfast::zip.mle(target)$loglik

      for (i in 1:D) { 
        x <- model.matrix(target ~ ., as.data.frame(dataset[, i]) )[1:n, ]
        mod <- glm.fit(x[ - poia, ], target1, family = poisson(log) ) 
        p1 <- ( n0 - sum( exp( - fitted(mod) ) ) ) / n
        g1 <- log( p1 / ( 1 - p1 ) ) 
        lik <- nlm( regzip, c(g1, coef(mod) ), y1 = target1, x = x, n1 = n1, poia = poia, iterlim = 10000 )
        lik2 <- nlm( regzip, lik$estimate, y1 = target1, x = x, n1 = n1, poia = poia, iterlim = 10000 )
        loglik[i] <-  - lik2$minimum - lgy
        dof[i] <- dim(x)[2] - 1
      }
    } else { 
      wei <- wei / sum(wei)
      w0 <- wei[poia]    ;   w1 <- wei[-poia] 
      lgy <- sum( lgamma(target1 + 1) )  
      ini <- 2 * zipmle.wei(target, wei)$loglik

      for (i in 1:D) {
        x <- model.matrix(target ~ ., as.data.frame(dataset[, i]) )[1:n, ]
        mod <- glm.fit( x[ -poia, ], target1, family = poisson(log), weights = w1 ) 
        p1 <- ( n0 - sum( exp( - fitted(mod) ) ) ) / n
        g1 <- log( p1 / ( 1 - p1 ) )
        lik <- nlm( regzipwei, c(g1, coef(mod) ), y1 = target1, x = x, n1 = n1, w1 = w1, w0 = w0, poia = poia, iterlim = 10000 )
        lik2 <- nlm( regzipwei, lik$estimate, y1 = target1, x = x, n1 = n1, w1 = w1, w0 = w0, poia = poia, iterlim = 10000 )  
      }

      loglik[i] <-  - lik2$minimum - lgy
      dof[i] <- dim(x)[2] - 1

    }
    
    bic <-  - 2 *loglik + dof * log(n)
    stat <- 2 * loglik - ini
    pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = logged) 
    
  } else {
    
    if ( is.null(wei) ) {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      poia <- which( target == 0 ) 
      D <- ncol(dataset)
      n <- length(target)
      n0 <- length(poia)    
      n1 <-  n - n0
      target1 <- target[ -poia ]   
      lgy <- sum( lgamma(target1 + 1) )  
      ini <- 2 * Rfast::zip.mle(target)[3]
      #ini <- 2 * Rfast::zip.mle(target)$loglik

      mod <- foreach( i = 1:D, .combine = rbind, .export = "regzip") %dopar% {
    
        x <- model.matrix(target ~ ., as.data.frame(dataset[, i]) )[1:n, ]
  
        mod2 <- glm.fit(x[ - poia, ], target1, family = poisson(log) ) 
        p1 <- ( n0 - sum( exp( - fitted(mod2) ) ) ) / n
        g1 <- log( p1 / ( 1 - p1 ) )  
        lik <- nlm( regzip, c(g1, coef(mod2) ), y1 = target1, x = x, n1 = n1, poia = poia, iterlim = 10000 )
        lik2 <- nlm( regzip, lik$estimate, y1 = target1, x = x, n1 = n1, poia = poia, iterlim = 10000 )
  
        return( c( - lik2$minimum - lgy, dim(x)[2] - 1 ) )
      }

      stopCluster(cl)

    } else {  
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      poia <- which( target == 0 ) 
      D <- ncol(dataset)
      n <- length(target)
      n0 <- length(poia)    
      n1 <-  n - n0
      target1 <- target[ -poia ]   
      wei <- wei / sum(wei)
      w0 <- wei[poia]    ;   w1 <- wei[-poia] 
      ini <- 2 * zipmle.wei(target, wei)$loglik
      lgy <- sum( w1 * lgamma(target1 + 1) )  
         
      mod <- foreach( i = 1:D, .combine = rbind, .export = "regzipwei" ) %dopar% {

        mod2 <- glm.fit( x[ -poia, ], target1, family = poisson(log), weights = w1 ) 
        p1 <- ( n0 - sum( exp( - fitted(mod2) ) ) ) / n
        g1 <- log( p1 / ( 1 - p1 ) )
        lik <- nlm( regzipwei, c(g1, coef(mod2) ), y1 = target1, x = x, n1 = n1, w1 = w1, w0 = w0, iterlim = 10000 )
        lik2 <- nlm( regzipwei, lik$estimate, y1 = target1, x = x, n1 = n1, w1 = w1, w0 = w0, iterlim = 10000 )  
   
        return( c( - lik2$minimum - lgy, dim(x)[2] - 1 ) )
      }

    }
    
    bic <-  - 2 * mod[, 1] + mod[, 2] * log(n) 
    stat <- 2 * mod[, 1] - ini
    pvalue <- pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = logged) 
  }  

  cbind(stat, pvalue, bic)

}





   
   ### zero inflated log-likelihood with weights but no covariates
   
   zipmle.wei <- function(y, wei) {
     
     n <- length(y)
     poia <- which(y == 0)
     n0 <- length(poia)   ;    n1 <-  n - n0
     y1 <- y[ -poia ]    
     
     wei <- wei / sum(wei)
     w0 <- wei[poia]     ;    w1 <- wei[-poia]
     
     pini <- n0 / n1
     pini <- log(pini / (1 - pini) )
     lik <- nlm( zipwei, c(pini, log( sum(y1) / n ) ), y1 = y1, w1 = w1, w0 = w0, iterlim = 10000 )
     lik2 <- nlm( zipwei, lik$estimate, y1 = y1, w1 = w1, w0 = w0, iterlim = 10000 )  
     
     prop <- exp(lik2$estimate[1]) / ( 1 + exp(lik2$estimate[1]) )
     
     list(iters = lik$iterations + lik2$iterations, prop = prop, lam = exp(lik2$estimate[2]),
          loglik = -lik2$minimum - sum( w1 * lgamma(y1 + 1) ) )
     
   }
  
   
   

   ## objective functions used in the optimisation procedure
   
   regzip <- function(pa, y1, x, n1, poia) {
     a <- exp(pa[1]) / ( 1 + exp(pa[1]) )    ;    b <- pa[-1]
     es <- tcrossprod(b, x)
     - sum( log( a + (1 - a) * exp( - exp( es[poia]) ) ) ) - n1 * log( 1 - a ) - sum( y1 * es[ -poia ] - exp(es[ -poia ]) ) 
   }  
   
   
   regzipwei <- function(pa, y1, x, n1, w1, w0, poia) {
     a <- exp(pa[1]) / ( 1 + exp(pa[1]) )    ;    b <- pa[-1]
     es <- tcrossprod(b, x)
     - sum( w0 * log( a + (1 - a) * exp( - exp( es[poia]) ) ) ) - sum(w1 * log( 1 - a ) ) - sum( w1 * y1 * es[ -poia ] - w1 * exp(es[ -poia ]) ) 
   }  
   
   
   zipwei <- function(pa, y1, w1, w0) {
     a <- exp(pa[1]) / ( 1 + exp(pa[1]) )      ;     lam <- exp(pa[2])
     - sum(w0) * log( a + (1 - a) * exp( - lam ) ) - sum(w1) * log( 1 - a ) - sum(w1 * y1) * pa[2] + sum(w1) * lam 
   }  


   
   

#zip.regs(y, dataset, logged = T, wei = NULL, ncores = 2)

#ela <- function(y, x) {
#  ini <- 2 * as.numeric( logLik( zeroinfl(y~1|1) ) )
#  p <- ncol(x)
#  stat <- numeric(p)
#  for (i in 1:p) { 
#    mod <- zeroinfl(y~x[, i]|1)
#    stat[i] <- as.numeric( logLik(mod) )
#  }
  
#  stat <- 2 * stat - ini
#  pval <- pchisq(stat, 1, lower.tail=F, log.p = T)
#  cbind(stat, pval)
#}
  
  

