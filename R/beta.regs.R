beta.regs <- function(target, dataset, wei = NULL, logged = FALSE, ncores = 1) {

  D <- ncol(dataset)

  if ( ncores <= 1 ) {
  
    if ( is.null(wei) ) {
      ly1 <- log(1 - target)     ;    ly <- log(target) - ly1           
      sly1 <- sum(ly1)           ;    sly2 <- sum( log( target ) ) + sly1  
      a <- Rfast::beta.mle(target)
      #ini <- 2 * a$loglik
      ini <- 2 * a$loglik
    } else {
      w <- wei / sum(wei)
      ly1 <- w * log(1 - target)     ;    ly <- w * log(target) - ly1
      sly1 <- sum(ly1)               ;    sly2 <- sum( w * log( target * (1 - target) ) ) 
      a <- betamle.wei(target, w)    
      ini <- 2 * a$loglik
    }
    
    n <- length(target)
    iniphi <- log( sum(a$param) )
    loglik <- dof <- numeric(D)

    if ( is.null(wei) ) {

      for (i in 1:D) {
        x <- model.matrix(target ~ ., as.data.frame(dataset[, i]) )
        options(warn = -1)
        mod1 <- nlm(regbeta, c(iniphi, numeric(dim(x)[2]) ), 
             ly = ly, sly1 = sly1, x = x, n = n, iterlim = 10000 )
        mod2 <- nlm(regbeta, mod1$estimate, ly = ly, sly1 = sly1, x = x, n = n, 
             iterlim = 10000 )
        loglik[i] <-  - mod2$minimum 
        dof[i] <- dim(x)[2] - 1
      }

    } else {

      for (i in 1:D) {
        x <- model.matrix(target ~ ., as.data.frame(dataset[, i]) )
        options(warn = -1)
        mod1 <- nlm(regbetawei, c(iniphi, numeric(dim(x)[2]) ), 
             ly = ly, sly1 = sly1, x = x, w = w, iterlim = 10000 )
        mod2 <- nlm(regbetawei, mod1$estimate, ly = ly, sly1 = sly1, x = x, w = w, 
             iterlim = 10000 )
        loglik[i] <-  - mod2$minimum 
        dof[i] <- dim(x)[2] - 1
      }

    }
    
    lik <- loglik - sly2
    bic <-  - 2 * lik + dof * log(n) 
    stat <- 2 * lik - ini
    pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = logged) 
    
  } else {

   if ( is.null(wei) ) {

     cl <- makePSOCKcluster(ncores)
     registerDoParallel(cl)

     ly1 <- log(1 - target)     ;    ly <- log(target) - ly1           
     sly1 <- sum(ly1)           ;    sly2 <- sum( log( target ) ) + sly1  
     a <- Rfast::beta.mle(target)
     ini <- 2 * a$loglik
     n <- length(target)
     iniphi <- log( sum(a$param) )

     mod <- foreach( i = 1:D, .combine = rbind, .export = "regbeta" ) %dopar% {

       x <- model.matrix(target ~ ., as.data.frame(dataset[, i]) )
       options(warn = -1)
       mod1 <- nlm(regbeta, c(iniphi, numeric(dim(x)[2]) ), 
            ly = ly, sly1 = sly1, x = x, n = n, iterlim = 10000 )
       mod2 <- nlm(regbeta, mod1$estimate, ly = ly, sly1 = sly1, x = x, n = n, 
            iterlim = 10000 )
       return( c( - mod2$minimum, dim(x)[2] - 1 ) )
     }
    
     stopCluster(cl)

   } else {

     cl <- makePSOCKcluster(ncores)
     registerDoParallel(cl)
     
     n <- length(target)
     w <- wei / sum(wei)
     ly1 <- w * log(1 - target)     ;    ly <- w * log(target) - ly1
     sly1 <- sum(ly1)               ;    sly2 <- sum( w * log( target * (1 - target) ) )    
     a <- betamle.wei(target ,w)
     ini <- 2 * a$loglik
     iniphi <- log( sum( a$param ) )

     mod <- foreach( i = 1:D, .combine = rbind, .export = "regbetawei" ) %dopar% {

       x <- model.matrix(target ~ ., as.data.frame(dataset[, i]) )
       options(warn = -1)
        mod1 <- nlm(regbetawei, c(iniphi, numeric(dim(x)[2]) ), 
             ly = ly, sly1 = sly1, x = x, w = w, iterlim = 10000 )
        mod2 <- nlm(regbetawei, mod1$estimate, ly = ly, sly1 = sly1, x = x, w = w, 
             iterlim = 10000 )
       return( c( - mod2$minimum, dim(x)[2] - 1 ) )
     }
     stopCluster(cl)
   }
    
    lik <- mod[, 1] - sly2
    bic <-  - 2 * lik + mod[, 2] * log(n)
    stat <- 2 * lik - ini
    pvalue <- pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = logged) 
  }  

  cbind(stat, pvalue, bic)
}






## objective functions used in the optimisation procedure

regbeta <- function(pa, ly, sly1, x, n) {
  phi <- exp(pa[1])    ;    b <- pa[-1]
  m <- exp( tcrossprod(b, x) )
  m <- m / ( 1 + m )
  a1 <- m * phi   ;   a2 <- phi - a1
  - n * lgamma(phi) + sum( lgamma(a1) ) + sum( lgamma(a2) ) - sum(a1 * ly) - phi * sly1
}


regbetawei <- function(pa, ly, sly1, x, w) {
  phi <- exp(pa[1])    ;    b <- pa[-1]
  m <- exp( tcrossprod(b, x) )
  m <- m / ( 1 + m )
  a1 <- m * phi   ;   a2 <- phi - a1
  - lgamma(phi) + sum( w * lgamma(a1) ) + sum( w * lgamma(a2) ) - sum(a1 * ly) - phi * sly1
} 


### log-likelihood with weights but no covariates

betamle.wei <- function(y, wei) {
  
  n <- length(y)
  w <- wei / sum(wei)
  ly1 <- sum( w * log(y) )      ;      ly2 <- sum( w * log(1 - y) )  
  
  betawei <- function(pa, ly1, ly2) {
    a <- exp(pa[1])     ;     b <- exp(pa[2])
    lbeta(a, b) - (a - 1) * ly1 - (b - 1) * ly2
  } 
  
  iniphi <- sum( y * (1 - y) ) / var(y) / n 
  a1 <- sum(y) * iniphi / n        ;        a2 <- iniphi - a1
  options(warn = -1)
  lik <- nlm( betawei, c( log(a1), log(a2) ), ly1 = ly1, ly2 = ly2, iterlim = 10000 )
  lik2 <- nlm( betawei, lik$estimate, ly1 = ly1, ly2 = ly2, iterlim = 10000 )  
  
  list(iters = lik$iterations + lik2$iterations, param = exp(lik2$estimate), loglik = -lik2$minimum )
}



#ela = function(y, x) {
#  ini = as.numeric( logLik( betareg(y ~ 1) ) )
#  D <- dim(x)[2]
#  lik = dof = numeric(D)
#  for (i in 1:D) {
#    mod <- betareg(y ~ x[, i])
#    lik[i] = as.numeric( logLik(mod) )
#    dof[i] = length( coef(mod) )
#  } 
#  stat <- 2 * (lik - ini)
#  pvalue <- pchisq(stat, dof - 2, lower.tail =F)
#  cbind(stat, pvalue)
#}