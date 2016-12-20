beta.mod <- function(target, dataset, wei = NULL, xnew = NULL) {

   n <- length(target)
   x <- model.matrix(target ~ ., data.frame(dataset) )
   iniphi <- log( sum( target * (1 - target) ) / var(target) / n )
  
   if ( is.null(wei) ) {
    ly1 <- log(1 - target)     ;    ly <- log(target) - ly1           
    sly1 <- sum(ly1)           ;    sly2 <- sum( log(target) ) + sly1   

    options(warn = -1)
    mod1 <- nlm(regbeta, c( iniphi, numeric(dim(x)[2]) ), 
            ly = ly, sly1 = sly1, x = x, n = n, iterlim = 10000 )
    mod2 <- optim(mod1$estimate, regbeta, ly = ly, sly1 = sly1, x = x, n = n, 
            control = list(maxit = 10000), hessian = TRUE  )
  } else {
    w <- wei / sum(wei)
    ly1 <- w * log(1 - target)     ;    ly <- w * log(target) - ly1
    sly1 <- sum(ly1)    ;    sly2 <- sum( w * log( target ) ) + sly1   

    options(warn = -1)
    mod1 <- nlm(regbetawei, c(iniphi, numeric(dim(x)[2]) ), 
            ly = ly, sly1 = sly1, x = x, w = w, iterlim = 10000 )
    mod2 <- optim(mod1$estimate, regbetawei, ly = ly, sly1 = sly1, x = x, w = w, 
            control = list(maxit = 10000), hessian = TRUE )
    
  }
  
  seb <- sqrt( diag( solve(mod2$hessian) ) )
  be <- cbind(mod2$par[-1], seb[-1] ) 
  be <- cbind(be, (be[, 1] / be[, 2] )^2 )   
  pv <- pchisq(be[, 3], 1, lower.tail = FALSE)
  be <- cbind(be, pv)
  rownames(be) <- colnames(x)
  colnames(be) <- c("coefficient", "std error", "chisq-statistic", "p-value")
  
  if ( !is.null(xnew) ) {
    xnew <- model.matrix(~., data.frame(xnew) )
    est <- exp( as.vector( xnew %*% be[, 1] ) )
  } else  est <- exp( as.vector( x %*% be[, 1]) )
    est <- est / (1 + est)
  
   list(be = be, phi = exp(mod2$par[1]), loglik = - mod2$value - sly2, est = est)      
}

  