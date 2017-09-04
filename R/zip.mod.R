zip.mod <- function(target, dataset, wei = NULL, xnew = NULL) {

  n <- length(target)
  x <- model.matrix(target ~ ., data.frame(dataset) )
  poia <- which(target == 0)
  n0 <- length(poia)   ;    n1 <-  n - n0
  target1 <- target[ -poia ]    
  
  if ( is.null(wei) ) {
    mod <- glm.fit( x[ -poia, ], target1, family = poisson(log) ) 
    p1 <- ( n0 - sum( exp( - fitted(mod) ) ) ) / n
    g1 <- log( p1 / ( 1 - p1 ) )
    lik <- nlm( regzip, c(g1, coef(mod) ), y1 = target1, x = x, n1 = n1, poia = poia, iterlim = 10000 )
    lik2 <- optim( lik$estimate, regzip, y1 = target1, x = x, n1 = n1, poia = poia, control = list(maxit = 10000), hessian = TRUE )
    lgy <- sum( lgamma(target1 + 1) )  
    
  } else {  
    wei <- wei / sum(wei)
    w0 <- wei[poia]     ;    w1 <- wei[-poia]
    mod <- glm.fit( x[ -poia, ], target1, family = poisson(log), weights = w1 ) 
    p1 <- ( n0 - sum( exp( - fitted(mod) ) ) ) / n
    g1 <- log( p1 / ( 1 - p1 ) )
    lik <- nlm( regzipwei, c(g1, coef(mod) ), y1 = target1, x = x, n1 = n1, w1 = w1, w0 = w0, poia = poia, iterlim = 10000 )
    lik2 <- optim( lik$estimate, regzipwei, y1 = target1, x = x, n1 = n1, w1 = w1, w0 = w0, poia = poia, control = list(maxit = 10000), hessian = TRUE )  
    lgy <- sum( w1 * lgamma(target1 + 1) )  
  }
  
  prop <- exp(lik2$par[1]) / ( 1 + exp(lik2$par[1]) )
    
  seb <- sqrt( diag( solve(lik2$hessian) ) )
  be <- cbind(lik2$par[-1], seb[-1] ) 
  be <- cbind(be, (be[, 1] / be[, 2] )^2 )   
  pv <- pchisq(be[, 3], 1, lower.tail = FALSE)
  be <- cbind(be, pv)
  rownames(be) <- colnames(x)
  colnames(be) <- c("coefficient", "std error", "chisq-statistic", "p-value")
  
  if ( !is.null(xnew) ) {
    xnew <- model.matrix(~., data.frame(xnew) )
    est <- exp( as.vector( xnew %*% be[, 1] ) )
  } else  est <- exp( as.vector( x %*% be[, 1]) )
      
   list(be = be, prop = prop, loglik = - lik2$value - lgy, est = est)      
}

  