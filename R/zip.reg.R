zip.reg <- function(target, dataset, wei = NULL, lgy = NULL) {
  
  n <- length(target)
  x <- model.matrix(target ~ ., as.data.frame(dataset) )[1:n, ]
  poia <- which(target == 0)
  n0 <- length(poia)   ;    n1 <-  n - n0
  target1 <- target[ -poia ]    
  
  if ( is.null(wei) ) {
    mod <- glm.fit( x[ -poia, ], target1, family = poisson(log) ) 
    p1 <- ( n0 - sum( exp( - fitted(mod) ) ) ) / n
    g1 <- log( p1 / ( 1 - p1 ) )
    
    lik <- nlm( regzip, c(g1, coef(mod) ), y1 = target1, x = x, n1 = n1, poia = poia, iterlim = 10000 )
    lik2 <- nlm( regzip, lik$estimate, y1 = target1, x = x, n1 = n1, poia = poia, iterlim = 10000 )
    if ( is.null(lgy) )  lgy <- sum( lgamma(target1 + 1) )  
    
  } else {  
    wei <- wei / sum(wei)
    w0 <- wei[poia]     ;    w1 <- wei[-poia]
    mod <- glm.fit( x[ -poia, ], target1, family = poisson(log), weights = w1 ) 
    p1 <- ( n0 - sum( exp( - fitted(mod) ) ) ) / n
    g1 <- log( p1 / ( 1 - p1 ) )
    
    lik <- nlm( regzipwei, c(g1, coef(mod) ), y1 = target1, x = x, n1 = n1, w1 = w1, w0 = w0, iterlim = 10000 )
    lik2 <- nlm( regzipwei, lik$estimate, y1 = target1, x = x, n1 = n1, w1 = w1, w0 = w0, iterlim = 10000 )  
    if ( is.null(lgy) )  lgy <- sum( w1 * lgamma(target1 + 1) )  
    
  }
  
  prop <- exp(lik2$estimate[1]) / ( 1 + exp(lik2$estimate[1]) )
  
  list(iters = lik$iterations + lik2$iterations, be = lik2$estimate[-1], prop = prop, loglik = -lik2$minimum - lgy)
  
}