auc <- function(group, preds, roc = FALSE, cutoffs = NULL) {
  
  group <- as.numeric( as.factor(group) ) - 1
  ri <- rank(preds)
  n <- length(preds)
  n1 <- sum( group )
  n0 <- n - n1
  s1 <- sum( ri[ group == 1 ] )
  auc <- ( s1 - 0.5 * n1 * (n1 + 1) ) / n0 / n1
  result <- auc
  
  if ( roc ) {
    
    if ( is.null(cutoffs) ) {
      lena <- seq(1, 0, by = -0.01)
    } else  lena <- cutoffs
    
    nu <- length(lena)
    sens <- spec <- numeric(nu)
    npos <- sum( group == 1 )
    nnegs <- sum( group == 0 )
    
    for ( i in 1:nu ) {
      estp <- as.numeric( preds >= lena[i] )      
      sens[i] <-  sum( group == 1  &  estp == 1 ) / npos
      spec[i] <-  sum( group == 0  &  estp == 0 ) / nnegs
    }
    
    sens[1] <- 0
    spec[1] <- 1
    sens[nu + 1] <- 1
    spec[nu + 1] <- 0
    cutoffs <- c(1, lena, 0)    

    plot( 1 - spec, sens, type = "l", lty = 5, lwd = 2, xlab = "1- specificity", ylab = "Sensitivity", col = 3,
          xlim = c(0, 1))
    abline(a = 0, b = 1, lty = 2, col = 2)
    qa <- which.max( sens + spec - 1 )
    points( 1 - spec[qa], sens[qa], lwd = 2)
    youden <- c( 1 - spec[qa], sens[qa] )
    names(youden) <- c("1-specificity", "sensitivity")
    result <- list(cutoffs = cutoffs, sensitivity = sens, specificity = spec, youden = youden, auc = auc)
  }
  
  result
  
}