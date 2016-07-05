################################
#### Permutation based hypothesis testing 
#### for a zero correlation coefficient 
####
################################

permcor <- function(x, R = 999) {
  ## x is a 2 column matrix containing the data
  ## type can be either "pearson" or "spearman"
  ## R is the number of permutations
  
  x <- as.matrix(x)
  n <- nrow(x)
  r <- cor(x)[2]
  test <- 0.5 * log( (1 + r)/(1 - r) )  ## the test statistic
  x1 <- x[, 1]      ;     x2 <- x[, 2]
  m1 <- sum(x1)     ;     m12 <- sum(x1^2)
  m2 <- sum(x2)     ;     m22 <- sum(x2^2)
  up <-  m1 * m2 / n
  down <- sqrt( (m12 - m1^2 / n) * (m22 - m2^2 / n) )

  sxy <- numeric(R)
  for (i in 1:R) {
    y1 <- sample(x1, n)
    sxy[i] <- sum(y1 * x2)
  }  
    rb <- (sxy - up) / down
    tb <- 0.5 * log( (1 + rb)/(1 - rb) )  ## the test statistic

  pvalue <- ( sum( abs(tb) > abs(test) ) + 1 ) / (R + 1)  ## bootstrap p-value
  res <- c( r, pvalue )
  names(res) <- c('correlation', 'p-value')
  res
}
