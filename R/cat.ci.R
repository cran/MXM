####################
#### G^2 (and X^2) test of (un)conditional independence
####################

cat.ci <- function(xi, yi, cs, dataset, type = type, rob = FALSE, R = 1) {
  ## the xi and yi are two numbers, 1 and 2 for example
  ## indicating the two variables whose conditional independence 
  ## will be tested
  ## xi, yi and cs must be different, non onverlapping numbers
  ## cs is one or more numbers indicating the conditioning variable(s)
  ## it is et to 0 by default. In this case an uncodntional test of 
  ## independence is  performed
  ## dataset is the whole dataset, and is expected to be a matrix
  ## type is set to NULL be default. This argument is not taken into consideration anywhere
  ## rob is FALSE by default, even if it is TRUE it is not taken into cosideration
  ## the type and rob arguments are put here so as to have the same signature as condi
  if ( sum(cs == 0) > 0 ) {  ## There are no conditioning variables
  
    a1 <- Rfast::g2Test_univariate(dataset[, c(xi, yi)], type[c(xi, yi)])
    stat <- as.numeric( a1$statistic )
    dof <- as.numeric( a1$df )
    pval <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    res <- c( as.numeric(stat), pval, dof )  

  } else {   ## There are conditioning variables
    a1 <- Rfast::g2Test(dataset, xi, yi, cs, type)           
    stat <- as.numeric( a1$statistic )
    dof <- as.numeric( a1$df )
    pval <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)

    if ( nrow(dataset) > 5 * dof ) {  ## condition to perform the test
       res <- c( as.numeric(stat), pval, dof )  
    } else  res <- c( 0, 0, 1 )  

  }
  
  names(res) <- c("Chi-squared test", "logged p-value", "df")
  res
}