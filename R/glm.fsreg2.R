glm.fsreg2 <- function(target, dataset, ini = NULL, threshold = 0.05, wei = NULL, tol = 2, heavy = FALSE, robust = FALSE, ncores = 1) {
  
  ## target can be Real valued (normal), binary (binomial) or counts (poisson)

  if ( !is.null(ini) ) {
    result <- glm.fsreg_2(target, dataset, iniset = ini, wei = wei, threshold = threshold, tol = tol, heavy = heavy, robust = FALSE, ncores = ncores) 
    
  } else {  ## else do the classical forward regression
    result <- glm.fsreg(target, dataset, wei = wei, threshold = threshold, tol = tol, heavy = heavy, robust = FALSE, ncores = ncores) 
  }  
    
  result
}    










