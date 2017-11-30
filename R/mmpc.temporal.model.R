# target: the target value
# sesObject: the outcome of the ses
# nisgnat: Number of ypografis and generated models. It could be numeric from 1 to total number of ypografis or "all" for all the 
## ypografis. Default is 1.
mmpc.temporal.model = function(target, dataset, reps = NULL, group, slopes = FALSE, wei = NULL, mmpctemporal.Object, test = NULL) {
  
  if ( sum( is.na(mmpctemporal.Object@selectedVars) ) > 0 ) {
    mod = paste("No associations were found, hence no model is produced.")
    ypografi = NULL
    bic = NULL
  }
  
  if ( any(is.na(dataset) ) ) {
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
  }
  
  if ( is.null(test) ) {  
    ci_test = mmpctemporal.Object@test
    slopes = mmpctemporal.Object@slope
  } else ci_test = test 
  
  rob <- FALSE
  ypografi <- mmpctemporal.Object@selectedVars  
  p <- length(ypografi)
  
  if ( ci_test == "testIndGLMM" ) {
    
    if ( length( unique(target) ) == 2 ) {
      if ( is.null(reps) ) {
        mod = lme4::glmer( target ~ dataset[, ypografi] + (1|group), weights = wei, REML = FALSE , family = binomial ) 
      } else {
        reps = reps 
        if (slopes ) {
          mod = lme4::glmer( target ~ reps + dataset[, ypografi] + (reps|group), weights = wei, REML = FALSE, family = binomial )
        } else  mod = lme4::glmer( target ~ reps + dataset[, ypografi] + (1|group), weights = wei, REML = FALSE, family = binomial ) 
      }
    } else if ( identical(floor(target), target) )  {  
      if ( is.null(reps) ) {
        mod = lme4::glmer( target ~ dataset[, ypografi] + (1|group), weights = wei, REML = FALSE , family = poisson ) 
      } else {
        reps = reps 
        if (slopes ) {
          mod = lme4::glmer( target ~ reps + dataset[, ypografi] + (reps|group), weights = wei, REML = FALSE, family = poisson )
        } else  mod = lme4::glmer( target ~ reps + dataset[, ypografi] + (1|group), weights = wei, REML = FALSE, family = poisson ) 
      }
    } else {
      if ( is.null(reps) ) {
        mod = lme4::lmer( target ~ dataset[, ypografi] + (1|group), weights = wei, REML = FALSE )
      } else {
        reps = reps 
        if ( slopes ) {
          mod = lme4::lmer( target ~ reps + dataset[, ypografi] + (reps|group), weights = wei, REML = FALSE ) 
        } else {
          reps = reps 
          mod = lme4::lmer( target ~ reps + dataset[, ypografi] + (1|group), weights = wei, REML = FALSE )        
        }
      }   
    }
    
    bic <- BIC(mod)
    if ( is.null( colnames(dataset) ) ) {
      names(ypografi) = paste("Var", ypografi, sep = " ")
    } else  names(ypografi) = colnames(dataset)[ypografi]
    ypografi <- c(ypografi, bic)
    names(ypografi)[length(ypografi)] = "bic"
    
  }  ## end if ( ci_test == "testIndGLMM" )
  
  list(mod = mod, ypografi = ypografi)  
}

