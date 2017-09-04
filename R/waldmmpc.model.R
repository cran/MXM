# target: the target value
# sesObject: the outcome of the ses
# nisgnat: Number of ypografis and generated models. It could be numeric from 1 to total number of ypografis or "all" for all the 
## ypografis. Default is 1.
waldmmpc.model = function(target, dataset, wei = NULL, wald.mmpcObject, test = NULL) {
  
  if ( sum( is.na(wald.mmpcObject@selectedVars) ) > 0 ) {
    mod = paste("No associations were found, hence no model is produced.")
    ypografi = NULL
    bic = NULL
  }
  
  if ( any(is.na(dataset) ) ) {
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    dataset = apply(dataset, 2, function(x) { x[ which(is.na(x)) ] = median(x,na.rm = TRUE) } );
  }
  
  if ( is.null(test) ) {  
    ci_test = wald.mmpcObject@test
  } else ci_test = test 
  
  rob <- FALSE
  ypografi <- wald.mmpcObject@selectedVars  
  p <- length(ypografi)

  if ( ci_test == "waldMMreg" ) {
    if ( min(target) > 0 & max(target) < 1 )   target <- log( target/(1 - target) ) 
      mod = MASS::rlm(target ~., data = as.data.frame(dataset[, ypografi ]), maxit = 2000, weights = wei, method = "MM" )
      bic = BIC( mod )        

  } else if ( ci_test == "waldBeta" ) {
    mod <- beta.mod( target, dataset[, ypografi ], wei = wei )
    bic <-  - 2 * mod$loglik + ( length( coef(mod$be) ) + 1 ) * log( length(target) )
    
  } else if ( ci_test == "waldSpeedglm" ) {
    la <- length( unique(target) )
    if ( la == 2 ) {
      mod <- speedglm::speedglm( target ~ ., data = as.data.frame(dataset[, ypografi ]), weights = wei, family = binomial() )
      bic <-  - 2 * as.numeric( logLik(mod) ) + length( coef(mod) ) * log( length(target) )
      
    } else if ( la  > 2  &  sum( floor(target) - target) == 0 ) {
      mod <- speedglm::speedglm( target ~ ., data = as.data.frame(dataset[, ypografi ]), weights = wei, family = poisson() )
      bic <-  - 2 * as.numeric( logLik(mod) ) + length( coef(mod) ) * log( length(target) )
      
    } else {
      mod <- speedglm::speedlm( target ~ ., data = as.data.frame(dataset[, ypografi ]), weights = wei )
      bic <-  - 2 * as.numeric( logLik(mod) ) + ( length( coef(mod) ) + 1 ) * log( length(target) )
    }  
    
  } else if ( ci_test == "waldPois") {
    mod <- glm( target ~ ., data = as.data.frame(dataset[, ypografi ]), family = poisson, weights = wei )
    bic <- BIC( mod )
    
  } else if ( ci_test == "waldNB" ) {
    mod <- MASS::glm.nb( target ~ ., data = as.data.frame(dataset[, ypografi ]), weights = wei )
    bic <- BIC(mod)
    
  } else if ( ci_test == "waldZIP" ) {
    mod <- zip.mod( target, dataset[, ypografi], wei = wei )
    bic <-  -2 * mod$loglik + ( p + 1) * log( length(target) )
    
  } else if ( ci_test == "waldBinom" ) {
    mod <- glm( target[, 1] /target[, 2] ~., data = as.data.frame(dataset[, ypografi ]), weights = target[, 2], family = binomial )
    bic <- BIC(mod)
    
  } else if ( ci_test == "waldGamma" ) {
    mod <- glm( target ~ ., data = as.data.frame(dataset[, ypografi ]), weights = wei, family = Gamma(link = log) )
    bic <- BIC(mod)
    
  } else if ( ci_test == "waldNormLog" ) {
    mod <- glm( target ~ ., data = as.data.frame(dataset[, ypografi ]), weights = wei, family = gaussian(link = log) )
    bic <- BIC(mod)
    
  } else if ( ci_test == "waldTobit" ) {
    mod <- survival::survreg( target ~ ., data = as.data.frame(dataset[, ypografi ]), weights = wei, dist = "gaussian" )
    bic <-  - 2 * as.numeric( logLik(mod) ) + ( p + 1 ) * log( dim(dataset)[1] )
    
  } else if (ci_test == "waldCR") {
    mod <- survival::coxph( target ~ ., data = as.data.frame(dataset[, ypografi ]), weights = wei )
    bic <- BIC(mod)
    
  } else if (ci_test == "waldWR") {
    mod <- survival::survreg( target ~ ., data = as.data.frame(dataset[, ypografi ]), weights = wei )
    bic <- BIC(mod)
    
  } else if (ci_test == "waldER") {
    mod <- survival::survreg( target ~ ., data = as.data.frame(dataset[, ypografi ]), weights = wei, distribution = "exponential" )
    bic <- BIC(mod)

  } else if ( ci_test == "waldBinary" ) {
    mod <- glm( target ~., data = as.data.frame(dataset[, ypografi ]), family = binomial, weights = wei ) 
    bic <- BIC(mod)
    
  } else if ( ci_test == "waldOrdinal" ) {
    mod <- ordinal::clm( target ~., data = as.data.frame(dataset[, ypografi ]), weights = wei ) 
    bic <- BIC(mod)
    
  } else if (ci_test == "waldIGreg") {
    mod <- glm( target ~ ., data = as.data.frame(dataset[, ypografi ]), family = inverse.gaussian(log), weights = wei )
    bic <- BIC(mod)
    
  }
  
  if ( is.null( colnames(dataset) ) ) {
    names(ypografi) = paste("Var", ypografi, sep = " ")
  } else  names(ypografi) = colnames(dataset)[ypografi]
  
  ypografi <- c(ypografi, bic)
  names(ypografi)[length(ypografi)] = "bic"
  
  list(mod = mod, ypografi = ypografi)  
  
}

