# target: the target value
# sesObject: the outcome of the ses
# nisgnat: Number of ypografis and generated models. It could be numeric from 1 to total number of ypografis or "all" for all the 
## ypografis. Default is 1.
mmpc.model = function(target, dataset, wei = NULL, mmpcObject, test = NULL) {
  
  if ( sum( is.na(mmpcObject@selectedVars) ) > 0 ) {
    mod = paste("No associations were found, hence no model is produced.")
    ypografi = NULL
    bic = NULL
  }
  
  if ( any(is.na(dataset) ) ) {
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    if ( is.matrix(dataset) )  {
      dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
    } else {
      poia <- which( is.na(dataset), arr.ind = TRUE )[2]
      for ( i in poia )  {
        xi <- dataset[, i]
        if ( is.numeric(xi) ) {                    
          xi[ which( is.na(xi) ) ] <- median(xi, na.rm = TRUE) 
        } else if ( is.factor( xi ) ) {
          xi[ which( is.na(xi) ) ] <- levels(xi)[ which.max( as.vector( table(xi) ) )]
        }
        dataset[, i] <- xi
      }
    }
  }
  
  if ( is.null(test) ) {  
    ci_test = mmpcObject@test
  } else ci_test = test 
  
    rob <- mmpcObject@rob
    ypografi <- mmpcObject@selectedVars  
    p <- length(ypografi)

    if ( ci_test == "testIndFisher" || ci_test == "testIndReg" ) {
      if ( min(target) > 0 & max(target) < 1 )   target <- log( target/(1 - target) ) 
      
      if ( rob ) {
        mod = MASS::rlm(target ~., data = data.frame(dataset[, ypografi ]), maxit = 2000, weights = wei )
        bic = BIC( mod )        
      } else {
        mod = lm( target ~ ., data = as.data.frame(dataset[, ypografi ]), weights = wei )
        bic = BIC(mod)      
      }  
      
    } else if ( ci_test == "testIndSpearman" || ci_test == "testIndRQ" ) {

      mod <- quantreg::rq( target ~., data = as.data.frame(dataset[, ypografi ]), weights = wei  )
  	  la <- logLik(mod)
      bic <-  - 2 * as.numeric( la ) + attr(la, "df") * log( length(target) )   

    } else if ( ci_test == "testIndBeta" ) {
      mod <- beta.mod( target, dataset[, ypografi ], wei = wei )
      bic <-  - 2 * mod$loglik + ( length( mod$be ) + 1 ) * log( length(target) )
      
    } else if ( ci_test == "testIndSpeedglm" ) {
      la <- length( unique(target) )
      if ( la == 2 ) {
        mod <- speedglm::speedglm( target ~ ., data = as.data.frame(dataset[, ypografi ]), weights = wei, family = binomial() )
        bic <-  - 2 * as.numeric( logLik(mod) ) + length( coef(mod) ) * log( length(target) )
      
      } else if ( la  > 2  &  sum( floor(target) - target) == 0 ) {
        mod <- speedglm::speedglm( target ~ ., data = as.data.frame(dataset[, ypografi ]), weights = wei, family = poisson() )
        bic <-  - 2 * as.numeric( logLik(mod) ) + length( coef(mod) ) * log( length(target) )
		
      } else {
        mod <- speedglm::speedlm( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei )
        bic <-  - 2 * as.numeric( logLik(mod) ) + ( length( coef(mod) ) + 1 ) * log( length(target) )
      }  
      
    } else if ( ci_test == "testIndPois") {
      mod <- glm( target ~ ., data = data.frame(dataset[, ypografi ]), family = poisson, weights = wei )
      bic <- BIC( mod )

    } else if ( ci_test == "testIndNB" ) {
      mod <- MASS::glm.nb( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei )
      bic <- BIC(mod)
      
    } else if ( ci_test == "testIndZIP" ) {
      mod <- zip.mod( target, dataset[, ypografi], wei = wei )
      bic <-  -2 * mod$loglik + ( length( coef(mod$be) ) + 1) * log( length(target) )
      
    } else if (ci_test == "testIndIGreg") {
      mod <- glm(target ~., data = data.frame( dataset[, ypografi] ), family = inverse.gaussian(log), weights = wei)
      bic <- BIC(mod)
      
    } else if ( is.matrix(target)  & ci_test == "testIndMVreg" ) {
      if ( all(target > 0 & target < 1)  &  Rfast::Var( Rfast::rowsums(target) ) == 0 )   target = log( target[, -1]/(target[, 1]) ) 
      mod <- lm( target ~., data = data.frame(dataset[, ypografi ]), weights = wei )
      bic <- NULL
      
    } else if ( is.matrix(target)   &  ci_test == "testIndBinom" ) {
      mod <- glm( target[, 1] /target[, 2] ~., data = data.frame(dataset[, ypografi ]), weights = target[, 2], family = binomial )
      bic <- BIC(mod)
      
    } else if ( ci_test == "testIndGamma" ) {
      mod <- glm( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei, family = Gamma(link = log) )
      bic <- BIC(mod)
      
    } else if ( ci_test == "testIndNormLog" ) {
      mod <- glm( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei, family = gaussian(link = log) )
      bic <- BIC(mod)
      
    } else if ( ci_test == "testIndTobit" ) {
      mod <- survival::survreg( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei, dist = "gaussian" )
      bic <-  - 2 * as.numeric( logLik(mod) ) + ( length( mod$coefficients ) + 1 ) * log( NROW(dataset) )
      
    } else if (ci_test == "censIndCR") {
      mod <- survival::coxph( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei )
      bic <- BIC(mod)
      
    } else if (ci_test == "censIndWR") {
      mod <- survival::survreg( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei )
      bic <-  - 2 * as.numeric( logLik(mod) ) + ( length( mod$coefficients ) + 1 ) * log( NROW(dataset) )
      
    } else if (ci_test == "censIndER") {
      mod <- survival::survreg( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei, distribution = "exponential" )
      bic <-  - 2 * as.numeric( logLik(mod) ) + ( length( mod$coefficients ) + 1 ) * log( NROW(dataset) )
      
    } else if (ci_test == "testIndClogit") {
      case <- as.logical(target[, 1]);  
      id <- target[, 2]
      mod <- survival::clogit(case ~ . + strata(id), data = data.frame( dataset[ , ypografi] ) )  ## wieghts are ignored here anyway
      bic <- BIC(mod)
      
    } else if ( ci_test == "testIndLogistic" || ci_test == "gSquare" ) {
      if ( length( unique(target)) == 2 ) {
        mod <- glm( target ~., data = as.data.frame(dataset[, ypografi ]), family = binomial, weights = wei ) 
        bic <- BIC(mod)

      } else if ( !is.ordered(target) ) { 
        target <- as.factor( as.numeric( as.vector(target) ) )
        mod <- nnet::multinom( target ~., data = as.data.frame(dataset[, ypografi ]), trace = FALSE, weights = wei )
        bic <- BIC(mod)
        
      } else if ( is.ordered(target) ) {
        mod <- ordinal::clm( target ~., data = as.data.frame(dataset[, ypografi ]), weights = wei  )
        bic <- BIC(mod)
      }    
      
    }
       
    if ( is.null( colnames(dataset) ) ) {
      names(ypografi) = paste("Var", ypografi, sep = " ")
    } else  names(ypografi) = colnames(dataset)[ypografi]
    
    ypografi <- c(ypografi, bic)
    names(ypografi)[length(ypografi)] = "bic"
    
    list(mod = mod, ypografi = ypografi)  
  
}

