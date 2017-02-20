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
    #dataset = as.matrix(dataset);
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    if (class(dataset) == "matrix")   {
      dataset = apply(dataset, 2, function(x) { x[ which(is.na(x)) ] = median(x,na.rm = TRUE) } );
    }else{
      for(i in 1:ncol(dataset))
      {
        if( any( is.na(dataset[,i]) ) )
        {
          xi = dataset[, i]
          if(class(xi) == "numeric")
          {                    
            xi[ which( is.na(xi) ) ] = median(xi, na.rm = TRUE) 
          }else if ( class(xi) == "factor" ) {
            xi[ which( is.na(xi) ) ] = levels(xi)[ which.max(xi) ]
          }
          dataset[, i] = xi
        }
      }
    }
  }
  
  if ( is.null(test) ) {  
    ci_test = mmpcObject@test
  } else ci_test = test 
  
    rob <- mmpcObject@rob
    ypografi <- mmpcObject@selectedVars  
    p <- length(ypografi)
    # mat1 <- mat2 <- numeric(p)
    
    if ( ci_test == "testIndFisher" || ci_test == "testIndReg" ) {
      if ( min(target) > 0 & max(target) < 1 )   target <- log( target/(1 - target) ) 
      
      if ( rob ) {
        mod = MASS::rlm(target ~., data = as.data.frame(dataset[, ypografi ]), maxit = 2000, weights = wei )
        bic = BIC( mod )        
      } else {
        mod = lm( target ~ ., data = as.data.frame(dataset[, ypografi ]), weights = wei )
        bic = BIC(mod)      
      }  
      
    } else if ( ci_test == "testIndSpearman" || ci_test == "testIndRQ" ) {
      if ( all( target>0 & target<1 ) )  target = log( target/(1 - target) ) 

      mod <- quantreg::rq( target ~., data = as.data.frame(dataset[, ypografi ]), weights = wei  )
  	  la <- logLik(mod)
      bic <-  - 2 * as.numeric( la ) + attr(la, "df") * log( length(target) )   
      # for (i in 1:p) {
      #   mi <- quantreg::rq( target ~ ., data = as.data.frame( dataset[, c(1:i) ] ) )
      #   es <- fitted(mi)
      #   mat1[i] <- cor(target, es)^2
      #   mi <- quantreg::rq( target ~ ., data = as.data.frame( dataset[, -i ] ) )
      #   es <- fitted(mi)
      #   mat2[i] <- cor(target, es)^2
      # }  
      
    } else if ( ci_test == "testIndBeta" ) {
      mod <- beta.mod( target, dataset[, ypografi ], wei = wei )
      bic <-  - 2 * mod$loglik + ( length( coef(mod$be) ) + 1 ) * log( length(target) )
      
    } else if ( ci_test == "testIndSpeedglm" ) {
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
      
    } else if ( ci_test == "testIndPois") {
      #if ( rob == TRUE ) {
      #  mod <- robust::glmRob( target ~ ., data = as.data.frame(dataset[, ypografi ]) , poisson, maxit = 100 )
      #  bic <- mod$deviance + length( coef(mod) ) * log( length(target) )
      #} else {
      mod <- glm( target ~ ., data = as.data.frame(dataset[, ypografi ]), family = poisson, weights = wei )
      bic <- BIC( mod )
      #}
      
    } else if ( ci_test == "testIndNB" ) {
      mod <- MASS::glm.nb( target ~ ., data = as.data.frame(dataset[, ypografi ]), weights = wei )
      bic <- BIC(mod)
      
    } else if ( ci_test == "testIndZIP" ) {
      mod <- zip.mod( target, dataset[, ypografi], wei = wei )
      bic <-  -2 * mod$loglik + ( length( coef(mod$be) ) + 1) * log( length(target) )
      
    } else if ( class(target) == "matrix" & ci_test == "testIndMVreg" ) {
      if ( all(target > 0 & target < 1) )    target <- log( target[, -1] / (target[, 1]) ) 
      mod <- lm( target ~., data = as.data.frame(dataset[, ypografi ]), weights = wei )
      bic <- NULL
      
    } else if ( class(target) == "matrix"  &  ci_test == "testIndBinom" ) {
      mod <- glm( target[, 1] /target[, 2] ~., data = as.data.frame(dataset[, ypografi ]), weights = target[, 2], family = binomial )
      bic <- BIC(mod)
      
    } else if (ci_test == "censIndCR") {
      mod <- survival::coxph( target ~ ., data = as.data.frame(dataset[, ypografi ]), weights = wei )
      bic <- BIC(mod)
      
    } else if (ci_test == "censIndWR") {
      mod <- survival::survreg( target ~ ., data = as.data.frame(dataset[, ypografi ]), weights = wei )
      bic <- BIC(mod)

    } else if (ci_test == "censIndER") {
      mod <- survival::survreg( target ~ ., data = as.data.frame(dataset[, ypografi ]), weights = wei, distribution = "exponential" )
      bic <- BIC(mod)
      
    } else if (ci_test == "testIndClogit") {
      case <- as.logical(target[, 1]);  
      id <- target[, 2]
      mod <- survival::clogit(case ~ . + strata(id), data = as.data.frame( dataset[ , ypografi] ) )  ## wieghts are ignored here anyway
      bic <- BIC(mod)
      
    } else if ( ci_test == "testIndLogistic" || ci_test == "gSquare" ) {
      if ( length( unique(target)) == 2 ) {
        #if ( rob == TRUE ) {
        #  mod <- robust::glmRob( target ~ ., data = as.data.frame(dataset[, ypografi ]) , binomial, maxit = 100 )
        #  bic <- mod$deviance + length( coef(mod) ) * log( length(target) )
        #} else { 
        mod <- glm( target ~., data = as.data.frame(dataset[, ypografi ]), family = binomial, weights = wei ) 
        bic <- BIC(mod)
        #}
        
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

