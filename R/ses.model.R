# target: the target value
# sesObject: the outcome of the ses
# nisgnat: Number of ypografis and generated models. It could be numeric from 1 to total number of ypografis or "all" for all the 
## ypografis. Default is 1.

ses.model = function(target, dataset, wei = NULL, sesObject, nsignat = 1, test = NULL) {
  
  if ( sum( is.na(sesObject@signatures) ) > 0 )  mod = paste("No associations were found, hence no model is produced.")
  
  if ( any(is.na(dataset) ) )
  {
    #dataset = as.matrix(dataset);
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    if (class(dataset) == "matrix")
    {
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
    ci_test = sesObject@test
  } else ci_test = test 

  rob = sesObject@rob
  
  if ( nsignat == 1 || ( nsignat > 1 & nrow(sesObject@signatures) == 1 ) ) {
    ypografi = sesObject@selectedVars  
    
    p <- length(ypografi)
    # mat1 <- mat2 <- numeric(p)
    
   if ( ci_test == "testIndFisher"  ||  ci_test == "testIndReg" ) {
     if ( min(target) > 0 & max(target) < 1 )  target = log( target/(1 - target) ) 

     if ( rob ) {
       # cont = robust::lmRob.control(mxr = 2000, mxf = 2000, mxs = 2000 )  ## only used in robust linear regression
       # mod = robust::lmRob( target ~ ., data = as.data.frame( dataset[, ypografi ] ), control = cont )
       mod = MASS::rlm(target ~., data = data.frame(dataset[, ypografi ]), maxit = 2000, weights = wei )
       bic = BIC( mod )
      
     } else {
       mod = lm( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei )
       bic = BIC(mod)     
     }  
     
   } else if ( ci_test == "testIndSpearman"  ||  ci_test == "testIndRQ" ) {
       if ( all( target>0 & target<1 ) )  target = log( target/(1 - target) ) 
      
      mod <- quantreg::rq( target ~., data = data.frame(dataset[, ypografi ]), weights = wei )
	    la <- logLik(mod)
      bic <-  - 2 * as.numeric( la ) + attr(la, "df") * log( length(target) )
      
    } else if ( ci_test == "testIndBeta" ) {     
      mod <- beta.mod( target, dataset[, ypografi ], wei = wei )
      bic <- - 2 * mod$loglik + ( length( coef(mod$be) ) + 1 ) * log( length(target) )
      
    } else if ( ci_test == "testIndSpeedglm" ) {
      la <- unique(target) 
      if ( length( la )  == 2 ) {
	      mod = speedglm::speedglm( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei, family = binomial() )
        bic =  - 2 * as.numeric( logLik(mod) ) + length( coef(mod) ) * log( length(target) )
		
      } else if ( length( la )  > 2  &  sum( floor(target) - target) == 0 ) {
	      mod = speedglm::speedglm( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei, family = poisson() )
        bic =  - 2 * as.numeric( logLik(mod) ) + length( coef(mod) ) * log( length(target) )	  
	  		
      } else {
        mod = speedglm::speedlm( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei )
        bic =  - 2 * as.numeric( logLik(mod) ) + ( length( coef(mod) ) + 1 ) * log( length(target) )
      }  
      
    } else if ( ci_test == "testIndPois") {
      #if ( rob == TRUE ) {
      #  mod <- robust::glmRob( target ~ ., data = as.data.frame(dataset[, ypografi ]) , poisson, maxit = 100 )
      #  bic <- mod$deviance + length( coef(mod) ) * log( length(target) )
      #} else {
        mod <- glm( target ~ ., data = data.frame(dataset[, ypografi ]) , family = poisson, weights = wei )
        bic <- BIC( mod )
      #}

    } else if ( ci_test == "testIndNB" ) {
      mod = MASS::glm.nb( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei )
      bic = BIC(mod)
      
    } else if ( ci_test == "testIndZIP" ) {
      mod <- zip.mod( target, dataset[, ypografi], wei = wei )
      bic <-  -2 * mod$loglik + ( length( coef(mod$be) ) + 1) * log( length(target) )
      ## bic = BIC(mod)
      
    } else if ( class(target) == "matrix"  &  ci_test == "testIndMVreg" ) {
      if ( all(target > 0 & target < 1)  &  sd( Rfast::rowsums(target) ) == 0 )   target = log( target[, -1]/(target[, 1]) ) 
       mod = lm( target ~., data = data.frame(dataset[, ypografi ]), weights = wei )
       bic = NULL
      
    } else if ( class(target) == "matrix"  &  ci_test == "testIndBinom" ) {
      mod = glm( target[, 1] /target[, 2] ~., data = data.frame(dataset[, ypografi ]), weights = target[, 2], family = binomial )
      bic = BIC(mod)
      
    } else if (ci_test == "censIndCR") {
      mod = survival::coxph( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei )
      bic = BIC(mod)
      
    } else if (ci_test == "censIndWR") {
      mod = survival::survreg( target ~ ., data = as.data.frame(dataset[, ypografi ]), weights = wei )
      bic = BIC(mod)

    } else if (ci_test == "censIndER") {
      mod = survival::survreg( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei, dist = "exponential" )
      bic = BIC(mod)
      
    } else if (ci_test == "testIndClogit") {
      case = as.logical(target[, 1]);  
      id = target[, 2]
      mod = survival::clogit(case ~ . + strata(id), data = data.frame( dataset[ , ypografi] ) )  ## weights are ignored here anyway
      bic = BIC(mod)
      
    } else if ( ci_test == "testIndLogistic" || ci_test == "gSquare" ) {
      if ( length( Rfast::sort_unique(target) ) == 2 ) {
        #if ( rob == TRUE ) {
        #  mod <- robust::glmRob( target ~ ., data = as.data.frame(dataset[, ypografi ]) , binomial, maxit = 100 )
        #  bic <- mod$deviance + length( coef(mod) ) * log( length(target) )
        #} else { 
          mod = glm( target ~., data = data.frame(dataset[, ypografi ]) , family = binomial, weights = wei ) 
          bic = BIC(mod)
        #}

      } else if ( !is.ordered(target) ) { 
        target = as.factor( as.numeric( as.vector(target) ) )
        mod = nnet::multinom( target ~., data = data.frame(dataset[, ypografi ]) , trace = FALSE, weights = wei )
        bic = BIC(mod)
        
      } else if ( is.ordered(target) ) {
        mod = ordinal::clm( target ~., data = data.frame(dataset[, ypografi ]), weights = wei )
        bic = BIC(mod)
      }
      
      
    }
    
    if ( is.null( colnames(dataset) ) ) {
      names(ypografi) = paste("Var", ypografi, sep = " ")
    } else {
      names(ypografi) = colnames(dataset)[ypografi]
    }
    
    ypografi = c(ypografi, bic)
    names(ypografi)[length(ypografi)] = "bic"
    
  } 

  #############  more than one signatures
 
   if ( nsignat > 1 & nrow(sesObject@signatures) > 1 ) {
      
     if ( nsignat > nrow(sesObject@signatures) )  nsignat = nrow(sesObject@signatures)

    
	 con <- log( length(target) )
     bic <- numeric(nsignat)
     ypografi <- sesObject@signatures[1:nsignat, ] 
     ypografi <- as.matrix(ypografi)
     mod <- list()
    
    for ( i in 1:nsignat ) {
	
     if ( ci_test == "testIndFisher" || ci_test == "testIndReg" ) {
      if ( min(target) > 0 & max(target) < 1 )  target = log( target/(1 - target) ) 

      if ( rob ) {
        # cont = robust::lmRob.control(mxr = 2000, mxf = 2000, mxs = 2000 )  ## only used in robust linear regression
        # mod[[ i ]] = robust::lmRob( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ), control = cont )
        mod[[ i ]] = MASS::rlm(target~., data = data.frame( dataset[, ypografi[i, ] ] ), maxit = 2000, weights = wei )
        bic[i] = cor( target, fitted(mod[[ i ]]) )^2
      } else {
        mod[[ i ]] = lm( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
        bic[i] = BIC( mod[[ i ]] )
      }

     } else if ( ci_test == "testIndSpearman"  ||  ci_test == "testIndRQ" ) {
       if ( all( target>0 & target<1 ) )  target = log( target/(1 - target) ) 
       mod[[ i ]] = quantreg::rq( target ~., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
	     la <- logLik( mod[[ i ]] ) 
       bic[i] <-  - 2 * as.numeric( la ) + attr(la, "df") * con
	  
     } else if ( ci_test == "testIndBeta" ) {
       mod[[ i ]] <- beta.mod( target, dataset[, ypografi[i, ] ], wei = wei )
       bic[i] <-  - 2 * mod[[ i ]]$loglik + ( length( coef(mod[[ i ]]$be) ) + 1 ) * log( length(target) )

     } else if ( ci_test == "testIndSpeedglm" ) {
       la <- length( Rfast::sort_unique(target) )  
       if ( la == 2 ) {
         mod[[ i ]] = speedglm::speedglm( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei, family = binomial() )
         bic[i] =  - 2 * as.numeric( logLik(mod) ) + length( coef(mod) ) * con

       } else if ( la > 2  &  sum( floor(target) - target) == 0 )  {
         mod[[ i ]] = speedglm::speedglm( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei, family = poisson() )
         bic[i] =  - 2 * as.numeric( logLik(mod) ) + length( coef(mod) ) * con
		 
       } else {
         mod[[ i ]] = speedglm::speedlm( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei )
         bic[i] =  - 2 * as.numeric( logLik(mod) ) + ( length( coef(mod) ) + 1 ) * con
       } 
       
     } else if ( ci_test == "testIndPois ") {
       #if ( rob == TRUE ) {
       #  mod[[ i ]] = glm( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ), poisson, maxit = 100 )
       #  bic[i] = mod[[ i ]]$deviance + length( coef(mod[[ i ]]) ) * log( length(target) )
       #} else {  
         mod[[ i ]] = glm( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), family = poisson, weights = wei )
         bic[i] = BIC( mod[[ i ]] )
       #}

     } else if ( ci_test == "testIndNB" ) {
       mod[[ i ]] = MASS::glm.nb( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
       bic[i] = BIC( mod[[ i ]] )

     } else if ( ci_test == "testIndZIP" ) {
       mod[[ i ]] <- zip.mod( target, dataset[, ypografi[i, ] ], wei = wei )
       bic[i] <-  -2 * mod[[ i ]]$loglik + ( length( coef(mod[[ i ]]$be) ) + 1) * log( length(target) )

     } else if ( class(target) == "matrix" || ci_test == "testIndMVreg" ) {
       if ( all(target > 0 & target < 1)  &  sd( Rfast::rowsums(target) ) == 0 )   target = log( target[, -1]/(target[, 1]) ) 
       mod[[ i ]] = lm( target ~.,  data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
       bic[i] = NULL
       
     } else if ( class(target) == "matrix"  &  ci_test == "testIndBinom" ) {
       mod[[ i ]] = glm( target[, 1] /target[, 2] ~., data = data.frame(dataset[, ypografi[i, ] ]), weights = target[, 2], family = binomial )
       bic[i] = BIC(mod[[ i ]])

     } else if (ci_test == "censIndCR") {
       mod[[ i ]] = survival::coxph( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
       bic[i] = BIC( mod[[ i ]] )

     } else if (ci_test == "censIndWR") {
       mod[[ i ]] = survival::survreg( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
       bic[i] = BIC( mod[[ i ]] )

     } else if (ci_test == "censIndER") {
       mod[[ i ]] = survival::survreg( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei, dist = "exponential" )
       bic[i] = BIC( mod[[ i ]] )
       
     } else if (ci_test == "testIndClogit") {
       case = as.logical(target[, 1]);  
       id = target[, 2]
       mod[[ i ]] = survival::clogit(case ~ . + strata(id), data = data.frame( dataset[ , ypografi[i, ] ] ) )  ## weights are ignored here anyway
       bic[i] = BIC( mod[[ i ]] )
       
     } else if ( ci_test == "testIndLogistic" || ci_test == "gSquare" ) {
       if ( length( Rfast::sort_unique(target)) == 2 ) {
        #if ( rob == TRUE ) {
        #  mod[[ i ]] <- robust::glmRob( target ~ ., data = as.data.frame(dataset[, ypografi[i, ] ]) , binomial, maxit = 100 )
        #  bic[[ i ]] <- mod[[ i ]]$deviance + length( coef(mod[[ i ]]) ) * log( length(target) )
        #} else { 
          mod[[ i ]] = glm( target ~., data = data.frame(dataset[, ypografi[i, ] ]) , family = binomial, weights = wei ) 
          bic[[ i ]] = BIC(mod[[ i ]])
        #}

       } else if ( !is.ordered(target) ) { 
         target = as.factor( as.numeric( as.vector(target) ) )
         mod[[ i ]] = nnet::multinom( target ~., data = data.frame( dataset[, ypografi[i, ] ] ), trace = FALSE, weights = wei )
         bic[i] = BIC( mod[[ i ]] )

       } else if ( is.ordered(target) ) {
         mod[[ i ]] = ordinal::clm( target ~., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
         bic[i] = BIC( mod[[ i ]] )
       }
     }
      
    }
  
    ypografi = cbind(ypografi, bic)
	
  }  
  ####### all signatures

  if ( nsignat == "all" ) { 
    ypografi = sesObject@signatures
    bic = numeric( nrow(ypografi) )
    mod = list()
    con <- log( length(target) )
	
    for ( i in 1:nrow(ypografi) ) {
	
     if ( ci_test == "testIndFisher" || ci_test == "testIndReg" ) {
      if ( min(target) > 0 & max(target) < 1 )  target = log( target / (1 - target) ) 
 
      if ( rob ) {
        # cont = robust::lmRob.control(mxr = 2000, mxf = 2000, mxs = 2000 )  ## only used in robust linear regression
        # mod[[ i ]] = robust::lmRob( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ), control = cont )
        mod[[ i ]] = MASS::rlm(target~., data = data.frame( dataset[, ypografi[i, ] ] ), maxit = 2000, weights = wei)
        bic[[ i ]] = cor( target, fitted(mod[[ i ]]) )^2
      } else {
        mod[[ i ]] = lm( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
        bic[i] = BIC( mod[[ i ]] )
      }

     } else if ( ci_test == "testIndSpearman"  ||  ci_test == "testIndRQ" ) {
       if ( all( target > 0 & target < 1 ) )  target = log( target/(1 - target) ) 
       mod[[ i ]] = quantreg::rq( target ~., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
	   la <- logLik( mod[[ i ]] )
       bic[i] <-  - 2 * as.numeric( la ) + attr(la, "df") * con

     } else if ( ci_test == "testIndBeta" ) {
       mod[[ i ]] <- beta.mod( target, dataset[, ypografi[i, ] ], wei = wei )
       bic[i] <-  - 2 * mod[[ i ]]$loglik + ( length( coef(mod[[ i ]]$be) ) + 1 ) * log( length(target) )
       
     } else if ( ci_test == "testIndSpeedglm" ) {
       if ( length( Rfast::sort_unique(target) )  == 2 ) {
         mod[[ i ]] = speedglm::speedglm( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei, family = binomial() )
         bic[i] =  - 2 * as.numeric( logLik(mod) ) + length( coef(mod) ) * con
		 
	      if ( sum( floor(target) - target) == 0 )  {
         mod[[ i ]] = speedglm::speedglm( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei, family = poisson() )
         bic[i] =  - 2 * as.numeric( logLik(mod) ) + length( coef(mod) ) * con
		 
       } else {
         mod[[ i ]] = speedglm::speedlm( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei )
         bic[i] =  - 2 * as.numeric( logLik(mod) ) + ( length( coef(mod) ) + 1 ) * con
       } 

     } else if ( ci_test == "testIndPois ") {
       #if ( rob ) {
       #  mod[[ i ]] = glm( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ), poisson, maxit = 100, weights = wei )
       #  bic[i] = mod[[ i ]]$deviance + length( coef(mod[[ i ]]) ) * log( length(target) )
       #} else {  
         mod[[ i ]] = glm( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), poisson, weights = wei )
         bic[i] = BIC( mod[[ i ]] )
       #}

     } else if ( ci_test == "testIndNB" ) {
       mod[[ i ]] = MASS::glm.nb( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
       bic[i] = BIC( mod[[ i ]] )

     } else if ( ci_test == "testIndZIP" ) {
       mod[[ i ]] <- zip.mod( target, dataset[, ypografi[i, ] ], wei = wei )
       bic[i] <-  -2 * mod[[ i ]]$loglik + ( length( coef(mod[[ i ]]$be) ) + 1) * log( length(target) )

     } else if (class(target) == "matrix" || ci_test == "testIndMVreg" ) {
       if ( all(target > 0 & target < 1)  &  sd( Rfast::rowsums(target) ) == 0 )   target = log( target[, -1]/(target[, 1]) ) 
       mod[[ i ]] = lm( target ~., data = dataset[, ypografi], weights = wei )
       bic[i] = NULL
       
     } else if ( class(target) == "matrix"  &  ci_test == "testIndBinom" ) {      
       mod[[ i ]] = glm( target[, 1] /target[, 2] ~., data = data.frame(dataset[, ypografi[i, ] ]), weights = target[, 2], family = binomial )
       bic[i] = BIC(mod[[ i ]])
       
     } else if (ci_test == "censIndCR") {
       mod[[ i ]] = survival::coxph( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
       bic[i] = BIC( mod[[ i ]] )

     } else if (ci_test == "censIndWR") {
       mod[[ i ]] = survival::survreg( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
       bic[i] = BIC( mod[[ i ]] )
	   
     } else if (ci_test == "censIndER") {
       mod[[ i ]] = survival::survreg( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei, dist = "exponential" )
       bic[i] = BIC( mod[[ i ]] )
       
     } else if (ci_test == "testIndClogit") {
       case = as.logical(target[, 1]);  
       id = target[, 2]
       mod[[ i ]] = survival::clogit(case ~ . + strata(id), data = data.frame( dataset[ , ypografi[i, ] ] ), weights = wei )
       bic[i] = BIC( mod[[ i ]] )

     } else if ( ci_test == "testIndLogistic" || ci_test == "gSquare" ) {
       if ( length( Rfast::sort_unique(target) ) == 2) {
        #if ( rob == TRUE ) {
        #  mod[[ i ]] <- robust::glmRob( target ~ ., data = as.data.frame(dataset[, ypografi[i, ] ]) , binomial, maxit = 100 )
        #  bic[[ i ]] <- mod[[ i ]]$deviance + length( coef(mod[[ i ]]) ) * log( length(target) )
        #} else { 
          mod[[ i ]] = glm( target ~., data = data.frame(dataset[, ypografi[i, ] ]), family = binomial, weights = wei ) 
          bic[[ i ]] = BIC(mod[[ i ]])
        #}

       } else if ( !is.ordered(target) ) { 
         target = as.factor( as.numeric( as.vector(target) ) )
         mod[[ i ]] = nnet::multinom( target ~., data = data.frame( dataset[, ypografi[i, ] ] ), trace = FALSE, weights = wei )
         bic[i] = BIC( mod[[ i ]] )

       } else if ( is.ordered(target) ) {
         mod[[ i ]] = ordinal::clm( target ~., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
         bic[i] = BIC( mod[[ i ]] )
       }
     }
      
    }
    
  }
    
    ypografi = cbind(ypografi, bic)
	
  }  
  
    # if ( is.null( colnames(dataset) ) ) {
    #   names(ypografi) = paste("X", ypografi, sep = "")
    #   colnames(mat1) = paste("+X", ypografi, sep = "")
    #   colnames(mat2) = paste("-X", ypografi, sep = "")
    # } else {
    #   nama = colnames(dataset)
    #   names(ypografi) = nama
    #   colnames(mat1) = paste("+", nama, sep = "")
    #   colnames(mat2) = paste("-", nama, sep = "")
    # } 

  
  list(mod = mod, ypografi = ypografi)  
}
  
 