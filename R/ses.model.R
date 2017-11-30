# target: the target value
# sesObject: the outcome of the ses
# nisgnat: Number of ypografis and generated models. It could be numeric from 1 to total number of ypografis or "all" for all the 
## ypografis. Default is 1.
ses.model = function(target, dataset, wei = NULL, sesObject, nsignat = 1, test = NULL) {
  
  if ( sum( is.na(sesObject@signatures) ) > 0 )  mod = paste("No associations were found, hence no model is produced.")
  
  if ( any(is.na(dataset) ) ) {
    #dataset = as.matrix(dataset);
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    if (class(dataset) == "matrix")  {
      dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
    } else {
      poia <- which( is.na(dataset), arr.ind = TRUE )[2]
      for( i in poia )  {
        xi <- dataset[, i]
        if(class(xi) == "numeric") {                    
          xi[ which( is.na(xi) ) ] <- median(xi, na.rm = TRUE) 
        } else if ( is.factor( xi ) ) {
          xi[ which( is.na(xi) ) ] <- levels(xi)[ which.max( as.vector( table(xi) ) )]
        }
        dataset[, i] <- xi
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

     if ( rob ) {
       mod = MASS::rlm(target ~., data = data.frame(dataset[, ypografi ]), maxit = 2000, weights = wei )
       bic = BIC( mod )
      
     } else {
       mod = lm( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei )
       bic = BIC(mod)     
     }  
     
   } else if ( ci_test == "testIndSpearman"  ||  ci_test == "testIndRQ" ) {
      mod <- quantreg::rq( target ~., data = data.frame(dataset[, ypografi ]), weights = wei )
	    la <- logLik(mod)
      bic <-  - 2 * as.numeric( la ) + attr(la, "df") * log( length(target) )
      
    } else if ( ci_test == "testIndBeta" ) {     
      mod <- beta.mod( target, dataset[, ypografi ], wei = wei )
      bic <- - 2 * mod$loglik + ( length( mod$be ) + 1 ) * log( length(target) )
      
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
        mod <- glm( target ~ ., data = data.frame(dataset[, ypografi ]) , family = poisson, weights = wei )
        bic <- BIC( mod )

    } else if ( ci_test == "testIndNB" ) {
      mod = MASS::glm.nb( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei )
      bic = BIC(mod)
      
    } else if ( ci_test == "testIndZIP" ) {
      mod <- zip.mod( target, dataset[, ypografi], wei = wei )
      bic <-  -2 * mod$loglik + ( length( coef(mod$be) ) + 1) * log( length(target) )
      ## bic = BIC(mod)
      
    } else if ( is.matrix(target) &  ci_test == "testIndMVreg" ) {
      if ( all(target > 0 & target < 1)  &  Rfast::Var( Rfast::rowsums(target) ) == 0 )   target = log( target[, -1]/(target[, 1]) ) 
       mod = lm( target ~., data = data.frame(dataset[, ypografi ]), weights = wei )
       bic = NULL
      
    } else if ( is.matrix(target) &  ci_test == "testIndBinom" ) {
      mod = glm( target[, 1] /target[, 2] ~., data = data.frame(dataset[, ypografi ]), weights = target[, 2], family = binomial )
      bic = BIC(mod)
      
    } else if ( ci_test == "testIndGamma" ) {
      mod <- glm( target ~ ., data = as.data.frame(dataset[, ypografi ]), weights = wei, family = Gamma(link = log) )
      bic <- BIC(mod)
      
    } else if ( ci_test == "testIndNormLog" ) {
      mod <- glm( target ~ ., data = as.data.frame(dataset[, ypografi ]), weights = wei, family = gaussian(link = log) )
      bic <- BIC(mod)
      
    } else if ( ci_test == "testIndTobit" ) {
      mod <- survival::survreg( target ~ ., data = as.data.frame(dataset[, ypografi ]), weights = wei, dist = "gaussian" )
      bic <-  - 2 * logLik(mod) + (length(mod$coefficients) + 1) * log( NROW(dataset)  )
      
    } else if (ci_test == "censIndCR") {
      mod = survival::coxph( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei )
      bic = BIC(mod)
      
    } else if (ci_test == "censIndWR") {
      mod = survival::survreg( target ~ ., data = as.data.frame(dataset[, ypografi ]), weights = wei )
      bic =  - 2 * logLik(mod) + (length(mod$coefficients) + 1) * log( NROW(dataset)  )

    } else if (ci_test == "censIndER") {
      mod = survival::survreg( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei, dist = "exponential" )
      bic = - 2 * logLik(mod) + (length(mod$coefficients) + 1) * log( NROW(dataset)  )
      
    } else if (ci_test == "testIndClogit") {
      case = as.logical(target[, 1]);  
      id = target[, 2]
      mod = survival::clogit(case ~ . + strata(id), data = data.frame( dataset[ , ypografi] ) )  ## weights are ignored here anyway
      bic = BIC(mod)
      
    } else if ( ci_test == "testIndLogistic" || ci_test == "gSquare" ) {
      if ( length( unique(target) ) == 2 ) {
        mod = glm( target ~., data = data.frame(dataset[, ypografi ]) , family = binomial, weights = wei ) 
        bic = BIC(mod)

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
    } else    names(ypografi) = colnames(dataset)[ypografi]
    
    ypografi = c(ypografi, bic)
    names(ypografi)[length(ypografi)] = "bic"
  } 

  #############  more than one signatures
 
   if ( nsignat > 1 & nrow(sesObject@signatures) > 1 ) {
      
     if ( nsignat > nrow(sesObject@signatures) )  nsignat = nrow(sesObject@signatures)

	   con <- log( NROW(dataset) )
     bic <- numeric(nsignat)
     ypografi <- sesObject@signatures[1:nsignat, ] 
     ypografi <- as.matrix(ypografi)
     mod <- list()
    
    for ( i in 1:nsignat ) {
	
     if ( ci_test == "testIndFisher" || ci_test == "testIndReg" ) {
      if ( rob ) {
        mod[[ i ]] = MASS::rlm(target~., data = data.frame( dataset[, ypografi[i, ] ] ), maxit = 2000, weights = wei )
        bic[i] = BIC( mod[[ i ]] ) 
      } else {
        mod[[ i ]] = lm( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
        bic[i] = BIC( mod[[ i ]] )
      }

     } else if ( ci_test == "testIndSpearman"  ||  ci_test == "testIndRQ" ) {
       mod[[ i ]] = quantreg::rq( target ~., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
	     la <- logLik( mod[[ i ]] ) 
       bic[i] <-  - 2 * as.numeric( la ) + attr(la, "df") * con
	  
     } else if ( ci_test == "testIndBeta" ) {
       mod[[ i ]] <- beta.mod( target, dataset[, ypografi[i, ] ], wei = wei )
       bic[i] <-  - 2 * mod[[ i ]]$loglik + ( length( mod[[ i ]]$be ) + 1 ) * con

     } else if ( ci_test == "testIndSpeedglm" ) {
       la <- length( unique(target) )  
       if ( la == 2 ) {
         mod[[ i ]] = speedglm::speedglm( target ~ ., data = data.frame(dataset[, ypografi[i, ] ]), weights = wei, family = binomial() )
         bic[i] =  - 2 * as.numeric( logLik(mod[[ i ]]) ) + length( coef(mod[[ i ]]) ) * con

       } else if ( la > 2  &  sum( floor(target) - target) == 0 )  {
         mod[[ i ]] = speedglm::speedglm( target ~ ., data = data.frame(dataset[, ypografi[i, ] ]), weights = wei, family = poisson() )
         bic[i] =  - 2 * as.numeric( logLik(mod[[ i ]]) ) + length( coef(mod[[ i ]]) ) * con
		 
       } else {
         mod[[ i ]] = speedglm::speedlm( target ~ ., data = data.frame(dataset[, ypografi[i, ] ]), weights = wei )
         bic[i] =  - 2 * as.numeric( logLik(mod[[ i ]]) ) + ( length( coef(mod[[ i ]]) ) + 1 ) * con
       } 
       
     } else if ( ci_test == "testIndPois ") {
       mod[[ i ]] = glm( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), family = poisson, weights = wei )
       bic[i] = BIC( mod[[ i ]] )

     } else if ( ci_test == "testIndNB" ) {
       mod[[ i ]] = MASS::glm.nb( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
       bic[i] = BIC( mod[[ i ]] )

     } else if ( ci_test == "testIndZIP" ) {
       mod[[ i ]] <- zip.mod( target, dataset[, ypografi[i, ] ], wei = wei )
       bic[i] <-  -2 * mod[[ i ]]$loglik + ( length( coef(mod[[ i ]]$be) ) + 1) * con

     } else if ( is.matrix(target)  || ci_test == "testIndMVreg" ) {
       if ( all(target > 0 & target < 1)  &  Rfast::Var( Rfast::rowsums(target) ) == 0 )   target = log( target[, -1]/(target[, 1]) ) 
       mod[[ i ]] = lm( target ~.,  data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
       bic[i] = NULL
       
     } else if ( is.matrix(target)  &  ci_test == "testIndBinom" ) {
       mod[[ i ]] = glm( target[, 1] /target[, 2] ~., data = data.frame(dataset[, ypografi[i, ] ]), weights = target[, 2], family = binomial )
       bic[i] = BIC(mod[[ i ]])
       
     } else if ( ci_test == "testIndGamma" ) {
       mod[[ i ]] <- glm( target ~ ., data = as.data.frame(dataset[, ypografi[i, ] ]), weights = wei, family = Gamma(link = log) )
       bic[i] <- BIC(mod[[ i ]])
       
     } else if ( ci_test == "testIndNormLog" ) {
       mod[[ i ]] <- glm( target ~ ., data = as.data.frame(dataset[, ypografi[i, ] ]), weights = wei, family = gaussian(link = log) )
       bic[i] <- BIC(mod[[ i ]])
       
     } else if ( ci_test == "testIndTobit" ) {
       mod[[ i ]] <- survival::survreg( target ~ ., data = as.data.frame(dataset[, ypografi[i, ] ]), weights = wei, dist = "gaussian" )
       bic[i] <-  - 2 * as.numeric( logLik(mod[[ i ]]) ) + ( length( coef(mod[[ i ]]) ) + 1 ) * con

     } else if (ci_test == "censIndCR") {
       mod[[ i ]] = survival::coxph( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
       bic[i] = BIC( mod[[ i ]] )

     } else if (ci_test == "censIndWR") {
       mod[[ i ]] = survival::survreg( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
       bic[i] = - 2 * logLik(mod[[ i ]]) + (length(mod[[ i ]]$coefficients) + 1) * con

     } else if (ci_test == "censIndER") {
       mod[[ i ]] = survival::survreg( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei, dist = "exponential" )
       bic[i] = - 2 * logLik(mod[[ i ]]) + (length(mod[[ i ]]$coefficients) + 1) * con
       
     } else if (ci_test == "testIndClogit") {
       case = as.logical(target[, 1]);  
       id = target[, 2]
       mod[[ i ]] = survival::clogit(case ~ . + strata(id), data = data.frame( dataset[ , ypografi[i, ] ] ) )  ## weights are ignored here anyway
       bic[i] = BIC( mod[[ i ]] )
       
     } else if ( ci_test == "testIndLogistic" || ci_test == "gSquare" ) {
       if ( length( unique(target)) == 2 ) {
          mod[[ i ]] = glm( target ~., data = data.frame(dataset[, ypografi[i, ] ]) , family = binomial, weights = wei ) 
          bic[[ i ]] = BIC(mod[[ i ]])

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
    con <- log( NROW(dataset) )
	
    for ( i in 1:nrow(ypografi) ) {
	
     if ( ci_test == "testIndFisher" || ci_test == "testIndReg" ) {
      if ( rob ) {
        mod[[ i ]] = MASS::rlm(target~., data = data.frame( dataset[, ypografi[i, ] ] ), maxit = 2000, weights = wei)
        bic[[ i ]] = BIC( mod[[ i ]] )
      } else {
        mod[[ i ]] = lm( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
        bic[i] = BIC( mod[[ i ]] )
      }

     } else if ( ci_test == "testIndSpearman"  ||  ci_test == "testIndRQ" ) {
       mod[[ i ]] = quantreg::rq( target ~., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
	     la <- logLik( mod[[ i ]] )
       bic[i] <-  - 2 * as.numeric( la ) + attr(la, "df") * con

     } else if ( ci_test == "testIndBeta" ) {
       mod[[ i ]] <- beta.mod( target, dataset[, ypografi[i, ] ], wei = wei )
       bic[i] <-  - 2 * mod[[ i ]]$loglik + ( length( coef(mod[[ i ]]$be) ) + 1 ) * con
       
     } else if ( ci_test == "testIndSpeedglm" ) {
       if ( length( unique(target) )  == 2 ) {
         mod[[ i ]] = speedglm::speedglm( target ~ ., data = data.frame(dataset[, ypografi[i, ] ]), weights = wei, family = binomial() )
         bic[i] =  - 2 * as.numeric( logLik(mod[[ i ]]) ) + length( coef(mod[[ i ]]) ) * con
		 
	      if ( sum( floor(target) - target) == 0 )  {
         mod[[ i ]] = speedglm::speedglm( target ~ ., data = data.frame(dataset[, ypografi[i, ] ]), weights = wei, family = poisson() )
         bic[i] =  - 2 * as.numeric( logLik(mod[[ i ]]) ) + length( coef(mod[[ i ]]) ) * con
		 
       } else {
         mod[[ i ]] = speedglm::speedlm( target ~ ., data = data.frame(dataset[, ypografi[i, ] ]), weights = wei )
         bic[i] =  - 2 * as.numeric( logLik(mod[[ i ]]) ) + ( length( coef(mod[[ i ]]) ) + 1 ) * con
       } 

     } else if ( ci_test == "testIndPois ") {
       mod[[ i ]] = glm( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), poisson, weights = wei )
       bic[i] = BIC( mod[[ i ]] )

     } else if ( ci_test == "testIndNB" ) {
       mod[[ i ]] = MASS::glm.nb( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
       bic[i] = BIC( mod[[ i ]] )

     } else if ( ci_test == "testIndZIP" ) {
       mod[[ i ]] <- zip.mod( target, dataset[, ypografi[i, ] ], wei = wei )
       bic[i] <-  -2 * mod[[ i ]]$loglik + ( length( coef(mod[[ i ]]$be) ) + 1) * log( length(target) )

     } else if ( is.matrix(target)  ||  ci_test == "testIndMVreg" ) {
       if ( all(target > 0 & target < 1)  &  Rfast::Var( Rfast::rowsums(target) ) == 0 )   target = log( target[, -1]/(target[, 1]) ) 
       mod[[ i ]] = lm( target ~., data = dataset[, ypografi[i, ]], weights = wei )
       bic[i] = NULL
       
     } else if ( is.matrix(target)  &  ci_test == "testIndBinom" ) {      
       mod[[ i ]] = glm( target[, 1] /target[, 2] ~., data = data.frame(dataset[, ypografi[i, ] ]), weights = target[, 2], family = binomial )
       bic[i] = BIC(mod[[ i ]])
       
     } else if ( ci_test == "testIndGamma" ) {
       mod[[ i ]] <- glm( target ~ ., data = as.data.frame(dataset[, ypografi[i, ] ]), weights = wei, family = Gamma(link = log) )
       bic[i] <- BIC(mod[[ i ]])
       
     } else if ( ci_test == "testIndNormLog" ) {
       mod[[ i ]] <- glm( target ~ ., data = as.data.frame(dataset[, ypografi[i, ] ]), weights = wei, family = gaussian(link = log) )
       bic[i] <- BIC(mod[[ i ]])
       
     } else if ( ci_test == "testIndTobit" ) {
       mod[[ i ]] <- survival::survreg( target ~ ., data = as.data.frame(dataset[, ypografi[i, ] ]), weights = wei, dist = "gaussian" )
       bic[i] <-  - 2 * as.numeric( logLik(mod[[ i ]]) ) + ( length( coef(mod[[ i ]]) ) + 1 ) * con
       
     } else if (ci_test == "censIndCR") {
       mod[[ i ]] = survival::coxph( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
       bic[i] = BIC( mod[[ i ]] )

     } else if (ci_test == "censIndWR") {
       mod[[ i ]] = survival::survreg( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
       bic[i] = - 2 * logLik(mod[[ i ]]) + (length(mod[[ i ]]$coefficients) + 1) * con
	   
     } else if (ci_test == "censIndER") {
       mod[[ i ]] = survival::survreg( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei, dist = "exponential" )
       bic[i] = BIC( mod[[ i ]] )
       
     } else if (ci_test == "testIndClogit") {
       case = as.logical(target[, 1]);  
       id = target[, 2]
       mod[[ i ]] = survival::clogit(case ~ . + strata(id), data = data.frame( dataset[ , ypografi[i, ] ] ), weights = wei )
       bic[i] = BIC( mod[[ i ]] )

     } else if ( ci_test == "testIndLogistic" || ci_test == "gSquare" ) {
       if ( length( unique(target) ) == 2) {
          mod[[ i ]] = glm( target ~., data = data.frame(dataset[, ypografi[i, ] ]), family = binomial, weights = wei ) 
          bic[[ i ]] = BIC(mod[[ i ]])
  
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
  
 