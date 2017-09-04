# target: the target value
# sesObject: the outcome of the ses
# nisgnat: Number of ypografis and generated models. It could be numeric from 1 to total number of ypografis or "all" for all the 
## ypografis. Default is 1.
waldses.model = function(target, dataset, wei = NULL, wald.sesObject, nsignat = 1, test = NULL) {
  
  if ( sum( is.na(wald.sesObject@signatures) ) > 0 )  mod = paste("No associations were found, hence no model is produced.")
  
  if ( any(is.na(dataset) ) ) {
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    dataset = apply(dataset, 2, function(x) { x[ which(is.na(x)) ] = median(x,na.rm = TRUE) } );
  }
  
  if ( is.null(test) ) {  
    ci_test = wald.sesObject@test
  } else ci_test = test 
  rob = wald.sesObject@rob
  
  if ( nsignat == 1 || ( nsignat > 1 & nrow(wald.sesObject@signatures) == 1 ) ) {
    ypografi = wald.sesObject@selectedVars  
    p <- length(ypografi)
    # mat1 <- mat2 <- numeric(p)
    
    if ( ci_test == "waldMMreg" ) {
      if ( min(target) > 0 & max(target) < 1 )  target = log( target/(1 - target) ) 
        mod = MASS::rlm(target ~., data = data.frame(dataset[, ypografi ]), maxit = 2000, weights = wei, method = "MM" )
        bic = BIC( mod )

    } else if ( ci_test == "waldBeta" ) {     
      mod <- beta.mod( target, dataset[, ypografi ], wei = wei )
      bic <- - 2 * mod$loglik + ( length( coef(mod$be) ) + 1 ) * log( length(target) )
      
    } else if ( ci_test == "waldSpeedglm" ) {
      la <- unique(target) 
      if ( length( la )  == 2 ) {
        mod = speedglm::speedglm( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei, family = binomial() )
        bic =  - 2 * as.numeric( logLik(mod) ) + p * log( length(target) )
        
      } else if ( length( la )  > 2  &  sum( floor(target) - target) == 0 ) {
        mod = speedglm::speedglm( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei, family = poisson() )
        bic =  - 2 * as.numeric( logLik(mod) ) + p * log( length(target) )	  
        
      } else {
        mod = speedglm::speedlm( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei )
        bic =  - 2 * as.numeric( logLik(mod) ) + ( p + 1 ) * log( length(target) )
      }  
      
    } else if ( ci_test == "waldPois") {
      mod <- glm( target ~ ., data = data.frame(dataset[, ypografi ]) , family = poisson, weights = wei )
      bic <- BIC( mod )

    } else if ( ci_test == "waldNB" ) {
      mod = MASS::glm.nb( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei )
      bic = BIC(mod)
      
    } else if ( ci_test == "waldZIP" ) {
      mod <- zip.mod( target, dataset[, ypografi], wei = wei )
      bic <-  -2 * mod$loglik + ( p + 1) * log( length(target) )
      ## bic = BIC(mod)
      
    } else if ( ci_test == "waldBinom" ) {
      mod = glm( target[, 1] /target[, 2] ~., data = data.frame(dataset[, ypografi ]), weights = target[, 2], family = binomial )
      bic = BIC(mod)
      
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
      mod = survival::coxph( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei )
      bic = BIC(mod)
      
    } else if (ci_test == "waldWR") {
      mod = survival::survreg( target ~ ., data = as.data.frame(dataset[, ypografi ]), weights = wei )
      bic = BIC(mod)
      
    } else if (ci_test == "waldER") {
      mod = survival::survreg( target ~ ., data = data.frame(dataset[, ypografi ]), weights = wei, dist = "exponential" )
      bic = BIC(mod)

    } else if ( ci_test == "waldBinary" ) {
      mod = glm( target ~., data = data.frame(dataset[, ypografi ]) , family = binomial, weights = wei ) 
      bic = BIC(mod)
        
    } else if ( ci_test == "waldOrdinal" ) { 
      mod = nnet::multinom( target ~., data = data.frame(dataset[, ypografi ]) , trace = FALSE, weights = wei )
      bic = BIC(mod)

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
  
  if ( nsignat > 1 & nrow(wald.sesObject@signatures) > 1 ) {
    
    if ( nsignat > nrow(wald.sesObject@signatures) )  nsignat = nrow(wald.sesObject@signatures)
    
    con <- log( dim(dataset)[1] )
    bic <- numeric(nsignat)
    ypografi <- wald.sesObject@signatures[1:nsignat, ] 
    ypografi <- as.matrix(ypografi)
    mod <- list()
    
    for ( i in 1:nsignat ) {
      
      if ( ci_test == "waldMMreg" ) {
        if ( min(target) > 0  &  max(target) < 1 )  target = log( target/(1 - target) ) 
          mod[[ i ]] = MASS::rlm(target~., data = data.frame( dataset[, ypografi[i, ] ] ), maxit = 2000, weights = wei )
          bic[i] = BIC( mod[[ i ]] )

      } else if ( ci_test == "waldBeta" ) {
        mod[[ i ]] <- beta.mod( target, dataset[, ypografi[i, ] ], wei = wei )
        bic[i] <-  - 2 * mod[[ i ]]$loglik + ( length( coef(mod[[ i ]]$be) ) + 1 ) * con
        
      } else if ( ci_test == "waldSpeedglm" ) {
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
        
      } else if ( ci_test == "waldPois ") {
        mod[[ i ]] = glm( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), family = poisson, weights = wei )
        bic[i] = BIC( mod[[ i ]] )

      } else if ( ci_test == "waldNB" ) {
        mod[[ i ]] = MASS::glm.nb( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
        bic[i] = BIC( mod[[ i ]] )
        
      } else if ( ci_test == "waldZIP" ) {
        mod[[ i ]] <- zip.mod( target, dataset[, ypografi[i, ] ], wei = wei )
        bic[i] <-  -2 * mod[[ i ]]$loglik + ( length( coef(mod[[ i ]]$be) ) + 1) * con
        
      } else if ( ci_test == "waldBinom" ) {
        mod[[ i ]] = glm( target[, 1] /target[, 2] ~., data = data.frame(dataset[, ypografi[i, ] ]), weights = target[, 2], family = binomial )
        bic[i] = BIC(mod[[ i ]])
        
      } else if ( ci_test == "waldGamma" ) {
        mod[[ i ]] <- glm( target ~ ., data = as.data.frame(dataset[, ypografi[i, ] ]), weights = wei, family = Gamma(link = log) )
        bic[i] <- BIC(mod[[ i ]])
        
      } else if ( ci_test == "waldNormLog" ) {
        mod[[ i ]] <- glm( target ~ ., data = as.data.frame(dataset[, ypografi[i, ] ]), weights = wei, family = gaussian(link = log) )
        bic[i] <- BIC(mod[[ i ]])
        
      } else if ( ci_test == "waldTobit" ) {
        mod[[ i ]] <- survival::survreg( target ~ ., data = as.data.frame(dataset[, ypografi[i, ] ]), weights = wei, dist = "gaussian" )
        bic[i] <-  - 2 * as.numeric( logLik(mod[[ i ]]) ) + ( length( coef(mod[[ i ]]) ) + 1 ) * con
        
      } else if (ci_test == "waldCR") {
        mod[[ i ]] = survival::coxph( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
        bic[i] = BIC( mod[[ i ]] )
        
      } else if (ci_test == "waldWR") {
        mod[[ i ]] = survival::survreg( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
        bic[i] = BIC( mod[[ i ]] )
        
      } else if (ci_test == "waldER") {
        mod[[ i ]] = survival::survreg( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei, dist = "exponential" )
        bic[i] = BIC( mod[[ i ]] )
        
      } else if ( ci_test == "waldBinary" ) {
        mod[[ i ]] = glm( target ~., data = data.frame(dataset[, ypografi[i, ] ]) , family = binomial, weights = wei ) 
        bic[[ i ]] = BIC(mod[[ i ]])

      } else if ( ci_test == "waldOrdinal" ) {
        mod[[ i ]] = ordinal::clm( target ~., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
        bic[i] = BIC( mod[[ i ]] )
      }
      
    }
    
    ypografi = cbind(ypografi, bic)
    
  }  
  ####### all signatures
  
  if ( nsignat == "all" ) { 
    ypografi = wald.sesObject@signatures
    bic = numeric( nrow(ypografi) )
    mod = list()
    con <- log( dim(dataset)[1] )
    
    for ( i in 1:nrow(ypografi) ) {
      
      if ( ci_test == "waldMMreg" ) {
        if ( min(target) > 0 & max(target) < 1 )  target = log( target / (1 - target) ) 
        
          mod[[ i ]] = MASS::rlm(target~., data = data.frame( dataset[, ypografi[i, ] ] ), maxit = 2000, weights = wei)
          bic[[ i ]] = BIC( mod[[ i ]] )

      } else if ( ci_test == "waldBeta" ) {
        mod[[ i ]] <- beta.mod( target, dataset[, ypografi[i, ] ], wei = wei )
        bic[i] <-  - 2 * mod[[ i ]]$loglik + ( length( coef(mod[[ i ]]$be) ) + 1 ) * con
        
      } else if ( ci_test == "waldSpeedglm" ) {
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
          
        } else if ( ci_test == "waldPois ") {
          mod[[ i ]] = glm( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), poisson, weights = wei )
          bic[i] = BIC( mod[[ i ]] )

        } else if ( ci_test == "waldNB" ) {
          mod[[ i ]] = MASS::glm.nb( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
          bic[i] = BIC( mod[[ i ]] )
          
        } else if ( ci_test == "waldZIP" ) {
          mod[[ i ]] <- zip.mod( target, dataset[, ypografi[i, ] ], wei = wei )
          bic[i] <-  -2 * mod[[ i ]]$loglik + ( length( coef(mod[[ i ]]$be) ) + 1) * con

        } else if ( ci_test == "waldBinom" ) {      
          mod[[ i ]] = glm( target[, 1] /target[, 2] ~., data = data.frame(dataset[, ypografi[i, ] ]), weights = target[, 2], family = binomial )
          bic[i] = BIC(mod[[ i ]])
          
        } else if ( ci_test == "waldGamma" ) {
          mod[[ i ]] <- glm( target ~ ., data = as.data.frame(dataset[, ypografi[i, ] ]), weights = wei, family = Gamma(link = log) )
          bic[i] <- BIC(mod[[ i ]])
          
        } else if ( ci_test == "waldNormLog" ) {
          mod[[ i ]] <- glm( target ~ ., data = as.data.frame(dataset[, ypografi[i, ] ]), weights = wei, family = gaussian(link = log) )
          bic[i] <- BIC(mod[[ i ]])
          
        } else if ( ci_test == "waldTobit" ) {
          mod[[ i ]] <- survival::survreg( target ~ ., data = as.data.frame(dataset[, ypografi[i, ] ]), weights = wei, dist = "gaussian" )
          bic[i] <-  - 2 * as.numeric( logLik(mod[[ i ]]) ) + ( length( coef(mod[[ i ]]) ) + 1 ) * con
          
        } else if (ci_test == "waldCR") {
          mod[[ i ]] = survival::coxph( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
          bic[i] = BIC( mod[[ i ]] )
          
        } else if (ci_test == "waldWR") {
          mod[[ i ]] = survival::survreg( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
          bic[i] = BIC( mod[[ i ]] )
          
        } else if (ci_test == "waldER") {
          mod[[ i ]] = survival::survreg( target ~ ., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei, dist = "exponential" )
          bic[i] = BIC( mod[[ i ]] )

        } else if ( ci_test == "waldBinary" ) {
            mod[[ i ]] = glm( target ~., data = data.frame(dataset[, ypografi[i, ] ]), family = binomial, weights = wei ) 
            bic[[ i ]] = BIC(mod[[ i ]])
            
        } else if ( ci_test == "waldOrdinal" ) {
          mod[[ i ]] = ordinal::clm( target ~., data = data.frame( dataset[, ypografi[i, ] ] ), weights = wei )
          bic[i] = BIC( mod[[ i ]] )
          
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

