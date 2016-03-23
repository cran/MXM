# target: the target value
# sesObject: the outcome of the ses
# nisgnat: Number of ypografis and generated models. It could be numeric from 1 to total number of ypografis or "all" for all the 
## ypografis. Default is 1.

model = function(target, dataset, sesObject, nsignat = 1, test = NULL) {
  
  if ( sum( is.na(sesObject@signatures) ) > 0 ) {
    mod = paste("No associations were found, hence no model is produced.")
  }
  
  if ( any(is.na(dataset) ) == TRUE )
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
  
   if ( ci_test == "testIndFisher" || ci_test == "testIndReg" ) {
     if ( min(target) > 0 & max(target) < 1 ) {  ## are they proportions?
       target = log( target/(1 - target) ) 
     }  
     if (rob == TRUE) {
      mod = rlm( target ~ ., data = as.data.frame( dataset[, ypografi ] )  )
      bic = BIC(mod)
     } else {
      mod = lm( target ~ ., data = as.data.frame(dataset[, ypografi ]) )
      bic = BIC(mod)
     }  
   } else if (ci_test == "testIndRQ") {
       if ( all( target>0 & target<1 ) ) {  ## are they proportions?
         target = log( target/(1 - target) ) 
       }
     mod = quantreg::rq( target ~., data = as.data.frame(dataset[, ypografi ])  )
     bic = BIC(mod)
    } else if ( ci_test == "testIndBeta" ) {
      mod = betareg::betareg( target ~ ., data = as.data.frame(dataset[, ypografi ])  )
      bic = BIC(mod)
    } else if ( ci_test == "testIndPois ") {
      mod = glm( target ~ ., data = as.data.frame(dataset[, ypografi ]) , poisson )
      bic = BIC(mod)
    } else if ( ci_test == "testIndNB" ) {
      mod = MASS::glm.nb( target ~ ., data = as.data.frame(dataset[, ypografi ])  )
      bic = BIC(mod)
    } else if ( ci_test == "testIndZIP" ) {
      mod = pscl::zeroinfl( target ~. | ., data = as.data.frame( dataset[, ypografi] ) )
      bic = BIC(mod)
    } else if ( class(target) == "matrix" || ci_test == "testIndMVreg" ) {
       if ( all(target > 0 & target < 1) ) {  ## are they compositional data?
         target = log( target[, -1]/(target[, 1]) ) 
       } 
      mod = lm( target ~., data = as.data.frame(dataset[, ypografi ]) )
      bic = BIC(mod)
    } else if (ci_test == "censIndCR") {
      mod = survival::coxph( target ~ ., data = as.data.frame(dataset[, ypografi ]) )
      bic = BIC(mod)
    } else if (ci_test == "censIndWR") {
      mod = survival::survreg( target ~ ., data = as.data.frame(dataset[, ypografi ]) )
      bic = BIC(mod)
    } else if ( ci_test == "testIndLogistic" || ci_test == "gSquare" ) {
      if ( length(unique(target)) == 2 ) {
        mod = glm( target ~., data = as.data.frame(dataset[, ypografi ]) , binomial ) 
        bic = BIC(mod)
      } else if ( is.ordered(target) == FALSE ) { 
        target = as.factor( as.numeric( as.vector(target) ) )
        mod = nnet::multinom( target ~., data = as.data.frame(dataset[, ypografi ]) , trace = FALSE )
        bic = BIC(mod)
      } else if ( is.ordered(target) == TRUE ) {
        mod = ordinal::clm( target ~., data = as.data.frame(dataset[, ypografi ])  )
        bic = BIC(mod)
      }
    }
    ypografi = c(ypografi, bic)
    names(ypografi)[length(ypografi)] = "bic"
  } 

    if ( nsignat > 1 & nrow(sesObject@signatures) > 1 ) {
      
      if ( nsignat > nrow(sesObject@signatures) ) {
        nsignat = nrow(sesObject@signatures)
      }
    
    bic = numeric(nsignat)
    ypografi = sesObject@signatures[1:nsignat, ] 
    ypografi = as.matrix(ypografi)
    mod = list()
    
    for ( i in 1:nsignat ) {
     if ( ci_test == "testIndFisher" || ci_test == "testIndReg" ) {
      if ( min(target) > 0 & max(target) < 1 ) {  ## are they proportions?
       target = log( target/(1 - target) ) 
      }  
      if (rob == TRUE) {
        mod[[ i ]] = rlm( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
        bic[i] = BIC( mod[[ i ]] )
      } else {
        mod[[ i ]] = lm( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
        bic[i] = BIC( mod[[ i ]] )
      }
     } else if (ci_test == "testIndRQ") {
        if ( all( target>0 & target<1 ) ) {  ## are they proportions?
         target = log( target/(1 - target) ) 
        }
       mod[[ i ]] = quantreg::rq( target ~., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
       bic[i] = BIC( mod[[ i ]] )
     } else if ( ci_test == "testIndBeta" ) {
       mod[[ i ]] = betareg::betareg( target ~ ., as.data.frame( dataset[, ypografi[i, ] ] ) )
       bic[i] = BIC( mod[[ i ]] )
     } else if ( ci_test == "testIndPois ") {
       mod[[ i ]] = glm( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ), poisson )
       bic[i] = BIC( mod[[ i ]] )
     } else if ( ci_test == "testIndNB" ) {
       mod[[ i ]] = MASS::glm.nb( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
       bic[i] = BIC( mod[[ i ]] )
     } else if ( ci_test == "testIndZIP" ) {
       mod[[ i ]] = pscl::zeroinfl( target ~. | ., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
       bic = BIC(mod)
     } else if ( class(target) == "matrix" || ci_test == "testIndMVreg" ) {
       if ( all(target > 0 & target < 1) ) {  ## are they compositional data?
         target = log( target[, -1] / (target[, 1]) ) 
       } 
       mod = lm( target ~.,  data = as.data.frame( dataset[, ypografi[i, ] ] ) )
       bic[i] = BIC( mod[[ i ]] )
     } else if (ci_test == "censIndCR") {
       mod[[ i ]] = survival::coxph( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
       bic[i] = BIC( mod[[ i ]] )
     } else if (ci_test == "censIndWR") {
       mod[[ i ]] = survival::survreg( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
       bic[i] = BIC( mod[[ i ]] )
     } else if ( ci_test == "testIndLogistic" || ci_test == "gSquare" ) {
       if ( length(unique(target)) == 2 ) {
         mod[[ i ]] = glm( target ~., data = as.data.frame( dataset[, ypografi[i, ] ] ), binomial ) 
         bic[i] = BIC( mod[[ i ]] )
       } else if ( is.ordered(target) == FALSE ) { 
         target = as.factor( as.numeric( as.vector(target) ) )
         mod[[ i ]] = nnet::multinom( target ~., data = as.data.frame( dataset[, ypografi[i, ] ] ), trace = FALSE )
         bic[i] = BIC( mod[[ i ]] )
       } else if ( is.ordered(target) == TRUE ) {
         mod[[ i ]] = ordinal::clm( target ~., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
         bic[i] = BIC( mod[[ i ]] )
       }
     }
    }
    ypografi = cbind(ypografi, bic)
  }
   
  if ( nsignat == "all" ) { 
    ypografi = sesObject@signatures
    bic = numeric( nrow(ypografi) )
    mod = list()
    for ( i in 1:nrow(ypografi) ) {
     if ( ci_test == "testIndFisher" || ci_test == "testIndReg" ) {
      if ( min(target) > 0 & max(target) < 1 ) {  ## are they proportions?
        target = log( target / (1 - target) ) 
      }  
      if (rob == TRUE) {
        mod[[ i ]] = rlm( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
        bic[i] = BIC( mod[[ i ]] )
      } else {
       mod[[ i ]] = lm( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
       bic[i] = BIC( mod[[ i ]] )
      }
     } else if (ci_test == "testIndRQ") {
        if ( all( target > 0 & target < 1 ) ) {  ## are they proportions?
           target = log( target/(1 - target) ) 
        }
       mod[[ i ]] = quantreg::rq( target ~., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
       bic[i] = BIC( mod[[ i ]] )
     } else if ( ci_test == "testIndBeta" ) {
       mod[[ i ]] = betareg::betareg( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
       bic[i] = BIC( mod[[ i ]] )
     } else if ( ci_test == "testIndPois ") {
       mod[[ i ]] = glm( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ), poisson )
       bic[i] = BIC( mod[[ i ]] )
     } else if ( ci_test == "testIndNB" ) {
       mod[[ i ]] = MASS::glm.nb( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
       bic[i] = BIC( mod[[ i ]] )
     } else if ( ci_test == "testIndZIP" ) {
       mod[[ i ]] = pscl::zeroinfl( target ~. | ., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
       bic = BIC(mod)
     } else if (class(target) == "matrix" || ci_test == "testIndMVreg" ) {
       if ( all(target > 0 & target < 1) ) {  ## are they compositional data?
         target = log( target[, -1] / (target[, 1]) ) 
       } 
       mod = lm( target ~., data = dataset[, ypografi] )
       bic[i] = BIC( mod[[ i ]] )
     } else if (ci_test == "censIndCR") {
       mod[[ i ]] = survival::coxph( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
       bic[i] = BIC( mod[[ i ]] )
     } else if (ci_test == "censIndWR") {
       mod[[ i ]] = survival::survreg( target ~ ., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
       bic[i] = BIC( mod[[ i ]] )
     } else if ( ci_test == "testIndLogistic" || ci_test == "gSquare" ) {
       if ( length(unique(target)) == 2) {
         mod[[ i ]] = glm( target ~., data = as.data.frame( dataset[, ypografi[i, ] ] ), binomial ) 
         bic[i] = BIC( mod[[ i ]] )
       } else if ( is.ordered(target) == FALSE ) { 
         target = as.factor( as.numeric( as.vector(target) ) )
         mod[[ i ]] = nnet::multinom( target ~., data = as.data.frame( dataset[, ypografi[i, ] ] ), trace = FALSE )
         bic[i] = BIC( mod[[ i ]] )
       } else if ( is.ordered(target) == TRUE ) {
         mod[[ i ]] = ordinal::clm( target ~., data = as.data.frame( dataset[, ypografi[i, ] ] ) )
         bic[i] = BIC( mod[[ i ]] )
       }
     }
    }
    ypografi = cbind(ypografi, bic)
  } 
  list(mod = mod, ypografi = ypografi)  
}
  
 
      