fs.reg <- function(target, dataset, threshold = 0.05, wei = NULL, test = NULL, user_test = NULL, stopping = "BIC", tol = 2, robust = FALSE, ncores = 1) {
  
  ## target can be Real valued (normal), binary (binomial) or counts (poisson)
  ## dataset is a matrix or a data.frame with the predictor variables
  ## test is the test used, but not really required because depending on the data this will be decided.
  ## there is no hrm is psecifying this though
  ## threshold is the level of significance
  ## method can be either BIC or adjrsq (for non-robust linear models only). The BIC is a consistent method for selecting
  ## models, thus we also use it to avoid overfitting
  ## stopping is based on "BIC"
  ## tol is the tolerance value for the method. If BIC is used as the stopping rule, the default is 2, but usually can be 2 or 4.
  ## If BIC is used as a way to proceed, the tol is 0.
  ## robust is for robust modelling. TRUE or FALSE
  ## ncores is for parallel processing 
  
  threshold <- log(threshold)  ## log of the significance level
  
  p <- ncol(dataset)  ## number of variables
  devi <- dof <- numeric( p )  
  moda <- list()
  k <- 1   ## counter
  n <- length(target)  ## sample size
  tool <- numeric( min(n, p) )
  result <- NULL
  con <- log(n)
  sela <- NULL
  
  #check for NA values in the dataset and replace them with the variable median or the mode
  if(any(is.na(dataset)) == TRUE)
  {

  
    #dataset = as.matrix(dataset);
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    
    if (class(dataset) == "matrix")  {
    
       dataset = apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
              
    }else{
	
    poia <- which( is.na(dataset), arr.ind = TRUE )[2]
 	for( i in poia )
      {
          xi = dataset[, i]
          if(class(xi) == "numeric")
          {                    
            xi[ which( is.na(xi) ) ] = median(xi, na.rm = TRUE) 
          } else if ( class(xi) == "factor" ) {
            xi[ which( is.na(xi) ) ] = levels(xi)[ which.max( as.vector( table(xi) ) )]
          }
          dataset[, i] = xi
        }
    }
    
  }
  
  ##################################
  # target checking and initialize #
  ##################################
  
  
  if ( is.null( colnames(dataset) ) )  {  ## checks for column names
    colnames(dataset) <- paste("X", 1:p, sep = "")
  }	
  
  ## dependent (target) variable checking if no test was given, 
  ## but other arguments are given. For some cases, these are default cases
  
  if ( is.null(test) & is.null(user_test) ) {
    
    ## percentages
    if ( min( target ) > 0 & max( target ) < 1 )  {  ## are they percentages?
      target <- log( target / (1 - target) ) 
      ci_test <- test <- "testIndReg"
    }
    
    ## surival data
    if ( sum( class(target) == "Surv" ) == 1 ) {
      ci_test <- test <- "censIndCR"
    }
    
    ## ordinal, multinomial or perhaps binary data
    if ( is.factor(target) ||  is.ordered(target) || length( unique(target) ) == 2 ) {
      ci_test <- test <- "testIndLogistic"
    }
    
    ## count data
    if ( is.vector(target) ) {
      
      if ( length( unique(target) ) > 2  &  sum( round(target) - target ) == 0 ) {
        ci_test <- test <- "testIndPois"
    
        ## binomial regression 
      
      } else if ( length( unique(target) ) == 2 ) {
        ci_test <- test <- "testIndLogistic"

        ## linear regression 
      } else if ( sum( class(target) == "numeric" ) > 0 ) {
        ci_test <- test <- "testIndReg"  
      }
    }  
    
  }
  
  #available conditional independence tests
  av_models = c("testIndReg", "testIndBinom", "testIndBeta", "censIndCR", "testIndRQ", "censIndWR", "testIndLogistic", "testIndPois", "testIndNB", "testIndZIP", "testIndSpeedglm");
  
  ci_test <- test
  #cat(test)
  
  if ( ( test == "testIndLogistic" & ( length( unique(target) ) == 2 ) )  ||  test == "testIndBinom"  || test == "testIndPois"  ||  ( sum( round(target) - target ) == 0  &  length(target) > 2 )  ) {
    
    result <- glm.fsreg( target, dataset, wei = wei, threshold = exp(threshold), tol = tol, robust = robust, ncores = ncores) 
    
  } else if ( test == "testIndReg"  &  !is.matrix(target) ) {
    
    result <- lm.fsreg( target, dataset, wei = wei, threshold = exp(threshold), stopping = stopping, tol = tol, robust = robust, ncores = ncores ) 
  

  } else if ( test == "testIndSpeedglm" ) {
    
    if ( length( unique(target) ) == 2  ||  ( sum( round(target) - target ) == 0  &  length(target) > 2 ) ) { 
       result <- glm.fsreg( target, dataset, wei = wei, threshold = exp(threshold), tol = tol, heavy = TRUE, robust = robust, ncores = ncores) 
    
    } else  result <- lm.fsreg( target, dataset, wei = wei, threshold = exp(threshold), stopping = stopping, tol = tol, heavy = TRUE, robust = robust, ncores = ncores ) 

    
  } else {
    
    test = match.arg(test, av_models, TRUE);
    #convert to closure type
    
    if ( test == "testIndBeta" ) {
      ci_test <- test 
      test = betareg::betareg 
      robust = FALSE
      stopping = "BIC"
      
    } else if ( test == "censIndCR" ) {
      ci_test <- test
      test = survival::coxph 
      robust = FALSE
      stopping = "BIC"
      
    } else if ( test == "testIndLogistic" ) {
	  
      ci_test <- test
      
	    if ( is.ordered(target) )  {
        test = ordinal::clm
        robust = FALSE
        stopping = "BIC"
      } else {
        test = nnet::multinom
        robust = FALSE
        stopping = "BIC"
      }
	  
    } else if ( test == "testIndNB" ) {
      ci_test <- test
      test = MASS::glm.nb
      robust = FALSE
      stopping = "BIC"

    } else if ( test == "testIndRQ" ) {
      ci_test <- test
      test = quantreg::rq
      robust = FALSE
      stopping = "BIC"
      
    } else if ( test == "testIndZIP" ) {
      ci_test <- test
      test = pscl::zeroinfl 
      robust = FALSE
      stopping = "BIC"
    
    } 
    
	
    if ( !is.null(user_test) )  {
	  test = user_test 
	  ci_test <- "user_test"
	}  
    	  
    
    runtime <- proc.time()
    
    devi = dof = numeric(p)
    ini = test( target ~ 1, weights = wei ) 
    ini =  2 * as.numeric( logLik(ini) )  ## initial 
    
    if (ncores <= 1) {
      for (i in 1:p) {
        mi <- test( target ~ dataset[, i], weights = wei )
        devi[i] <-  2 * as.numeric( logLik(mi) )
        dof[i] = length( coef( mi ) ) 
      }
      
      stat = abs( devi - ini )
      pval = pchisq( stat, dof - 1, lower.tail = FALSE, log.p = TRUE )
      
    } else {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      mata <- matrix(0, p, 2)
      mod <- foreach( i = 1:p, .combine = rbind, .packages = c("MASS", "betareg", "quantreg", "nnet", "survival", "ordinal", "pscl") ) %dopar% {
        ww <- test( target ~ dataset[, i], weights = wei )
        mata[i, ] <- c( 2 * as.numeric( logLik(ww) ), length( coef( ww ) )  )
      }
      
      stopCluster(cl)
      
      stat =  abs( mod[, 1] - ini )
      pval = pchisq( stat, mod[, 2] - 1, lower.tail = FALSE, log.p = TRUE )
    }
    
    mat <- cbind(1:p, pval, stat) 
    
    colnames(mat)[1] <- "variables"
    rownames(mat) <- 1:p
    
    sel <- which.min(pval)
    info <- matrix( numeric(3), ncol = 3 )

    if ( mat[sel, 2] < threshold ) {
      info[1, ] <- mat[sel, ]
      mat <- mat[-sel, ] 
      if ( !is.matrix(mat) )  mat <- matrix(mat, ncol = 3) 
      mat <- mat[ order( mat[, 2] ), ]
      sela <- sel
      mi <- test( target ~ dataset[, sel], weights = wei )
      la <- logLik(mi)
      tool[1] <-  - 2 * as.numeric( la ) +  attr(la, "df") * con
    
      moda[[ 1 ]] <- mi
    }
    
    ############
    ###       k equals 2
    ############ 
    
    if ( info[k, 2] < threshold  &  nrow(mat) > 0 ) {
      
      k <- 2
      pn <- p - k + 1   

      ini =  2 * as.numeric( logLik( moda[[ 1 ]] ) ) 
      do = length( coef( moda[[ 1 ]]  ) ) 
      
      if ( ncores <= 1 ) {
        devi = dof = numeric(pn)
        for ( i in 1:pn ) {
          ww <- test( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), weights = wei )
          devi[i] <-  2 * as.numeric( logLik(ww) )
          dof[i] = length( coef( ww ) )          
        }
        
        stat = abs( devi - ini )
        pval = pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
        
      } else {
        
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mata = matrix(0, pn, 2)  
        mod <- foreach( i = 1:pn, .combine = rbind, .packages = c("MASS", "betareg", "quantreg", "nnet", "survival", "ordinal", "pscl") ) %dopar% {
          ww <- test( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), weights = wei )
          mata[i, ] <- c( 2 * as.numeric( logLik(ww) ), length( coef( ww ) )  )
        }
        
        stopCluster(cl)
        
        stat = abs( mod[, 1] - ini )
        pval = pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
        
      }

    }
    
    mat[, 2:3] <- cbind(pval, stat)
    
    ina <- which.min(mat[, 2])
    sel <- mat[ina, 1]    
    
    if ( mat[ina, 2] < threshold ) {
      ma <- test( target ~ dataset[, sela] + dataset[, sel], weights = wei )
      la <- logLik(ma)
      tool[2]  - 2 * as.numeric( la ) +  attr(la, "df") * con
      
      if ( tool[ 1 ] - tool[ 2 ] <= tol ) {
        info <- info
        
      } else { 
        info <- rbind(info, c( mat[ina, ] ) )
        sela <- info[, 1]
        mat <- mat[-ina , ] 
        if ( !is.matrix(mat) )    mat <- matrix(mat, ncol = 3) 
        mat <- mat[ order( mat[, 2] ), ]
        moda[[ 2 ]] <- ma
      }
      
    } else  info <- info
  
  ############
  ###       k greater than 2
  ############ 
  
  
  if ( nrow(info) > 1  &  nrow(mat) > 0 ) {
    while ( ( info[k, 2] < threshold )  &  ( k < n )  &  ( tool[ k - 1 ] - tool[ k ] > tol )  &  nrow(mat) > 0 )  {
      
      ini =  2 * as.numeric( logLik( moda[[ k ]] ) ) 
      do = length( coef( moda[[ k ]]  ) ) 
      
      k <- k + 1   
      pn <- p - k  + 1
      
      if (ncores <= 1) {  
        devi = dof = numeric(pn) 
        for ( i in 1:pn ) {
          ma <- test( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1] ) ] ), weights = wei )
          devi[i] <-  2 * as.numeric( logLik(ma) ) 
          dof[i] = length( coef( ma ) ) 
        }
        
        stat = abs( devi - ini )
        pval = pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
        
      } else {
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          mata = matrix(0, pn, 2)  
          mod <- foreach( i = 1:pn, .combine = rbind, .packages = c("MASS", "betareg", "quantreg", "nnet", "survival", "ordinal", "pscl") ) %dopar% {
            ww <- test( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), weights = wei )
            mata[i, ] <- c( 2 * as.numeric( logLik(ww) ), length( coef( ww ) )  )
          }
          
          stopCluster(cl)
          
          stat = abs( mod[, 1] - ini )
          pval = pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
          
        }
        
        mat[, 2:3] <- cbind(pval, stat)
        
        ina <- which.min(mat[, 2])
        sel <- mat[ina, 1]    
        
        if ( mat[ina, 2] < threshold ) {
            ma <- test( target ~., data = as.data.frame( dataset[, c(sela, sel) ] ), weights = wei )
            la <- logLik(ma)
            tool[k]  - 2 * as.numeric( la ) +  attr(la, "df") * con
    
 
          if ( tool[ k - 1 ] - tool[ k  ] < tol ) {
            info <- rbind(info, c( 1e300, 0, 0 ) )
            
          } else { 
            info <- rbind( info, mat[ina, ] )
            sela <- info[, 1]
            mat <- mat[-ina , ]
            if ( !is.matrix(mat) ) {
              mat <- matrix(mat, ncol = 3) 
            }
            mat <- mat[ order( mat[, 2] ), ]
            
            moda[[ k ]] <- ma
          } 
          
        } else {
          info <- rbind(info, c( 1e300, 0, 0 ) )
        }
      
    } 
    
  }  
    
    runtime <- proc.time() - runtime
    
    d <- length(sela)
    final <- NULL
    models <- NULL
    
    if ( d >= 1 ) {
      models <- NULL
      xx <- as.data.frame( dataset[, sela] )
      colnames(xx) <- paste("V", sela, sep = "") 
      
      if ( d == 1 ) {
        
        models <- NULL
        xx <- as.data.frame( dataset[, sela] )
        colnames(xx) <- paste("V", sela, sep = "") 
        
        models[[ 1 ]] <- test( target ~., data = data.frame( xx ), weights = wei )

        
      } else {
        for (i in 1:d) {
          models[[ i ]] <- test( target ~., data = data.frame( xx[, 1:i] ), weights = wei )
        }
      }
      
      final <- summary( models[[ d ]] )
      
      info <- info[1:d, ]
      if ( d == 1 )  info <- matrix(info, nrow = 1)
      info <- cbind( info, tool[ 1:d ] ) 
      colnames(info) <- c( "variables", "p-value", "stat", "BIC" )
      rownames(info) <- info[, 1]
    }
    
    result = list(mat = t(mat), info = info, models = models, final = final, ci_test = ci_test, runtime = runtime ) 
    
  }         
  
  result
  
}    










