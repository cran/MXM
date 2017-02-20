beta.fsreg <- function(target, dataset, threshold = 0.05, wei = NULL, tol = 2, ncores = 1) {
  
  ## target is a proportion in the (0, 1) interval
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
  if( any(is.na(dataset)) ) {
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    if (class(dataset) == "matrix")  {
      dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
    } else {
      poia <- which( is.na(dataset), arr.ind = TRUE )[2]
      for( i in poia )  {
        xi <- dataset[, i]
        if ( is.numeric(xi) ) {                    
          xi[ which( is.na(xi) ) ] <- median(xi, na.rm = TRUE) 
        } else if ( is.factor( xi ) )     xi[ which( is.na(xi) ) ] <- levels(xi)[ which.max( as.vector( table(xi) ) )]
        dataset[, i] <- xi
      }
    }
  }
  
  if ( is.null( colnames(dataset) ) )  colnames(dataset) <- paste("X", 1:p, sep = "") 

    runtime <- proc.time()
    
    if ( is.null(wei) ) {
      ini <-  - 2 * Rfast::beta.mle(target)$loglik + 2 * con 
    } else ini <-  - 2 * betamle.wei(target, wei)$loglik + 2 * con 
       
    mod <- beta.regs(target, dataset, wei, logged = TRUE, ncores = ncores)[, 1:2]
    mat <- cbind(1:p, mod[, 2], mod[, 1])
    rownames(mat) <- 1:p
    colnames(mat) <- c( "variables", "log.p-value", "stat" )
    
    sel <- which.min(mat[, 2])
    info <- matrix( numeric(3), ncol = 3 )
    
    if ( mat[sel, 2] < threshold ) {
      info[1, ] <- mat[sel, ]
      mat <- mat[-sel, ] 
      if ( !is.matrix(mat) )  mat <- matrix(mat, ncol = 3) 
      sela <- sel
      mi <- beta.reg( target, dataset[, sel], wei = wei )
      tool[1] <-  - 2 * mi$loglik + ( length(mi$be) + 1 ) * con
      moda[[ 1 ]] <- mi
    }  else  {
      info <- info  
      sela <- NULL
    }
    
    ############
    ###       k equals 2
    ############ 
    
    if ( info[k, 2] < threshold  &  nrow(mat) > 0 ) {
      
      k <- 2
      pn <- p - k + 1   
      
      ini <- 2 * mi$loglik
      do <- length(mi$be)

      if ( ncores <= 1 ) {
        devi <- dof <- numeric(pn)
        for ( i in 1:pn ) {
          ww <- beta.reg( target, dataset[, c(sela, mat[i, 1]) ], wei = wei )
          devi[i] <-  2 * ww$loglik
          dof[i] <- length( ww$be )          
        }
        
        stat <- abs( devi - ini )
        pval <- pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
        
      } else {
        
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mata <- matrix(0, pn, 2)  
        mod <- foreach( i = 1:pn, .combine = rbind, .export = "beta.reg" ) %dopar% {
          ww <- beta.reg( target, dataset[, c(sela, mat[i, 1]) ], wei = wei )
          mata[i, ] <- c( 2 * ww$loglik, length( ww$be )  )
        }
        
        stopCluster(cl)
        
        stat <- abs( mod[, 1] - ini )
        pval <- pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
        
      }
    
    mat[, 2:3] <- cbind(pval, stat)
    ina <- which.min(mat[, 2])
    sel <- mat[ina, 1]    
    
    if ( mat[ina, 2] < threshold ) {
      ma <- beta.reg( target, dataset[, c(sela, sel)], wei = wei )
      tool[2] <-  - 2 * ma$loglik + ( length(ma$be) + 1 ) * con
      
      if ( tool[ 1 ] - tool[ 2 ] <= tol ) {
        info <- info
        
      } else { 
        info <- rbind(info, c( mat[ina, ] ) )
        sela <- info[, 1]
        mat <- mat[-ina , ] 
        if ( !is.matrix(mat) )    mat <- matrix(mat, ncol = 3) 
        moda[[ 2 ]] <- ma
      }
      
    } else  info <- info
    
  }  
    ############
    ###       k greater than 2
    ############ 
    
    if ( nrow(info) > 1  &  nrow(mat) > 0 ) {
      while ( ( info[k, 2] < threshold )  &  ( k < n )  &  ( tool[ k - 1 ] - tool[ k ] > tol )  &  nrow(mat) > 0 )  {
        
        ini <-  2 * moda[[ k ]]$loglik 
        do <- length( moda[[ k ]]$be ) 
        
        k <- k + 1   
        pn <- p - k  + 1
        
        if (ncores <= 1) {  
          devi <- dof <- numeric(pn) 
          for ( i in 1:pn ) {
            ma <- beta.reg( target, dataset[, c(sela, mat[i, 1] ) ], wei = wei )
            devi[i] <-  2 * ma$loglik
            dof[i] <- length( ma$be ) 
          }
          
          stat <- abs( devi - ini )
          pval <- pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
          
        } else {
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          mata <- matrix(0, pn, 2)  
          mod <- foreach( i = 1:pn, .combine = rbind, .export = "beta.reg" ) %dopar% {
            ww <- beta.reg( target, dataset[, c(sela, mat[i, 1]) ], wei = wei )
            mata[i, ] <- c( 2 * ww$loglik, length( ww$be )  )
          }
          stopCluster(cl)
          
          stat <- abs( mod[, 1] - ini )
          pval <- pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
          
        }
        
        mat[, 2:3] <- cbind(pval, stat)
        
        ina <- which.min(mat[, 2])
        sel <- mat[ina, 1]    
        
        if ( mat[ina, 2] < threshold ) {
          ma <- beta.reg( target, dataset[, c(sela, sel) ], wei = wei )
          tool[k] <-  - 2 * ma$loglik + ( length(ma$be) + 1 ) * con
          
          
          if ( tool[ k - 1 ] - tool[ k  ] <= tol ) {
            info <- rbind(info, c( 1e300, 0, 0 ) )
            
          } else { 
            info <- rbind( info, mat[ina, ] )
            sela <- info[, 1]
            mat <- mat[-ina , ]
            if ( !is.matrix(mat) )  mat <- matrix(mat, ncol = 3) 
            moda[[ k ]] <- ma
          } 
          
        } else    info <- rbind(info, c( 1e300, 0, 0 ) )
        
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
        models[[ 1 ]] <- beta.reg( target, xx, wei = wei )
        
      } else    for (i in 1:d)    models[[ i ]] <- beta.reg( target, xx[, 1:i], wei = wei )
      
      final <- models[[ d ]] 
      info <- info[1:d, ]
      if ( d == 1 )  info <- matrix(info, nrow = 1)
      info <- cbind( info, tool[ 1:d ] ) 
      colnames(info) <- c( "variables", "log.p-value", "stat", "BIC" )
      rownames(info) <- info[, 1]
    }
    
    list(mat = t(mat), info = info, models = models, final = final, ci_test = "testIndBeta", runtime = runtime ) 
    
  }         
  
  










