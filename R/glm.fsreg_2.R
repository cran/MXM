glm.fsreg_2 <- function(target, dataset, iniset = NULL, wei = NULL, threshold = 0.05, tol = 2, heavy, robust = FALSE, ncores = 1) {
  
  ## target can be Real valued (normal), binary (binomial) or counts (poisson)
  dm <- dim(dataset) 
  if ( is.null(dm) ) {
    n <- length(target)
    p <- 1
  } else {
    n <- dm[1]  ## sample size 
    p <- dm[2]  ## number of variables
  }  
  devi <- dof <- numeric( p )  
  moda <- list()
  k <- 1   ## counter
  tool <- numeric( min(n, p) )
  threshold <- log(threshold)
  
  #########
  ## if it is binomial or poisson regression
  #########
  
  if ( is.null(iniset) ) {
    da <- 0
    pa <- 0
  } else {
    pa <- NCOL(iniset)
    da <- 1:pa
    dataset <- cbind(iniset, dataset)
  }  
  
  if ( is.matrix(target)  &  NCOL(target) == 2 )  {
    
    ci_test <- "testIndBinom"
    y <- target[, 1]
    wei <- target[, 2]
    ywei <- y / wei
    
    runtime <- proc.time()
    
    devi = dof = numeric(p)
    #if ( robust == FALSE ) {
    if ( pa == 0 ) {
      mi <- glm( ywei ~ 1, weights = wei, family = binomial, y = FALSE, model = FALSE )
      do <- 1
      ini <- mi$deviance  ## residual deviance
    } else  
      mi <- glm(ywei ~., data = as.data.frame( iniset ), weights = wei, family = binomial, y = FALSE, model = FALSE )
      do <- length( coef(mi) )
      ini <- mi$deviance  ## residual deviance
    #} else {
    #  ini = robust::glmRob( target ~ 1, family = oiko, maxit = maxit )$deviance  ## residual deviance
    #}
    
    if (ncores <= 1) {
      #if ( robust == FALSE ) {  ## Non robust
      for (i in 1:p) {
        mi <- glm( ywei ~ . , as.data.frame( dataset[, c(da, pa + i)] ), weights = wei, family = binomial, y = FALSE, model = FALSE )
        devi[i] <- mi$deviance
        dof[i] = length( coef( mi ) ) 
      }
      
      #} else {  ## Robust
      #  for (i in 1:p) { 
      #    mi <- robust::glmRob( target ~ dataset[, i], family = oiko, maxit = maxit )
      #    devi[i] <- mi$deviance
      #    dof[i] = length( coef( mi ) )          
      #  }
      #}
      
      stat = ini - devi
      pval = pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
      
    } else {
      #if ( robust == FALSE ) {  ## Non robust
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      mod <- foreach( i = 1:p, .combine = rbind) %dopar% {
        ww <- glm( ywei ~., data = as.data.frame( dataset[, c(da, pa + i)] ), weights = wei, family = binomial )
        return( c( ww$deviance, length( coef( ww ) ) ) )
      }
      
      stopCluster(cl)
      
      #       } else {  ## Robust
      #         cl <- makePSOCKcluster(ncores)
      #         registerDoParallel(cl)
      #         mata <- matrix(0, p, 2)
      #         mod <- foreach( i = 1:p, .combine = rbind, .export = "glmRob", .packages = "robust" ) %dopar% {
      #           ww <- robust::glmRob( target ~ dataset[, i], family = oiko, maxit = maxit  )
      #           mata[i, ] <- c( ww$deviance, length( coef( ww ) )  )
      #         }
      # 
      #         stopCluster(cl)
      #       }
      
      stat <- ini - mod[, 1]
      pval <- pchisq( stat, mod[, 2] - 1, lower.tail = FALSE, log.p = TRUE )
    }
    
    mat <- cbind(1:p, pval, stat) 
    
    colnames(mat)[1] <- "variables"
    rownames(mat) <- 1:p
    
    sel <- which.min(pval)
    info <- matrix( c( 1e300, 0, 0 ), ncol = 3 )
    sela <- sel
    
    if ( mat[sel, 2] < threshold ) {
      info[k, ] <- mat[sel, ]
      mat <- mat[-sel, , drop = FALSE] 
      #if ( robust == FALSE ) {
      mi <- glm( ywei ~., data = as.data.frame( dataset[, c(da, sel) ] ), weights = wei, family = binomial, y = FALSE, model = FALSE )
      tool[k] <- BIC( mi )
      #} else {
      #  mi <- robust::glmRob( target ~ dataset[, sel], family = oiko, maxit = maxit )
      #  tool[1] <- mi$deviance + length( coef( mi ) ) * log(n)
      #}
      moda[[ k ]] <- mi
    }
    
    
    ############
    ###       k equals 2
    ############ 
    
    
    if ( info[k, 2] < threshold  &  nrow(mat) > 0 ) {
      
      k <- k + 1
      pn <- p - k + 1   
      ini <- moda[[ 1 ]]$deviance  ## residual deviance
      do <- length( coef( moda[[ 1 ]]  ) ) 
      devi <- dof <- numeric( pn )  
      
      if ( ncores <= 1 ) {
        #if ( robust == FALSE ) {  ## Non robust
        for ( i in 1:pn ) {
          ww <- glm( ywei ~., data = as.data.frame( dataset[, c(da, sela, mat[pa + i, 1]) ] ), weights = wei, family = binomial, y = FALSE, model = FALSE )
          devi[i] <- ww$deviance
          dof[i] <- length( coef( ww ) )          
        }
        
        # } else {  ## Robust
        #   for ( i in 1:pn ) {
        #     ww <- robust::glmRob( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), family = oiko, maxit = maxit )
        #     devi[i] <- ww$deviance
        #     dof[i] = length( coef( ww ) ) 
        #   }   
        # }
        
        stat <- ini - devi
        pval <- pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
        
      } else {
        
        #if ( robust == FALSE ) {  ## Non robust
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
          ww <- glm( ywei ~., data = as.data.frame( dataset[, c(da, sela, mat[pa + i, 1]) ] ), weights = wei, family = binomial )
          return( c( ww$deviance, length( coef( ww ) ) ) )
        }
        
        stopCluster(cl)
        
        # } else {  ## Robust
        #   cl <- makePSOCKcluster(ncores)
        #   registerDoParallel(cl)
        #   mata = matrix(0, pn, 2)  
        #   mod <- foreach( i = 1:pn, .combine = rbind, .export = "glmRob", .packages = "robust" ) %dopar% {
        #     ww <- robust::glmRob( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), family = oiko, maxit = maxit )
        #     mata[i, ] <- c( ww$deviance, length( coef( ww ) ) )
        #   }
        # 
        #   stopCluster(cl)
        # }
        
        stat <- ini - mod[, 1]
        pval <- pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
        
      }
      
      mat[, 2:3] <- cbind(pval, stat)
      ina <- which.min(mat[, 2])
      sel <- mat[ina, 1]    
      
      if ( mat[ina, 2] < threshold ) {
        #if ( robust == FALSE ) {
        ma <- glm( ywei ~., data=as.data.frame( dataset[, c(da, sela, sel) ] ), weights = wei, family = binomial, y = FALSE, model = FALSE )
        tool[k] <- BIC( ma )
        #} else {
        #  ma <- robust::glmRob( target ~  dataset[, sela] + dataset[, sel], family = oiko, maxit = maxit )
        #  tool[k] <- ma$deviance + length( coef( ma ) ) * log(n)
        #}
        
        if ( tool[ k - 1 ] - tool[ k ] <= tol ) {
          info <- info
          
        } else { 
          info <- rbind(info, c( mat[ina, ] ) )
          sela <- info[, 1]
          mat <- mat[-ina , , drop = FALSE] 
          moda[[ k ]] <- ma
        }
        
      } else   info <- info
      
    }
    ############
    ###       k greater than 2
    ############ 
    if ( nrow(info) > 1  &  nrow(mat) > 0 ) {
      while ( info[k, 2] < threshold &  k < n - 15  &  tool[ k - 1 ] - tool[ k ] > tol  &  nrow(mat) > 0  )  {
        
        ini <- moda[[ k ]]$deviance  ## residual deviance
        do <- length( coef( moda[[ k ]]  ) ) 
        k <- k + 1   
        pn <- p - k  + 1
        devi <- dof <- numeric( pn )  
        
        if (ncores <= 1) {  
          #if ( robust == FALSE ) {  ## Non robust
          for ( i in 1:pn ) {
            ma <- glm( ywei ~., data = as.data.frame( dataset[, c(da, sela, mat[pa + i, 1] ) ] ), weights = wei, family = binomial, y = FALSE, model = FALSE )
            devi[i] <- ma$deviance
            dof[i] <- length( coef( ma ) ) 
          }
          # } else {  ## Robust
          #   for ( i in 1:pn ) {
          #     ma <- robust::glmRob( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), family = oiko, maxit = maxit )
          #     devi[i] <- ma$deviance
          #     dof[i] = length( coef( ma ) ) 
          #   }
          # }
          stat <- ini - devi
          pval <- pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
          
        } else {
          #if ( robust == FALSE ) {  ## Non robust
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
            ww <- glm( ywei ~., data = as.data.frame( dataset[, c(da, sela, mat[pa + i, 1]) ] ), weights = wei, family = binomial, y = FALSE, model = FALSE )
            return( c( ww$deviance, length( coef( ww ) ) ) )
          }
          
          stopCluster(cl)
          
          # } else {  ## Robust
          #   cl <- makePSOCKcluster(ncores)
          #   registerDoParallel(cl)
          #   mata = matrix(0, pn, 2)  
          #   mod <- foreach( i = 1:pn, .combine = rbind, .export = "glmRob", .packages = "robust" ) %dopar% {
          #     ww <- robust::glmRob( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), family = oiko, maxit = maxit )
          #     mata[i, ] <- c( ww$deviance, length( coef( ww ) ) )
          #   }
          # 
          #   stopCluster(cl)
          # }
          
          stat <- ini - mod[, 1]
          pval <- pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
          
        }
        
        mat[, 2:3] <- cbind(pval, stat)
        
        ina <- which.min(mat[, 2])
        sel <- mat[ina, 1]    
        
        if ( mat[ina, 2] < threshold ) {
          #if ( robust == FALSE ) {
          ma <- glm( ywei ~., data = as.data.frame( dataset[, c(da, sela, sel) ] ), weights = wei, family = binomial, y = FALSE, model = FALSE )
          tool[k] <- BIC( ma )
          #} else {
          #  ma <- robust::glmRob( target ~., data = as.data.frame( dataset[, c(sela, sel) ] ), family = oiko, maxit = maxit )
          #  tool[k] <- ma$deviance + length( coef( ma ) ) * log(n)
          #} 
          if ( tool[ k - 1 ] - tool[ k  ] < tol ) {
            info <- rbind(info, c( 1e300, 0, 0 ) )
            
          } else { 
            info <- rbind( info, mat[ina, ] )
            sela <- info[, 1]
            mat <- mat[-ina , , drop = FALSE]
            moda[[ k ]] <- ma
          } 
          
        } else   info <- rbind(info, c( 1e300, 0, 0 ) )
        
      }
      
    } 
    
    runtime <- proc.time() - runtime
    
    d <- length(moda)
    final <- NULL
      
      if ( d >= 1 ) {
        final <- glm( ywei ~., data = as.data.frame( dataset[, c(da, sela) ] ), weights = wei, family = binomial, y = FALSE, model = FALSE )
        info <- info[1:d, , drop = FALSE]
        info <- cbind( info, tool[ 1:d ] ) 
        colnames(info) <- c( "variables", "log.p-value", "stat", "BIC" )
        rownames(info) <- info[, 1]
      }
    
    result <- list(mat = t(mat), info = info, final = final, runtime = runtime ) 
   
    #############################
    #############################
    
  } else {
    
    ####################
    ### Logistic or poisson regression
    ####################
    
  if ( length( unique(target) ) == 2 ) {
    oiko <- binomial(logit)  ## binomial regression
    ci_test <- "testIndLogistic"
  } else  {
    ci_test <- "testIndLogistic"
    oiko <- poisson(log)  ## poisson regression
  }
    
  runtime <- proc.time()
  
  devi = dof = numeric(p)
  
  #if ( robust == FALSE ) {
  if ( pa == 0 ) {
  
    if ( !heavy ) {
      ini <- glm( target ~ 1, family = oiko, weights = wei, y = FALSE, model = FALSE )$deviance  ## residual deviance
      do <- 1
	  } else {
	    ci_test <- "testIndSpeedglm"
	    ini <- speedglm::speedglm( target ~ 1, data = data.frame(dataset), family = oiko, weights = wei )$deviance  ## residual deviance
	    do <- 1
	  }   
	
  } else {
    if ( !heavy ) {
      mi <- glm(target ~., data = data.frame( iniset ), family = oiko, y = FALSE, model = FALSE )
      ini <- mi$deviance  ## residual deviance
	  } else {
	    ci_test <- "testIndSpeedglm"
	    mi <- speedglm::speedglm(target ~., data = data.frame( iniset ), family = oiko, weights = wei )
	    ini <- mi$deviance  ## residual deviance
	  }
	  do <- length( coef(mi) )
  }
  
  #} else {
  #  ini = robust::glmRob( target ~ 1, family = oiko, maxit = maxit )$deviance  ## residual deviance
  #}
  
  if (ncores <= 1) {
    #if ( robust == FALSE ) {  ## Non robust
	  if ( !heavy  ) {
      for (i in 1:p) {
        mi <- glm( target ~ ., data.frame( dataset[, c(da, pa + i)] ), family = oiko, weights= wei, y = FALSE, model = FALSE )
        devi[i] <- mi$deviance
        dof[i] <- length( coef( mi ) ) 
      }
    
	} else {
	  for (i in 1:p) {
        mi <- speedglm::speedglm( target ~ ., data.frame( dataset[, c(da, pa + i)] ), family = oiko, weights = wei )
        devi[i] <- mi$deviance
        dof[i] <- length( coef( mi ) ) 
      }
	  
	}
    #} else {  ## Robust
    #  for (i in 1:p) { 
    #    mi <- robust::glmRob( target ~ dataset[, i], family = oiko, maxit = maxit )
    #    devi[i] <- mi$deviance
    #    dof[i] = length( coef( mi ) )          
    #  }
    #}
    
    stat <- ini - devi
    pval <- pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
    
  } else {
    #if ( robust == FALSE ) {  ## Non robust
	if ( !heavy ) { 
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      mod <- foreach( i = 1:p, .combine = rbind) %dopar% {
        ww <- glm( target ~., data = data.frame( dataset[, c(da, pa + i)] ), family = oiko, weights = wei, y = FALSE, model = FALSE )
        return( c( ww$deviance, length( coef(ww) )  ) )
      }
    
      stopCluster(cl)
    
	} else {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      mod <- foreach( i = 1:p, .combine = rbind, .export = "speedglm", .packages = "speedglm") %dopar% {
        ww <- speedglm::speedglm( target ~., data = data.frame( dataset[, c(da, pa + i)] ), family = oiko, weights = wei )
        return( c( 2 * logLik(ww), length( coef(ww) )  ) )
      }
    
      stopCluster(cl)	
	}
	
    #       } else {  ## Robust
    #         cl <- makePSOCKcluster(ncores)
    #         registerDoParallel(cl)
    #         mata <- matrix(0, p, 2)
    #         mod <- foreach( i = 1:p, .combine = rbind, .export = "glmRob", .packages = "robust" ) %dopar% {
    #           ww <- robust::glmRob( target ~ dataset[, i], family = oiko, maxit = maxit  )
    #           mata[i, ] <- c( ww$deviance, length( coef( ww ) )  )
    #         }
    # 
    #         stopCluster(cl)
    #       }
    
    stat <- ini - mod[, 1]
    pval <- pchisq( stat, mod[, 2] - 1, lower.tail = FALSE, log.p = TRUE )
  }
  
  mat <- cbind(1:p, pval, stat) 
  
  colnames(mat)[1] <- "variables"
  rownames(mat) <- 1:p
  
  sel <- which.min(pval)
  info <- matrix( c( 1e300, 0, 0 ), ncol = 3 )
  sela <- sel
  
  if ( mat[sel, 2] < threshold ) {
    info[k, ] <- mat[sel, ]
    mat <- mat[-sel, , drop= FALSE] 

    #if ( robust == FALSE ) {
	   if ( !heavy ) {
      mi <- glm( target ~., data = as.data.frame( dataset[, c(da, sel) ] ), family = oiko, weights = wei, y = FALSE, model = FALSE )
      tool[k] <- BIC( mi )
	  } else {
	    mi <- speedglm::speedglm( target ~., data = as.data.frame( dataset[, c(da, sel) ] ), family = oiko, weights = wei )
      tool[k] <- BIC(mi)
  	}
    #} else {
    #  mi <- robust::glmRob( target ~ dataset[, sel], family = oiko, maxit = maxit )
    #  tool[1] <- mi$deviance + length( coef( mi ) ) * log(n)
    #}
    moda[[ k ]] <- mi
  }
  
  
  ############
  ###       k equals 2
  ############ 
  
  
  if ( info[k, 2] < threshold  &  nrow(mat) > 0 ) {
    
    k <- k + 1
    pn <- p - k + 1   
    
    ini <- mi$deviance  ## residual deviance
    do <- length( coef( mi ) ) 
    
    if ( ncores <= 1 ) {
      devi <- dof <- numeric(pn)
      #if ( robust == FALSE ) {  ## Non robust
	  if ( !heavy ) {
        for ( i in 1:pn ) {
          ww <- glm( target ~., data = as.data.frame( dataset[, c(da, sela, mat[pa + i, 1]) ] ), family = oiko, weights = wei, y = FALSE, model = FALSE )
          devi[i] <- ww$deviance
          dof[i] <- length( coef( ww ) )          
        }
		
	  } else {
        for ( i in 1:pn ) {
          ww <- speedglm::speedglm( target ~., data = as.data.frame( dataset[, c(da, sela, mat[pa + i, 1]) ] ), family = oiko, weights = wei )
          devi[i] <- ww$deviance
          dof[i] <- length( coef( ww ) )          
        }	  
	  }	
      
      # } else {  ## Robust
      #   for ( i in 1:pn ) {
      #     ww <- robust::glmRob( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), family = oiko, maxit = maxit )
      #     devi[i] <- ww$deviance
      #     dof[i] = length( coef( ww ) ) 
      #   }   
      # }
      
      stat <- ini - devi
      pval <- pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
      
    } else {
      
      #if ( robust == FALSE ) {  ## Non robust
	  if ( !heavy ) {
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
          ww <- glm( target ~., data = as.data.frame( dataset[, c(da, sela, mat[pa + i, 1]) ] ), family = oiko, weights = wei, y = FALSE, model = FALSE )
          return( c( ww$deviance, length( coef( ww ) ) ) )
        }
      
        stopCluster(cl)
	  
      } else {
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mod <- foreach( i = 1:pn, .combine = rbind, .export = "speedglm", .packages = "speedglm") %dopar% {
          ww <- speedglm::speedglm( target ~., data = as.data.frame( dataset[, c(da, sela, mat[pa + i, 1]) ] ), family = oiko, weights = wei )
          return( c( ww$deviance, length( coef( ww ) ) ) )
        }
      
        stopCluster(cl)	  
	  }
	  
	  
      # } else {  ## Robust
      #   cl <- makePSOCKcluster(ncores)
      #   registerDoParallel(cl)
      #   mata = matrix(0, pn, 2)  
      #   mod <- foreach( i = 1:pn, .combine = rbind, .export = "glmRob", .packages = "robust" ) %dopar% {
      #     ww <- robust::glmRob( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), family = oiko, maxit = maxit )
      #     mata[i, ] <- c( ww$deviance, length( coef( ww ) ) )
      #   }
      # 
      #   stopCluster(cl)
      # }
      
      stat <- ini - mod[, 1]
      pval <- pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
      
    }
    
    mat[, 2:3] <- cbind(pval, stat)
    
    ina <- which.min(mat[, 2])
    sel <- mat[ina, 1]    
    
    if ( mat[ina, 2] < threshold ) {
      #if ( robust == FALSE ) {
	    if ( !heavy )  {
        ma <- glm( target ~., data=as.data.frame( dataset[, c(da, sela, sel) ] ), family = oiko, weights = wei, y = FALSE, model = FALSE )
        tool[k] <- BIC( ma )
	    } else {
        ma <- speedglm::speedglm( target ~., data=as.data.frame( dataset[, c(da, sela, sel) ] ), family = oiko, weights = wei )
        tool[k] <-  BIC(ma)
	    }	
      #} else {
      #  ma <- robust::glmRob( target ~  dataset[, sela] + dataset[, sel], family = oiko, maxit = maxit )
      #  tool[k] <- ma$deviance + length( coef( ma ) ) * log(n)
      #}
      
      if ( tool[ k - 1 ] - tool[ k ] <= tol ) {
        info <- rbind(info, c( 1e300, 0, 0 ) )
        
      } else { 
        info <- rbind(info, c( mat[ina, ] ) )
        sela <- info[, 1]
        mat <- mat[-ina , , drop = FALSE] 
        moda[[ k ]] <- ma
      }
      
    } else    info <- info
    
  }
  
  ############
  ###       k greater than 2
  ############ 
   
  if ( nrow(info) > 1  &  nrow(mat) > 0 ) {
    while ( ( info[k, 2] < threshold ) &  ( k < n ) & ( tool[ k - 1 ] - tool[ k ] > tol ) & ( nrow(mat) > 0 ) )  {
      
      ini <- moda[[ k ]]$deviance  ## residual deviance
      do <- length( coef( moda[[ k ]]  ) ) 
      
      k <- k + 1   
      pn <- p - k  + 1
      
      if (ncores <= 1) {  
        devi <- dof <- numeric(pn) 
        #if ( robust == FALSE ) {  ## Non robust
		    if ( !heavy ) {
          for ( i in 1:pn ) {
            ma <- glm( target ~., data = as.data.frame( dataset[, c(da, sela, mat[pa + i, 1] ) ] ), family = oiko, weights = wei, y = FALSE, model = FALSE )
            devi[i] <- ma$deviance
            dof[i] <- length( coef( ma ) ) 
          }
		
        } else {
          for ( i in 1:pn ) {
            ma <- speedglm:speedglm( target ~., data = as.data.frame( dataset[, c(da, sela, mat[pa + i, 1] ) ] ), family = oiko, weights = wei )
            devi[i] <- ma$deviance
            dof[i] <- length( coef( ma ) ) 
          }		
		    }
        # } else {  ## Robust
        #   for ( i in 1:pn ) {
        #     ma <- robust::glmRob( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), family = oiko, maxit = maxit )
        #     devi[i] <- ma$deviance
        #     dof[i] = length( coef( ma ) ) 
        #   }
        # }
        stat <- ini - devi
        pval <- pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
        
      } else {
        #if ( robust == FALSE ) {  ## Non robust
		    if ( !heavy )  {
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          devi <- dof <- numeric(pn)
          mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
            ww <- glm( target ~., data = as.data.frame( dataset[, c(da, sela, mat[pa + i, 1]) ] ), family = oiko, weights = wei, y = FALSE, model = FALSE )
            return( c( ww$deviance, length( coef( ww ) ) ) )
          }
        
          stopCluster(cl)
        } else {
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          mod <- foreach( i = 1:pn, .combine = rbind ,.export = "speedglm", .packages = "speedglm") %dopar% {
            ww <- speedglm( target ~., data = as.data.frame( dataset[, c(da, sela, mat[pa + i, 1]) ] ), family = oiko, weights = wei )
            return( c( ww$deviance, length( coef( ww ) ) ) )
          }
        
          stopCluster(cl)		
	    	}
        # } else {  ## Robust
        #   cl <- makePSOCKcluster(ncores)
        #   registerDoParallel(cl)
        #   mata = matrix(0, pn, 2)  
        #   mod <- foreach( i = 1:pn, .combine = rbind, .export = "glmRob", .packages = "robust" ) %dopar% {
        #     ww <- robust::glmRob( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), family = oiko, maxit = maxit )
        #     mata[i, ] <- c( ww$deviance, length( coef( ww ) ) )
        #   }
        # 
        #   stopCluster(cl)
        # }
        stat <- ini - mod[, 1]
        pval <- pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
      }
      
      mat[, 2:3] <- cbind(pval, stat)
      ina <- which.min(mat[, 2])
      sel <- mat[ina, 1]    
      
      if ( mat[ina, 2] < threshold ) {
        #if ( robust == FALSE ) {
		    if ( !heavy ) {
          ma <- glm( target ~., data = as.data.frame( dataset[, c(da, sela, sel) ] ), family = oiko, weights = wei, y = FALSE, model = FALSE )
          tool[k] <- BIC( ma )
        } else {
          ma <- speedglm::speedglm( target ~., data = as.data.frame( dataset[, c(da, sela, sel) ] ), family = oiko, weights = wei )
          tool[k] <-  BIC(ma)	
		    }
		#} else {
        #  ma <- robust::glmRob( target ~., data = as.data.frame( dataset[, c(sela, sel) ] ), family = oiko, maxit = maxit )
        #  tool[k] <- ma$deviance + length( coef( ma ) ) * log(n)
        #} 
        if ( tool[ k - 1 ] - tool[ k  ] <= tol ) {
          info <- rbind(info, c( 1e300, 0, 0 ) )
          
        } else { 
          info <- rbind( info, mat[ina, ] )
          sela <- info[, 1]
          mat <- mat[-ina , , drop = FALSE]
          moda[[ k ]] <- ma
        } 
        
      } else   info <- rbind(info, c( 1e300, 0, 0 ) )
      
    }
    
  } 
  
  runtime <- proc.time() - runtime
  
  d <- length(moda)

  if ( d == 0 ) {
    if ( heavy ) {
      final <- speedglm::speedglm( target ~., data = as.data.frame( iniset ), weights = wei, family = oiko, y = FALSE, model = FALSE )
    }  else  final <- glm( target ~., data = as.data.frame( iniset ), weights = wei, family = oiko, y = FALSE, model = FALSE )
    
  } else {
    
    if ( !heavy ) {
      final <- glm( target ~., data = as.data.frame( dataset[, c(da, sela) ] ), family = oiko, weights = wei, y = FALSE, model = FALSE )
	} else  final <- speedglm::speedglm( target ~., data = as.data.frame( dataset[, c(da, sela) ] ), family = oiko, weights = wei )
    info <- info[1:d, , drop = FALSE]
    info <- cbind( info, tool[ 1:d ] ) 
    colnames(info) <- c( "variables", "log.p-value", "stat", "BIC" )
    rownames(info) <- info[, 1]
  }    
  
  result <- list(runtime = runtime, mat = t(mat), info = info, ci_test = ci_test, final = final ) 

  }
  
  result   
}    










