lm.fsreg_2 <- function(target, dataset, iniset = NULL, threshold = 0.05, wei = NULL, stopping = "BIC", tol = 2, heavy = FALSE, robust = FALSE, ncores = 1 ) {
  
  threshold = log(threshold)
  
  p <- dim(dataset)[2]  ## number of variables
  pval <- stat <- dof <- numeric( p )  
  moda <- list()
  
    #check for NA values in the dataset and replace them with the variable median or the mode
  if(any(is.na(dataset)) == TRUE)
  {

  
    #dataset = as.matrix(dataset);
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    
    if (class(dataset) == "matrix")  {
    
       dataset = apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x)} ) 
              
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
  
  ## is there already an initial set of variables to start with?
  if ( is.null(iniset) ) {
    da <- 0
    pa <- 0
    
  } else {
    pa <- NCOL(iniset)
    da <- 1:pa
    dataset <- cbind(iniset, dataset)
  }  

  
  n <- length(target)  ## sample size
  con <- log(n)
  tool <- numeric( length( min(n, p) ) )
  info <- matrix( c( 1e300, 0, 0 ), ncol = 3 )
  
  k <- 1   ## counter k is 1, step 1
  
  runtime <- proc.time()
  
  if (ncores <= 1) {
    
    if ( robust == FALSE ) {  ## Non robust
	
	  if ( heavy == FALSE ) {
	      ci_test <- "testsIndReg"
        for (i in 1:p) {
          ww = lm( target ~., data = as.data.frame( dataset[, c(da, pa + i)] ), weights = wei )
          tab = anova(ww)
          stat[i] = tab[pa + k, 4] 
          df1 = tab[pa + k, 1]   ;  df2 = tab[pa + k + 1, 1]
          pval[i] = pf( stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE )
        } 
		
	  } else {
	    ci_test <- "testIndSpeedglm"
	  	fit1 <- speedglm::speedlm( target ~., data = as.data.frame(iniset), weights = wei )
		  fit1rss <- fit1$RSS 
      d1 <- length( coef(fit1) )
		  for (i in 1:pn) {
        fit2 = speedglm::speedlm( target ~., data = as.data.frame( dataset[, c(da, pa + i)] ), weights= wei )
        d2 = length( coef(fit2) )
        df1 = d2 - d1
        df2 = n - d2
        stat[i] = (fit1rss - fit2$RSS)/df1 / ( fit2$RSS /df2 )
        pval[i] = pf( stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE )
      } 	
	  
	  }	
      mat <- cbind(1:p, pval, stat)
      
    } else {  ## Robust
      
      for (i in 1:p) {
        ww = MASS::rlm( target ~., data = as.data.frame( dataset[, c(da, pa + i)] ), maxit = 2000, weights = wei ) 
        stat[i] = 2 * as.numeric( logLik(ww) )
        dof[i] = length( coef(ww) )
      }
	  
      fit0 = MASS::rlm(target ~ ., data = as.data.frame(iniset), maxit = 2000, weights = wei) 
      stat0 = 2 * as.numeric( logLik(fit0) )
      difa = abs( stat - stat0 )
      pval = pchisq(difa, dof - 1, lower.tail = FALSE, log.p = TRUE)
      mat <- cbind(1:p, pval, difa)
    }
    
  } else {
    
    if ( robust == FALSE ) {  ## Non robust
	
	  if ( heavy == FALSE ) {
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mat <- matrix(0, p, 2)
        mod <- foreach( i = 1:p, .combine = rbind ) %dopar% {
          ww <- lm( target ~., data = as.data.frame( dataset[, c(da, pa + i)] ), weights = wei )
          tab <- anova( ww )
          stat[i] <- tab[pa + k, 4] 
          df1 <- tab[pa + k, 1]   ;  df2 = tab[pa + k + 1, 1]
          pval[i] <- pf( stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE )
          mat[i, ] <- c(pval[i], stat[i]) 
        }
		
        stopCluster(cl)
		
	  } else {
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
	      fit1 <- speedglm::speedlm( target ~., data = as.data.frame(iniset), weights = wei )
		    fit1rss <- fit1$RSS 
        d1 <- length( coef(fit1) )
        mat <- matrix(0, p, 2)
        mod <- foreach( i = 1:p, .combine = rbind, .export = "speedlm", .packages = "speedglm" ) %dopar% {
          fit2 <- speedlm( target ~., data = as.data.frame( dataset[, c(da, pa + i)] ), weights = wei )
          d2 = length( coef(fit2) )
          df1 = d2 - d1
          df2 = n - d2
          stat[i] = (fit1rss - fit2$RSS)/df1 / ( fit2$RSS /df2 )
          pval[i] = pf( stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE )
          mat[i, ] <- c(pval[i], stat[i]) 
        }
		
        stopCluster(cl)	 	  
	  }
	       
    } else {  ## Robust
      
      fit0 = MASS::rlm(target ~., data = as.data.frame(iniset), maxit = 2000, weights = wei) 
      stat0 = 2 * logLik(fit0)
      
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      mat <- matrix(0, p, 2)
      mod <- foreach( i = 1:p, .combine = rbind, .export = "rlm", .packages = "MASS" ) %dopar% {
        ww = rlm( target ~., data = as.data.frame( dataset[, c(da, pa + i)] ), maxit = 2000, weights = wei ) 
        mat[i, ] = c( 2 * as.numeric( logLik(ww) ), length( coef(ww) ) )
      }
      stopCluster(cl)
      
      difa = abs( mod[, 1] - stat0 )
      pval = pchisq(difa, mod[, 2] - 1, lower.tail = FALSE, log.p = TRUE)
      mod = cbind( pval, difa)
    }
    
    mat <- cbind(1:p, mod)      
    
  }
  
  colnames(mat) <- c( "variables", "p-value", "stat" )
  rownames(mat) <- 1:p
  sel <- which.min(mat[, 2])
  sela <- sel
  
  if ( mat[sel, 2] < threshold ) {
    
    info[k, ] <- mat[sel, ]
    mat <- mat[-sel, ] 
    if ( !is.matrix(mat) ) {
      mat <- matrix(mat, ncol = 3) 
    }
    mat <- mat[ order( mat[, 2] ), ] 
    
    if ( stopping == "adjrsq" ) {
      
      if ( robust == FALSE ) {
	  
	      if ( heavy == TRUE ) {
          mi = lm( target ~., data = as.data.frame( dataset[, c(da, sel) ] ), weights = wei )
          tool[k] <- as.numeric( summary( mi )[[ 9 ]] )
		    } else {
		      mi = speedglm::speedlm( target ~., data = as.data.frame( dataset[, c(da, sel) ] ), weights = wei )
          tool[k] <- as.numeric( summary( mi )[[ 11 ]] )
		    }
		
      } else {
        mi = MASS::rlm( target ~., data = as.data.frame( dataset[, c(da, sel) ] ), maxit = 2000, weights = wei )
        r2 = cor( target, fitted(mi) )^2
        tool[k] = 1 - (1 - r2) * (n - 1) / ( n - length( coef(mi) ) - 1 )
      } 
      
    } else if ( stopping  == "BIC" ) { 
      if ( robust == FALSE ) {
	      if ( heavy == FALSE ) {
          mi = lm( target ~., data = as.data.frame( dataset[, c(da, sel) ] ) )
          tool[k] <- BIC( mi )
		    } else {
          mi = speedglm::speedlm( target ~., data = as.data.frame( dataset[, c(da, sel) ] ), weights = wei )
          tool[k] <-  - 2 * mi$logLik + length( coef(mi) ) * con 
	   	  }
		
      } else {
        mi = MASS::rlm( target ~., data=as.data.frame( dataset[, c(da, sel) ] ), maxit = 2000, weights = wei )
        tool[k] <-  BIC(mi)
      }
    }
    moda[[ k ]] <- mi
    
  }
  
  
  ######
  ###### k equal to 2
  ######
  
  if ( info[k, 2] < threshold  &  nrow(mat) > 0 ) {
    
    k <- k + 1
    pn <- p - k + 1   
    
    if ( ncores <= 1 ) {
	
      if ( robust == FALSE ) { 
	  
        if ( heavy == FALSE ) {
          for (i in 1:pn) {
            ww = lm( target ~., data = as.data.frame( dataset[, c(da, sel, mat[pa + i, 1]) ] ), weights = wei )
            tab = anova( ww )
            mat[i, 3] = tab[pa + k, 4] 
            df1 = tab[pa + k, 1]   ;  df2 = tab[pa + k + 1, 1]
            mat[i, 2] = pf( mat[i, 3], df1, df2, lower.tail = FALSE, log.p = TRUE )
          }
		  
        } else {
		      fit1rss <- mi$RSS 
          d1 <- length( coef(fit1) )
		      for (i in 1:pn) {
            fit2 = speedglm::speedlm( target ~., data = as.data.frame( dataset[, c(da, sel, mat[pa + i, 1]) ] ), weights= wei )
            d2 = length( coef(fit2) )
            df1 = d2 - d1
            df2 = n - d2
            mat[i, 3] = (fit1rss - fit2$RSS)/df1 / ( fit2$RSS /df2 )
            mat[i, 2] = pf( mat[i, 3], df1, df2, lower.tail = FALSE, log.p = TRUE )
		     }	  
		   }
		
      } else {
        do = 2 * as.numeric( logLik( mi ) )
        fr = length( coef( mi ) )
        sta = dof = numeric(pn)
        
        for (i in 1:pn) {
          ww = MASS::rlm( target ~., data = as.data.frame( dataset[, c(da, sel, mat[pa + i, 1]) ] ), maxit = 2000, weights = wei )
          sta[i] = 2 * as.numeric( logLik(ww) )
          dof[i] = length( coef(ww) )
        } 
        mat[, 3] = abs( sta - do )
        mat[, 2] = pchisq(mat[, 3], dof - fr, lower.tail = FALSE, log.p = TRUE)
      }
      
    } else {
      
      if ( robust == FALSE ) {  ## Non robust
	  
	      if ( heavy == FALSE ) {
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          stat <- pval <- numeric(pn) 
          mata <- matrix(0, pn, 2)
          mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
            ww <- lm( target ~., data = as.data.frame( dataset[, c(da, sel, mat[pa + i, 1]) ] ), weights = wei )
            tab <- anova( ww )
            stat[i] <- tab[pa + k, 4] 
            df1 <- tab[pa + k, 1]   ;  df2 = tab[pa + k + 1, 1]
            pval[i] <- pf( stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE )
            mata[i, ] <- c(pval[i], stat[i]) 
          }
          stopCluster(cl)
		  
        } else {
		      cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          fit1rss <- mi$RSS 
          stat <- pval <- numeric(pn) 
          mata <- matrix(0, pn, 2)
          mod <- foreach( i = 1:pn, .combine = rbind, .export = "speedlm", .packages = "speedglm") %dopar% {
            fit2 <- speedlm( target ~., data = as.data.frame( dataset[, c(da, sel, mat[pa + i, 1]) ] ), weights = wei )
            d2 = length( coef(fit2) )
            df1 = d2 - d1
            df2 = n - d2
            stat[i] = (fit1rss - fit2$RSS)/df1 / ( fit2$RSS /df2 )
            pval[i] = pf( mat[i, 3], df1, df2, lower.tail = FALSE, log.p = TRUE )
            mata[i, ] <- c( pval[i], stat[i] ) 
          }
          stopCluster(cl)
		   }
		
      } else {  ## Robust
        do = 2 * as.numeric( logLik( mi ) )
        fr = length( coef( mi ) )
        
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mata <- matrix(0, pn, 2)
        mod <- foreach( i = 1:pn, .combine = rbind, .export = c("rlm"), .packages = "MASS" ) %dopar% {
          ww <- rlm( target ~., data = as.data.frame( dataset[, c(da, sel, mat[pa + i, 1]) ] ), maxit = 2000, weights = wei )
          mata[i, ] = c( 2 * as.numeric( logLik(ww) ), length( coef(ww) ) )
        }
        stopCluster(cl)
        
        difa = abs( mod[, 1] - do )
        pval = pchisq(difa, mod[, 2] - fr, lower.tail = FALSE, log.p = TRUE)
        mod = cbind( pval, difa)
      } 
      mat <- cbind(mat[, 1], mod)   
      
    }
    
  }
  
  ina <- which.min(mat[, 2])
  sel <- mat[ina, 1]     
  
  if ( stopping == "adjrsq" ) {
  
    if ( mat[ina, 2] < threshold ) {
      if ( robust == FALSE ) {
	    if ( heavy == FALSE ) {
          ma = lm( target ~., data=as.data.frame( dataset[, c(da, sela, sel) ] ), weights = wei )
          tool[k] <- as.numeric( summary( ma )[[ 9 ]] )
		} else {
          ma = speedglm::speedlm( target ~., data=as.data.frame( dataset[, c(da, sela, sel) ] ), weights = wei )
          tool[k] <- as.numeric( summary( ma )[[ 11 ]] )		
		}
		
      } else {        
        ma= MASS::rlm( target ~., data=as.data.frame( dataset[, c(da, sela, sel) ] ), maxit = 2000, weights = wei )
        r2 = cor( target, fitted(ma) )^2
        tool[k] = 1 - (1 - r2) * (n - 1) / ( n - length( coef(ma) ) - 1)
      }
      
      if ( tool[ k ] - tool[ k - 1 ] <= tol ) {
        info <- info
        
      } else {  
        info <- rbind(info, mat[ina, ] )
        sela <- info[, 1]
        mat <- mat[-ina , ] 
        if ( is.matrix(mat) ) {
          mat <- mat[ order( mat[, 2] ), ]
        }
        
        moda[[ k ]] <- ma
      }
      
    } else  info <- info
    
  } else if ( stopping == "BIC" ) {
  
    if ( mat[ina, 2] < threshold ) {
      if ( robust == FALSE ) {            
	    if ( heavy == FALSE ) {
          ma = lm( target ~., data = as.data.frame( dataset[, c(da, sela, sel) ] ), weights = wei )
          tool[k] <- BIC( ma )
		} else {
          ma = speedglm::speedlm( target ~., data = as.data.frame( dataset[, c(da, sela, sel) ] ), weights = wei )
          tool[k] <-  - 2 *ma$logLik + length( coef(ma ) ) * con
		}
		
      } else {
        ma = MASS::rlm(target ~ target ~., data = as.data.frame( dataset[, c(da, sela, sel) ] ), maxit = 2000, weights = wei )
        tool[k] = BIC(ma)
      }
      if ( tool[ k - 1] - tool[ k ] <= tol ) {
        info <- rbind(info, c( 1e300, 0, 0 ) )
        
      } else {  
        info <- rbind(info, mat[ina, ] )
        sela <- info[, 1]
        mat <- mat[-ina , ]
        if ( !is.matrix(mat) ) {
          mat <- matrix(mat, ncol = 3) 
        }
        mat <- mat[ order( mat[, 2] ), ]
        
        moda[[ k ]] <- ma
      }
    } else  info <- info
    
  }
  

  
  ###########
  ######   k greater than 2
  ###########
  
  
  if ( nrow(info) > 1  &  nrow(mat) > 0 ) {
    
    while ( ( info[k, 2] < threshold ) &  ( k < n ) & ( abs( tool[ k ] - tool[ k - 1 ] ) > tol ) & ( nrow(mat) > 0 ) )  {
      
      k <- k + 1   
      pn <- p - k + 1 
      
      if ( ncores <= 1 ) {
	  
        if ( robust == FALSE ) {
		
		      if ( heavy == FALSE ) {
            for ( i in 1:pn ) {
              ww <- lm( target ~., data = as.data.frame( dataset[, c(da, sela, mat[pa + i, 1]) ] ), weights = wei )
              tab <- anova( ww )
              mat[i, 3] <- tab[pa + k, 4] 
              df1 <- tab[pa + k, 1]   ;  df2 = tab[pa + k + 1, 1]
              mat[i, 2] <- pf( mat[i, 3], df1, df2, lower.tail = FALSE, log.p = TRUE )
            }
			
          } else {
		        fit1rss <- ma$RSS 
            d1 <- length( coef(fit1) )
		        for (i in 1:pn) {
              fit2 = speedglm::speedlm( target ~., data = as.data.frame( dataset[, c(da, sel, mat[pa + i, 1]) ] ), weights= wei )
              d2 = length( coef(fit2) )
              df1 = d2 - d1
              df2 = n - d2
              mat[i, 3] = (fit1rss - fit2$RSS)/df1 / ( fit2$RSS /df2 )
              mat[i, 2] = pf( mat[i, 3], df1, df2, lower.tail = FALSE, log.p = TRUE )
		        }			  
		      }
		  
        } else {
          
          do = 2 * as.numeric( logLik( moda[[ k - 1 ]] ) )
          fr = length( coef( moda[[ k - 1 ]] ) )
          sta = dof = numeric(pn)
          
          for (i in 1:pn) {
            ww = MASS::rlm( target ~., data = as.data.frame( dataset[, c(da, sela, mat[pa + i, 1]) ] ), maxit = 2000, weights = wei )
            sta[i] = 2 * as.numeric( logLik(ww) )
            dof[i] = length( coef(ww) )
          } 
          mat[, 3] = abs( sta - do )
          mat[, 2] = pchisq(mat[, 3], dof - fr, lower.tail = FALSE, log.p = TRUE)
        }
        
      } else {
        if ( robust == FALSE ) {  ## Non robust
		
		     if ( heavy == FALSE ) {
            cl <- makePSOCKcluster(ncores)
            registerDoParallel(cl)
            stat <- pval <- numeric(pn)
            mata <- matrix(0, pn, 2)
            mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
              ww <- lm( target ~., data = as.data.frame( dataset[, c(da, sela, mat[pa + i, 1]) ] ), weights = wei )
              tab <- anova(ww)
              stat[i] <- tab[pa + k, 4] 
              df1 <- tab[pa + k, 1]   ;  df2 = tab[pa + k + 1, 1]
              pval[i] <- pf( stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE )
              mata[i, ] <- c( pval[i], stat[i] ) 
            }
            stopCluster(cl)
		  
		    } else {
		      cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          fit1rss <- ma$RSS 
          stat <- pval <- numeric(pn) 
          mata <- matrix(0, pn, 2)
          mod <- foreach( i = 1:pn, .combine = rbind, .export = "speedlm", .packages = "speedglm") %dopar% {
            fit2 <- speedlm( target ~., data = as.data.frame( dataset[, c(da, sela, mat[pa + i, 1]) ] ), weights = wei )
            d2 = length( coef(fit2) )
            df1 = d2 - d1
            df2 = n - d2
            stat[i] = (fit1rss - fit2$RSS)/df1 / ( fit2$RSS /df2 )
            pval[i] = pf( mat[i, 3], df1, df2, lower.tail = FALSE, log.p = TRUE )
            mata[i, ] <- c( pval[i], stat[i] ) 
          }
          stopCluster(cl)
		    }
          
       } else {  ## Robust
          do = 2 * as.numeric( logLik( ma) )
          fr = length( coef( ma ) )
          
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          mata <- matrix(0, pn, 2)
          mod <- foreach( i = 1:pn, .combine = rbind, .export = "rlm", .packages = "MASS" ) %dopar% {
            ww <- rlm( target ~., data = as.data.frame( dataset[, c(da, sela, mat[pa + i, 1]) ] ), maxit = 2000, weights = wei )
            mata[i, ] = c( 2 * as.numeric( logLik(ww) ), length( coef(ww) ) )
          }
          stopCluster(cl)  
          difa = abs( mod[, 1] - do )
          pval = pchisq(difa, mod[, 2] - fr, lower.tail = FALSE, log.p = TRUE)
          mod = cbind( pval, difa)
        }
        
        mat <- cbind( mat[, 1], mod )   
        
      }
      
      ina <- which.min(mat[, 2])
      sel <- mat[ina, 1]   
      
      if ( stopping == "BIC" ) {
	  
        if ( mat[ina, 2] < threshold ) {
		
          if ( robust == FALSE ) {
		        if ( heavy == FALSE ){
              ma <- lm( target ~., data = as.data.frame( dataset[, c(da, sela, sel)] ), weights = wei )
              tool[k] <- BIC( ma )
			      } else {
              ma <- speedglm::speedlm( target ~., data = as.data.frame( dataset[, c(da, sela, sel)] ), weights = wei )
              tool[k] <-  - 2 *ma$logLik + length( coef(ma ) ) * con		
			      }  
			
          } else {
            # ma = robust::lmRob( target ~., data = as.data.frame( dataset[, c(sela, sel)] ), control = cont )
            ma = MASS::rlm( target ~., data = as.data.frame( dataset[, c(da, sela, sel)] ), maxit = 2000, weights = wei )
            tool[k] =  BIC(ma)
          }
          
          if ( tool[ k - 1] - tool[ k ] <= tol ) {
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
        
      } else if ( stopping == "adjrsq" ) {
	    if ( mat[ina, 2] < threshold ) {
          
          if ( robust == FALSE ) {
		    if ( heavy == FALSE ) {
              ma <- lm( target ~., data = as.data.frame( dataset[, c(da, sela, sel)] ), weights = wei )
              tool[k] <- as.numeric( summary(ma)[[ 9 ]] )
			} else {
			  ma <- speedglm::speedlm( target ~., data=as.data.frame( dataset[, c(da, sela, sel) ] ), weights = wei )
              tool[k] <- as.numeric( summary( ma )[[ 11 ]] )	
			}
			
          } else {
            ma = MASS::rlm( target ~., data = as.data.frame( dataset[, c(da, sela, sel)] ), maxit = 2000, weights = wei )
            r2 = cor(target, fitted(ma) )^2
            tool[k] = 1 - (1 - r2) * (n - 1) / ( n - length( coef(ma) ) - 1)
          }               
          if ( tool[ k ] - tool[ k - 1 ] <= tol ) {
            info <- rbind(info, c( 1e300, 0, 0 ) )
            
          } else { 
            info <- rbind( info, mat[ina, ] )
            sela <- info[, 1]
            mat <- mat[-ina , ]
            if ( is.matrix(mat) ) {
              mat <- mat[ order( mat[, 2] ), ]
            }
            
            moda[[ k ]] <- ma
          }
          
        } else   info <- rbind(info, c( 1e300, 0, 0 ) )
        
      }
    }
    
  }
  
  runtime <- proc.time() - runtime
  
  d <- length(moda)
  final <- NULL
  models <- NULL

  if ( d == 0 ) {
    if ( robust == FALSE ) {
	
	  if ( heavy == FALSE ) {
        final <- lm( target ~., data = as.data.frame( iniset ) ) 
	  } else {
	    final <- speedglm::speedlm( target ~., data = as.data.frame( iniset ) ) 
	  }
	  
	} else {
      final <- MASS::rlm( target ~., data = as.data.frame( iniset ), maxit = 2000, weights = wei ) 
	} 
    info = NULL
    
  } else {
    
  if ( d >= 1 ) {
    models <- NULL
    xx <- as.data.frame( dataset[, c(da, sela) ] )
    if ( pa == 0 ) {
      colnames(xx) <- paste("V", sela, sep = "") 
    } else {
      colnames(xx) <- c( paste("X", da, sep = ""),  paste("V", sela, sep = "") )
    }
    
    if ( d == 1 ) {
      
      xx <- as.data.frame( dataset[, c(da, sela) ] )
      if ( pa == 0 ) {
        colnames(xx) <- paste("V", sela, sep = "") 
      } else {
        colnames(xx) <- c( paste("X", da, sep = ""),  paste("V", sela, sep = "") )
      }
      
      if ( robust == FALSE ) {
	    if ( heavy == FALSE ) {
          models[[ 1 ]] <- lm( target ~., data = as.data.frame( xx ), weights = wei )
		} else {
		  models[[ 1 ]] <- speedglm::speedlm( target ~., data = as.data.frame( xx ), weights = wei )
		}
		
      } else {
        # models[[ 1]] <- final <- robust::lmRob( target ~., data = as.data.frame( xx ), control = cont )
        models[[ 1 ]] <- MASS::rlm( target ~., data = as.data.frame( xx ), maxit = 2000, weights = wei )
      }
      
    } else {
      for ( i in 1:c( pa + d ) ) {
        if ( robust == FALSE ) {
		  if ( heavy == FALSE ) {
            models[[ i ]] <- lm( target ~., data = as.data.frame( xx[, 1:i] ), weights = wei )
		  } else {
		    models[[ i ]] <- speedglm::speedlm( target ~., data = as.data.frame( xx[, 1:i] ), weights = wei )
		  }	
		  
        } else { 
          # models[[ i ]] <- lmRob( target ~., data = as.data.frame( xx[, 1:i]), control = cont )
          models[[ i ]] <- MASS::rlm( target ~., data = data.frame( xx[, 1:i] ), maxit = 2000, weights = wei )
          
        }
      }
    }
    
    final <- summary( models[[ pa + d ]] )
    
    info <- info[1:d, ]
    if ( d == 1 )  info <- matrix(info, nrow = 1)
    info <- cbind( info, tool[ 1:d ] ) 
    colnames(info) <- c( "variables", "p-value", "stat", stopping )
    rownames(info) <- info[, 1]
  }
  
  }
  
  list(mat = t(mat), info = info, models = models, final = final, ci_test <- ci_test, runtime = runtime ) 
}
