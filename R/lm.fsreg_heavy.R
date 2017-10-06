lm.fsreg_heavy <- function(target, dataset, ini = NULL, threshold = 0.05, wei = NULL, stopping = "BIC", tol = 2, ncores = 1 ) {  
  
  ###### If there is an initial set of variables do this function
  if ( !is.null(ini) ) {
     result <- lm.fsreg_2.heavy(target, dataset, iniset = ini, threshold = threshold, wei = NULL, stopping = stopping, tol = tol, ncores = ncores) 
    
  } else {  ## else do the classical forward regression
    
    threshold <- log(threshold)
    p <- dim(dataset)[2]  ## number of variables
    pval <- stat <- dof <- numeric( p )  
    moda <- list()
    k <- 1   ## counter
    n <- length(target)  ## sample size
    con <- log(n)
    tool <- numeric( min(n, p) )
    ## percentages
    if ( min( target ) > 0  &  max( target ) < 1 )  target <- log( target / (1 - target) ) 
    ci_test <- "testIndSpeedglm"
  
    runtime <- proc.time()
    
    if (ncores <= 1) {
      
        if ( is.matrix(dataset)  &  is.null(wei) ) {
          mat <- Rfast::univglms( target, dataset, oiko = "normal", logged = TRUE ) 
          mat <- cbind(1:p, mat[, 2], mat[, 1])
          
        } else if ( is.data.frame(dataset)  &  is.null(wei) ) {
          #mod <- Rfast::regression(dataset, target)
          #pval <- pf(mod[1, ], mod[2, ], n - mod[2, ] - 1, lower.tail = FALSE, log.p = TRUE)
          #mat <- cbind(1:p, pval, mod[1, ])
		  stat <- numeric(p)
		  pval <- numeric(p)  
	      for (i in 1:p) {
            mod <- anova( lm(y~dataset[, i]) )
		    stat[i] <- mod[1, 5]
            pval[i] <- pf(stat[i], mod[1, 1], mod[2, 1], lower.tail = FALSE, log.p = TRUE)  		   
          } 
          mat <- cbind(1:p, pval, stat)	 
          
        } else {
          for (i in 1:p) {
            ww = speedglm::speedlm( target ~ dataset[, i], data = as.data.frame(dataset), weights = wei )
            suma = summary(ww)[[ 13 ]]
            stat[i] = suma[1]
            df1 = suma[2]
			      df2 <- suma[3]  
            pval[i] = pf(stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE)
          }
          mat <- cbind(1:p, pval, stat)
        }
      
    } else {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      mat <- matrix(0, p, 2)
      mod <- foreach( i = 1:p, .combine = rbind, .export = "speedlm", .packages = "speedglm" ) %dopar% {
        ww = speedlm( target ~ dataset[, i], data = as.data.frame(dataset), weights = wei ) 
        suma = summary(ww)[[ 13 ]]
        stat[i] = suma[1]
        df1 = suma[2]
		    df2 <- suma[3] 
        pval[i] = pf(stat, df1, df2, lower.tail = FALSE, log.p = TRUE)
        mat[i, ] = c(pval[i], stat[i]) 
      }
      stopCluster(cl)
      mat <- cbind(1:p, mod)      
    }
    
    colnames(mat) <- c( "variables", "log.p-value", "stat" )
    rownames(mat) <- 1:p
    sel <- which.min(mat[, 2])
    info <- matrix( numeric(3), ncol = 3 )
    sela <- sel
    
    if ( mat[sel, 2] < threshold ) {
      info[1, ] <- mat[sel, ]
      mat <- mat[-sel, , drop = FALSE] 
      mi = speedglm::speedlm( target ~ dataset[, sel], data = data.frame(dataset), weights = wei )
      if ( stopping == "adjrsq" ) {
        tool[1] <- as.numeric( summary( mi )[[ 11 ]] )		  
      } else if ( stopping  == "BIC" )   tool[1] <-  BIC(mi)
    
      moda[[ 1 ]] <- mi
      
    }  else  {
      info <- info  
      sela <- NULL
    }    
    ######
    ###### k equal to 2
    ######
    if ( info[k, 2] < threshold  &  nrow(mat) > 0 ) {
      
      k <- k + 1
      pn <- p - k + 1   
      
      if ( ncores <= 1 ) {
        fit1rss <- mi$RSS 
        d1 <- length( coef(mi) )
        for (i in 1:pn) {
          fit2 = speedglm::speedlm( target ~ dataset[, sel] + dataset[, mat[i, 1] ], data = as.data.frame(dataset), weights = wei )
          d2 = length( coef(fit2) )
          df1 = d2 - d1
          df2 = n - d2
          mat[i, 3] = ( (fit1rss - fit2$RSS)/df1 ) / ( fit2$RSS /df2 )
          mat[i, 2] = pf( mat[i, 3], df1, df2, lower.tail = FALSE, log.p = TRUE )
        }  		  
        
      } else {
        fit1rss <- mi$RSS 
        d1 <- length( coef(mi) ) 
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        stat <- pval <- numeric(pn) 
        mata <- matrix(0, pn, 2)
        mod <- foreach( i = 1:p, .combine = rbind, .export = "speedlm", .packages = "speedglm" ) %dopar% {
          fit2 <- speedlm( target ~ dataset[, sel] + dataset[, mat[i, 1] ], data = as.data.frame(dataset), weights = wei )
          d2 = length( coef(fit2) )
          df1 = d2 - d1
          df2 = n - d2
          stat[i] = (fit1rss - fit2$RSS)/df1 / ( fit2$RSS /df2 )
          pval[i] = pf(stat, df1, df2, lower.tail = FALSE, log.p = TRUE) 
          mata[i, ] <- c(pval[i], stat[i]) 
        }
        stopCluster(cl) 
        mat <- cbind(mat[, 1], mod)   
        
      }
      
      ina <- which.min(mat[, 2])
      sel <- mat[ina, 1]     
      ma = speedglm::speedlm( target ~ dataset[, sela] + dataset[, sel], data = as.data.frame(dataset), weights = wei )
            
      if ( stopping == "BIC" ) {
        
        if ( mat[ina, 2] < threshold ) {
          ma = speedglm::speedlm( target ~., data = as.data.frame( dataset[, c(sela, sel)] ), weights = wei )
          tool[k] <-  BIC(ma)	  
          
          if ( tool[ k - 1] - tool[ k ] <= tol ) {
            info <- info
            
          } else { 
            info <- rbind( info, mat[ina, ] )
            sela <- info[, 1]
            mat <- mat[-ina , ] 
            if ( !is.matrix(mat) )  mat <- matrix(mat, ncol = 3) 
            moda[[ k ]] <- ma
          }
          
        } else   info <- info
        
      } else if ( stopping == "adjrsq" ) {
        if ( mat[ina, 2] < threshold ) {
           ma = speedglm::speedlm( target ~., data = as.data.frame( dataset[, c(sela, sel)] ), weights = wei )
           tool[k] <- as.numeric( summary(ma)[[ 11 ]] )
           if ( tool[ k ] - tool[ k - 1 ] <= tol ) {
            info <- info
            
          } else { 
            info <- rbind( info, mat[ina, ] )
            sela <- info[, 1]
            mat <- mat[-ina , , drop = FALSE]
            moda[[ k ]] <- ma
          }
          
        } else  info <- info
        
      }
    }
    ######
    ###### k greater than 2
    ######
    if ( nrow(info) > 1  &  nrow(mat) > 0 ) {
      while ( info[k, 2] < threshold &  k < n - 15 & abs( tool[ k ] - tool[ k - 1 ] ) > tol &  nrow(mat) > 0 )  {
        
        k <- k + 1   
        pn <- p - k + 1 
        
        if ( ncores <= 1 ) {
           fit1rss <- ma$RSS 
           d1 <- length( coef(ma) )
           for (i in 1:pn) {
              fit2 = speedglm::speedlm( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), weights= wei )
              d2 = length( coef(fit2) )
              df1 = d2 - d1
              df2 = n - d2
              mat[i, 3] = (fit1rss - fit2$RSS)/df1 / ( fit2$RSS /df2 )
              mat[i, 2] = pf( mat[i, 3], df1, df2, lower.tail = FALSE, log.p = TRUE )
           }  

        } else {
          fit1rss <- ma$RSS 
          d1 <- length( coef(ma) ) 
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          stat <- pval <- numeric(pn) 
          mata <- matrix(0, pn, 2)
          mod <- foreach( i = 1:p, .combine = rbind, .export = "speedlm", .packages = "speedglm" ) %dopar% {
            fit2 <- speedlm( target ~., data = as.data.frame( dataset[, c(sela, mat[ i, 1]) ] ), weights = wei )
            d2 = length( coef(fit2) )
            df1 = d2 - d1
            df2 = n - d2
            stat[i] = (fit1rss - fit2$RSS)/df1 / ( fit2$RSS /df2 )
            pval[i] = pf(stat, df1, df2, lower.tail = FALSE, log.p = TRUE) 
            mata[i, ] <- c(pval[i], stat[i]) 
          }
          stopCluster(cl) 
          mat <- cbind( mat[, 1], mod )   
        }
        
        ina <- which.min(mat[, 2])
        sel <- mat[ina, 1]   
        ma = speedglm::speedlm( target ~., data = as.data.frame( dataset[, c(sela, sel)] ), weights = wei )
        
        if ( stopping == "BIC" ) {
          
          if ( mat[ina, 2] < threshold ) {
            ma = speedglm::speedlm( target ~., data = as.data.frame( dataset[, c(sela, sel)] ), weights = wei )
            tool[k] <-  BIC(ma)	  
            
            if ( tool[ k - 1] - tool[ k ] <= tol ) {
              info <- rbind(info, c( Inf, 0, 0 ) )
            } else { 
              info <- rbind( info, mat[ina, ] )
              sela <- info[, 1]
              mat <- mat[-ina , , drop = FALSE] 
              moda[[ k ]] <- ma
            }
            
          } else   info <- rbind(info, c( Inf, 0, 0 ) )
          
        } else if ( stopping == "adjrsq" ) {
          if ( mat[ina, 2] < threshold ) {
            ma = speedglm::speedlm( target ~., data = as.data.frame( dataset[, c(sela, sel)] ), weights = wei )
            tool[k] <- as.numeric( summary(ma)[[ 11 ]] )
            if ( tool[ k ] - tool[ k - 1 ] <= tol ) {
              info <- rbind(info, c( Inf, 0, 0 ) )
              
            } else { 
              info <- rbind( info, mat[ina, ] )
              sela <- info[, 1]
              mat <- mat[-ina , , drop = FALSE]
              moda[[ k ]] <- ma
            }
            
          } else  info <- rbind(info, c( Inf, 0, 0 ) )
          
        }

      }
      
    }
    
    runtime <- proc.time() - runtime
    
    d <- length(sela)
    final <- NULL

    if ( d >= 1 ) {
      final <- speedglm::speedlm( target ~., data = as.data.frame( dataset[, sela] ), weights = wei )
      info <- info[1:d, , drop = FALSE]
      info <- cbind( info, tool[ 1:d ] ) 
      colnames(info) <- c( "variables", "log.p-value", "stat", stopping )
      rownames(info) <- info[, 1]
    }
    
    result <- list(runtime = runtime, mat = t(mat), info = info, ci_test = ci_test, final = final ) 
    
  }  
  
  result
}
