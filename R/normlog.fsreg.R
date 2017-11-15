normlog.fsreg <- function(target, dataset, ini = NULL, threshold = 0.05, wei = NULL, tol = 2, heavy = FALSE, robust = FALSE, ncores = 1) {

  if ( !is.null(ini) ) {
    result <- normlog.fsreg_2(target, dataset, iniset = ini, wei = wei, threshold = threshold, tol = tol, heavy = heavy, robust = FALSE, ncores = ncores) 
    
  } else {  ## else do the classical forward regression
    
    threshold <- log(threshold)  ## log of the significance level
    p <- dim(dataset)[2]  ## number of variables
    devi <- dof <- numeric( p )  
    moda <- list()
    k <- 1   ## counter
    n <- length(target)  ## sample size
    con <- log(n)
    tool <- numeric( min(n, p) )
    dataset <- as.data.frame(dataset)
    
    runtime <- proc.time()
    if (heavy) {
      ini <- 2 * logLik( speedglm::speedglm( target ~ 1, data = dataset, family = gaussian(link = log), weights = wei ) )
    } else  ini <- 2 * logLik( glm( target ~ 1, family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE ) )
    
         if (ncores <= 1) {

           if ( !heavy ) {  		
              for (i in 1:p) {
                mi <- glm( target ~ dataset[, i], family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE )
                devi[i] <- 2 * as.numeric( logLik(mi) )
                dof[i] <- length( mi$coefficients ) 
              }
              
            } else {
              for (i in 1:p) {
                mi <- speedglm::speedglm( target ~ dataset[, i], data = dataset, family = gaussian(link = log), weights = wei )
                devi[i] <- 2 * as.numeric( logLik(mi) )
                dof[i] <- length( mi$coefficients ) 
              }
            }  
            stat <- devi - ini
            pval <- pchisq( stat, dof - 1, lower.tail = FALSE, log.p = TRUE )
        } else {
          #if ( robust == FALSE ) {  ## Non robust
          if ( !heavy )	 {
            cl <- makePSOCKcluster(ncores)
            registerDoParallel(cl)
            mata <- matrix(0, p, 2)
            mod <- foreach( i = 1:p, .combine = rbind) %dopar% {
              ww <- glm( target ~ dataset[, i], family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE )
              mata[i, ] <- c( 2 * as.numeric( logLik(ww) ), length( coef( ww ) )  )
            }
            stopCluster(cl)
            
          } else {
            cl <- makePSOCKcluster(ncores)
            registerDoParallel(cl)
            mod <- foreach( i = 1:p, .combine = rbind, .export = "speedglm", .packages = "speedglm") %dopar% {
              ww <- speedglm( target ~ dataset[, i], data = dataset, family = gaussian(link = log), weights = wei )
              return( c( 2 * as.numeric( logLik(ww) ), length( coef( ww ) ) ) )
            }
            
            stopCluster(cl)	  	  
            stat = mod[, 1] - ini
            pval = pchisq( stat, mod[, 2] - 1, lower.tail = FALSE, log.p = TRUE )
          }
        }  
        mat <- cbind(1:p, pval, stat) 
        colnames(mat) <- c( "variables", "log.p-value", "stat" )
        rownames(mat) <- 1:p
        sel <- which.min(pval)
        info <- matrix( numeric(3), ncol = 3 )
        sela <- sel
        
        if ( mat[sel, 2] < threshold ) {
          info[1, ] <- mat[sel, ]
          mat <- mat[-sel, , drop = FALSE] 
          if ( !heavy )  {
            mi <- glm( target ~ dataset[, sel], family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE )
            tool[1] <- BIC( mi )
          } else {
            mi <- speedglm::speedglm( target ~ dataset[, sel], data = dataset, family = gaussian(link = log), weights = wei )
            tool[1] <-  - 2 * as.numeric( logLik(mi) ) + length( coef(mi) ) * con 		
          }
          moda[[ 1 ]] <- mi
        }  else  {
          info <- info  
          sela <- NULL
        }
        ##########
        #####   k equals 2
        ########## 
        if ( info[k, 2] < threshold  &  nrow(mat) > 0 ) {
          
          k <- 2
          pn <- p - k + 1   
          ini <- 2 * logLik(mi)
          do <- length( mi$coefficients ) 
          devi <- dof <- numeric( pn )  
          
          if ( ncores <= 1 ) {
            devi <- dof <- numeric(pn)
            if ( !heavy ) {
              for ( i in 1:pn ) {
                ww <- glm( target ~., data = dataset[, c(sela, mat[i, 1]) ], family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE )
                devi[i] <- 2 * as.numeric( logLik(ww) )
                dof[i] <- length( ww$coefficients )          
              }
              
            } else {
              for ( i in 1:pn ) {
                ww <- speedglm::speedglm( target ~., data = dataset[, c(sela, mat[i, 1]) ], family = gaussian(link = log), weights = wei )
                devi[i] <- 2 * as.numeric( logLik(ww) )
                dof[i] <- length( ww$coefficients )          
              }		  
            }
            stat <- devi - ini
            pval <- pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
            
          } else {
            
            if ( !heavy ) {
              cl <- makePSOCKcluster(ncores)
              registerDoParallel(cl)
              mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
                ww <- glm( target ~., data = dataset[, c(sela, mat[i, 1]) ], family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE )
                return( c( 2 * as.numeric( logLik(ww) ), length( coef( ww ) ) ) )
              }
              stopCluster(cl)
              
            }	else {
              cl <- makePSOCKcluster(ncores)
              registerDoParallel(cl)
              mod <- foreach( i = 1:pn, .combine = rbind, .export = "speedglm", .packages = "speedglm") %dopar% {
                ww <- speedglm( target ~., data = dataset[, c(sela, mat[i, 1]) ], family = gaussian(link = log), weights = wei )
                return( c( 2 * as.numeric( logLik(ww) ), length( coef( ww ) ) ) )
              }
              stopCluster(cl)		  
            }
            stat <- mod[, 1] - ini
            pval <- pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
          }
          
          mat[, 2:3] <- cbind(pval, stat)
          ina <- which.min(mat[, 2])
          sel <- mat[ina, 1]    
          
          if ( mat[ina, 2] < threshold ) {
            if ( !heavy ) {
              ma <- glm( target ~ dataset[, sela] + dataset[, sel], family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE )
              tool[k] <- BIC( ma )
            } else {
              ma <- speedglm::speedglm( target ~ dataset[, sela] + dataset[, sel], data = as.data.frame(dataset), family = gaussian(link = log), weights = wei )
              tool[k] <-  - 2 * as.numeric( logLik(ma) ) + length( coef(ma) ) * con		  
            }	

            if ( tool[ k - 1 ] - tool[ k ] <= tol ) {
              info <- rbind(info, c( 1e300, 0, 0 ) )
              
            } else { 
              info <- rbind(info, c( mat[ina, ] ) )
              sela <- info[, 1]
              mat <- mat[-ina, , drop = FALSE] 
              moda[[ k ]] <- ma
            }
            
          } else   info <- rbind(info, c( 1e300, 0, 0 ) )
        }
        #######
        ####   k greater than 2
        ####### 
        if ( info[k, 2] < threshold  &  nrow(mat) > 0 ) {
          
          while ( ( info[k, 2] < threshold ) &  ( k < n ) & ( tool[ k - 1 ] - tool[ k ] > tol ) & ( nrow(mat) > 0 ) )  {
            
            ini = 2 * as.numeric( logLik(moda[[ k ]]) )
            do = length( coef( moda[[ k ]]  ) ) 
            k <- k + 1   
            pn <- p - k  + 1
            devi = dof = numeric(pn) 
            
            if (ncores <= 1) {  
              #if ( robust == FALSE ) {  ## Non robust
              if ( !heavy ) {
                for ( i in 1:pn ) {
                  ma <- glm( target ~., data = dataset[, c(sela, mat[i, 1] ) ], family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE )
                  devi[i] <- 2 * as.numeric( logLik(ma) )
                  dof[i] = length( ma$coefficients ) 
                }
              } else {
                for ( i in 1:pn ) {
                  ma <- speedglm( target ~., data = dataset[, c(sela, mat[i, 1] ) ], family = gaussian(link = log), weights = wei )
                  devi[i] <- 2 * as.numeric( logLik(ma) )
                  dof[i] = length( coef( ma ) ) 
                }
                
              }
              stat = devi - ini
              pval = pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
              
            } else {
              #if ( robust == FALSE ) {  ## Non robust
              if ( !heavy ) {
                cl <- makePSOCKcluster(ncores)
                registerDoParallel(cl)
                mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
                  ww <- glm( target ~., data = dataset[, c(sela, mat[i, 1]) ], family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE )
                  return( c( 2 * as.numeric( logLik(ww) ), length( coef( ww ) ) ) )
                }
                stopCluster(cl)
                
              } else {
                cl <- makePSOCKcluster(ncores)
                registerDoParallel(cl)
                mod <- foreach( i = 1:pn, .combine = rbind, .export = "speedglm", .packages = "speedglm") %dopar% {
                  ww <- speedglm( target ~., data = dataset[, c(sela, mat[i, 1]) ], family = gaussian(link = log), weights = wei )
                  return( c( 2 * as.numeric( logLik(ww) ), length( coef( ww ) ) ) )
                }
                stopCluster(cl)
              }
              stat <- mod[, 1] - ini
              pval <- pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
            }
            
            mat[, 2:3] <- cbind(pval, stat)
            ina <- which.min(mat[, 2])
            sel <- mat[ina, 1]    
            
            if ( mat[ina, 2] < threshold ) {
              #if ( robust == FALSE ) {
              if ( !heavy ) {
                ma <- glm( target ~., data = dataset[, c(sela, sel) ], family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE )
                tool[k] <- BIC( ma )
              } else {
                ma <- speedglm::speedglm( target ~., data = dataset[, c(sela, sel) ], family = gaussian(link = log), weights = wei )
                tool[k] <-  - 2 * as.numeric( logLik(ma) ) + length( coef(ma) ) * con			  
              }
              if ( tool[ k - 1 ] - tool[ k  ] <= tol ) {
                info <- rbind(info, c( 1e300, 0, 0 ) )
                
              } else { 
                info <- rbind( info, mat[ina, ] )
                sela <- info[, 1]
                mat <- mat[-ina, , drop = FALSE]
                moda[[ k ]] <- ma
              } 
              
            } else  info <- rbind(info, c( 1e300, 0, 0 ) )
          }
        } 
        
        runtime <- proc.time() - runtime
        final <- NULL
    
      d <- p - dim(mat)[1]
    
      if ( d >= 1 ) {
         if ( !heavy ) {
           final <- glm( target ~., data = dataset[, sela, drop = FALSE], family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE )
         } else   final <- speedglm::speedglm( target ~., data = dataset[, sela, drop = FALSE], family = gaussian(link = log), weights = wei )
          info <- info[1:d, , drop = FALSE]
          info <- cbind( info, tool[ 1:d ] ) 
          colnames(info) <- c( "variables", "log.p-values", "stat", "BIC" )
          rownames(info) <- info[, 1]
          mat <- cbind( mat, exp(mat[, 2]) )
          colnames(mat)[4] <- "p-value"
      }
      result = list( runtime = runtime, mat = mat, info = info, ci_test = "testIndNormLog", final = final ) 

  }
  
  result
}    










