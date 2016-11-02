lm.fsreg <- function(target, dataset,ini = NULL, threshold = 0.05, wei = NULL, stopping = "BIC", tol = 2, heavy = FALSE, robust = FALSE, ncores = 1 ) {
  
  
  ###### If there is an initia set of variables do this function
  
  if ( !is.null(ini) ) {
      result <- lm.fsreg_2(target, dataset, iniset = ini, threshold = threshold, wei = NULL, stopping = stopping, tol = tol, heavy = heavy, robust = robust, ncores = ncores) 
      
  } else {  ## else do the classical forward regression
    
  
  threshold = log(threshold)
  
  p <- dim(dataset)[2]  ## number of variables
  pval <- stat <- dof <- numeric( p )  
  moda <- list()
  k <- 1   ## counter
  n <- length(target)  ## sample size
  con <- log(n)
  tool <- numeric( min(n, p) )
  
  ci_test <- "testIndReg"
  
  if ( heavy == TRUE )   robust = FALSE    ;   ci_test <- "testIndSpeedglm"
   
  # cont = robust::lmRob.control(mxr = 2000, mxf = 2000, mxs = 2000 )

  runtime <- proc.time()
    
    if (ncores <= 1) {

     if ( robust == FALSE ) {  ## Non robust
       
       if ( is.matrix(dataset)  &  is.null(wei) ) {
         mat <- Rfast::univglms( target, dataset, oiko = "normal", logged = TRUE ) 
         mat <- cbind(1:p, mat[, 2], mat[, 1])
         
       } else if ( is.data.frame(dataset)  &  is.null(wei) ) {
         mod <- Rfast::regression(dataset, target)
         pval <- pf(mod[1, ], mod[2, ], n - mod[2, ] - 1, lower.tail = FALSE, log.p = FALSE)
         mat <- cbind(1:p, mod[1, ], pval)
       
       } else {

         if ( heavy == FALSE ) {
           for (i in 1:p) {
             ww = lm( target ~ dataset[, i], weights = wei )
             tab = anova(ww)
             stat[i] = tab[1, 4] 
             df1 = tab[1, 1]    ;   df2 = tab[2, 1]
             pval[i] = pf( stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE )
           }
		   
         } else {
		   for (i in 1:p) {
             ww = speedglm:speedlm( target ~ dataset[, i], data = as.data.frame(dataset), weights = wei )
             suma = summary(ww)[[ 13 ]]
             stat = suma[1]
             dof = suma[3]
             pvalue = pf(stat, 1, dof, lower.tail = FALSE, log.p = TRUE)
           }
		   
		 } 
           mat <- cbind(1:p, pval, stat)
		 
       }
      
      
     } else {  ## Robust
       for (i in 1:p) {
        # ww = robust::lmRob( target ~ dataset[, i], control = cont )
        # tab = aovlmrob( ww )
        ww = MASS::rlm(target ~ dataset[, i ], maxit = 2000) 
        stat[i] = 2 * as.numeric( logLik(ww) )
        dof[i] = length( coef(ww) )
       }

       fit0 = MASS::rlm(target ~ 1, maxit = 2000, weights = wei) 
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
            ww <- lm( target ~ dataset[, i], weights = wei )
            tab <- anova( ww )
            stat[i] <- tab[1, 4] 
            df1 <- tab[1, 1]   ;  df2 = tab[2, 1]
            pval[i] <- pf( stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE )
            mat[i, ] <- c(pval[i], stat[i]) 
         }
          stopCluster(cl)
		 
         } else {
		  cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          mat <- matrix(0, p, 2)
          mod <- foreach( i = 1:p, .combine = rbind, .export = "speedlm", .packages = "speedglm" ) %dopar% {
            ww = speedlm( target ~ dataset[, i], data = as.data.frame(dataset), weights = wei ) 
			      suma = summary(ww)[[ 13 ]]
            stat[i] = suma[1]
            dof = suma[3]
            pval[i] = pf(stat, 1, dof, lower.tail = FALSE, log.p = TRUE)
            mat[i, ] = c(pval[i], stat[i]) 
          }
          stopCluster(cl)
		}  
		
      } else {  ## Robust
        fit0 = MASS::rlm(target ~ 1, maxit = 2000, weights = wei) 
        stat0 = 2 * logLik(fit0)
        
         cl <- makePSOCKcluster(ncores)
         registerDoParallel(cl)
         mat <- matrix(0, p, 2)
         mod <- foreach( i = 1:p, .combine = rbind, .export = "rlm", .packages = "MASS" ) %dopar% {
           # ww <- lmRob( target ~ dataset[, i], control = cont )
           # tab <- aovlmrob( ww )
           # mat[i, ] <- c(tab[2, 3], tab[2, 2] ) 
           ww = rlm( target ~ dataset[, i], maxit = 2000, weights = wei ) 
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
    info <- matrix( numeric(3), ncol = 3 )
    sela <- sel
    
    if ( mat[sel, 2] < threshold ) {
      info[1, ] <- mat[sel, ]
      mat <- mat[-sel, ] 
      if ( !is.matrix(mat) )   mat <- matrix(mat, ncol = 3) 
      mat <- mat[ order( mat[, 2] ), ] 

      if ( stopping == "adjrsq" ) {

        if ( robust == FALSE ) {
		      if ( heavy == FALSE ) {
            mi = lm( target ~ dataset[, sel], weights = wei )
            tool[1] <- as.numeric( summary( mi )[[ 9 ]] )
		      } else {
            mi = speedglm::speedlm( target ~ dataset[, sel], data = as.data.frame(dataset), weights = wei )
            tool[1] <- as.numeric( summary( mi )[[ 11 ]] )		  
		      }
		   
        } else {
          mi = MASS::rlm( target ~ dataset[, sel], maxit = 2000, weights = wei )
          r2 = cor( target, fitted(mi) )^2
          tool[1] = 1 - (1 - r2) * (n - 1) / ( n - length( coef(mi) ) - 1 )
        }         
      } else if ( stopping  == "BIC" ) { 
        if ( robust == FALSE ) {
		      if ( heavy == FALSE ) {
            mi = lm( target ~ dataset[, sel], weights = wei )
            tool[1] <- BIC( mi )
		      } else {
            mi = speedglm::speedlm( target ~ dataset[, sel], data = as.data.frame(dataset), weights = wei )
            tool[1] <-  - 2 * mi$logLik + length( coef(mi) ) * con
		      }	 
          
        } else {
          mi = MASS::rlm( target ~ dataset[, sel], maxit = 2000, weights = wei )
          tool[1] <-  BIC(mi)
        }
      }
      moda[[ 1 ]] <- mi
      
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
             ww = lm( target ~ dataset[, sel] + dataset[, mat[i, 1] ], weights = wei )
             tab = anova( ww )
             mat[i, 3] = tab[2, 4] 
             df1 = tab[2, 1]   ;  df2 = tab[3, 1]
             mat[i, 2] = pf( mat[i, 3], df1, df2, lower.tail = FALSE, log.p = TRUE )
           }
		   
          } else {
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
		  
		      }
		  
        } else {
          do = 2 * as.numeric( logLik( moda[[ 1 ]] ) )
          fr = length( coef( moda[[ 1 ]] ) )
          sta = dof = numeric(pn)
          
          for (i in 1:pn) {
            ww = MASS::rlm( target ~ dataset[, sel] + dataset[, mat[i, 1] ], maxit = 2000, weights = wei )
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
             ww <- lm( target ~ dataset[, sel] + dataset[, mat[i, 1] ], weights = wei )
             tab <- anova( ww )
             stat[i] <- tab[2, 4] 
             df1 <- tab[2, 1]   ;  df2 = tab[3, 1]
             pval[i] <- pf( stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE )
             mata[i, ] <- c(pval[i], stat[i]) 
           }
           stopCluster(cl)
		   
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
		  } 

        } else {  ## Robust
          do = 2 * as.numeric( logLik( moda[[ 1 ]] ) )
          fr = length( coef( moda[[ 1 ]] ) )
          
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          mata <- matrix(0, pn, 2)
          mod <- foreach( i = 1:pn, .combine = rbind, .export = c("rlm"), .packages = "MASS" ) %dopar% {
            ww <- rlm( target ~ dataset[, sel] + dataset[, mat[i, 1] ], maxit = 2000, weights = wei )
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
                ma = lm( target ~ dataset[, sela] + dataset[, sel], weights = wei )
                tool[k] <- as.numeric( summary( ma )[[ 9 ]] )
			        } else {
			          ma = speedglm::speedlm( target ~ dataset[, sela] + dataset[, sel], data = as.data.frame(dataset), weights = wei )
                tool[k] <- as.numeric( summary( ma )[[ 11 ]] )
			        }
			  
            } else {         
              ma = MASS::rlm( target ~ dataset[, sela] + dataset[, sel], maxit = 2000, weights = wei )
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

          } else {
           info <- rbind(info, c( 1e300, 0, 0 ) )
          }

        } else if ( stopping == "BIC" ) {
		
          if ( mat[ina, 2] < threshold ) {
		  
            if ( robust == FALSE ) { 
              if ( heavy == FALSE ) {			
                ma = lm( target ~ dataset[, sela] + dataset[, sel], weights = wei )
                tool[2] <- BIC( ma )
		      } else {
                ma = speedglm::speedlm( target ~ dataset[, sela] + dataset[, sel], data = as.data.frame(dataset), weights = wei )
                tool[2] <-  - 2 * logLik(ma) + length( coef(ma) ) * con			  
			  }
			  
            } else {
              ma = MASS::rlm( target ~ dataset[, sela] + dataset[, sel], maxit = 2000, weights = wei )
              tool[2] = BIC(ma)
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


 
 
     ######
     ###### k greater than 2
     ######

	 
    if ( nrow(info) > 1  &  nrow(mat) > 0 ) {
      while ( ( info[k, 2] < threshold ) &  ( k < n ) & ( abs( tool[ k ] - tool[ k - 1 ] ) > tol ) & ( nrow(mat) > 0 ) )  {
         
        k <- k + 1   
        pn <- p - k + 1 
        
        if ( ncores <= 1 ) {
		
          if ( robust == FALSE ) {
		        if ( heavy == FALSE ) {
              for ( i in 1:pn ) {
                 ww <- lm( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), weights = wei )
                 tab <- anova( ww )
                 mat[i, 3] <- tab[k, 4] 
                 df1 <- tab[k, 1]   ;  df2 = tab[k + 1, 1]
                 mat[i, 2] <- pf( mat[i, 3], df1, df2, lower.tail = FALSE, log.p = TRUE )
              }
			 
			      } else {	
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
           } 		  
			
          } else {
            
            do = 2 * as.numeric( logLik( moda[[ k - 1 ]] ) )
            fr = length( coef( moda[[ k - 1 ]] ) )
            sta = dof = numeric(pn)
            
              for (i in 1:pn) {
                ww = MASS::rlm( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), maxit = 2000, weights = wei )
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
                ww <- lm( target ~., data = as.data.frame( dataset[, c(sela, mat[ i, 1]) ] ), weights = wei )
                tab <- anova(ww)
                stat[i] <- tab[k, 4] 
                df1 <- tab[k, 1]   ;  df2 = tab[k + 1, 1]
                pval[i] <- pf( stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE )
                mata[i, ] <- c( pval[i], stat[i] ) 
              }
              stopCluster(cl)
			
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
		   
		       }
		   
          } else {  ## Robust
            do = 2 * as.numeric( logLik( moda[[ k - 1 ]] ) )
            fr = length( coef( moda[[ k - 1 ]] ) )
            
            cl <- makePSOCKcluster(ncores)
            registerDoParallel(cl)
            mata <- matrix(0, pn, 2)
            mod <- foreach( i = 1:pn, .combine = rbind, .export = c("rlm"), .packages = "MASS" ) %dopar% {
              ww <- MASS::rlm( target ~., data = data.frame( dataset[, c(sela, mat[ i, 1]) ] ), maxit = 2000, weights = wei )
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
			    if ( heavy == FALSE ) {
                  ma = lm( target ~., data = as.data.frame( dataset[, c(sela, sel)] ), weights = wei )
                  tool[k] <- BIC( ma )
				} else {
                  ma = speedglm::speedlm( target ~., data = as.data.frame( dataset[, c(sela, sel)] ), weights = wei )
                  tool[k] <-  - 2 * ma$logLik + length( coef(ma) ) * con			  
				}
				
              } else {
                # ma = robust::lmRob( target ~., data = as.data.frame( dataset[, c(sela, sel)] ), control = cont )
                ma = MASS::rlm( target ~., data = as.data.frame( dataset[, c(sela, sel)] ), maxit = 2000, weights = wei )
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
                  ma = lm( target ~., data = as.data.frame( dataset[, c(sela, sel)] ), weights = wei )
                  tool[k] <- as.numeric( summary(ma)[[ 9 ]] )
				} else {
				  ma = speedglm::speedlm( target ~., data = as.data.frame( dataset[, c(sela, sel)] ), weights = wei )
                  tool[k] <- as.numeric( summary(ma)[[ 11 ]] )
				}
				
              } else {
                ma = MASS::rlm( target ~., data = as.data.frame( dataset[, c(sela, sel)] ), maxit = 2000, weights = wei )
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

            } else {
              info <- rbind(info, c( 1e300, 0, 0 ) )
            }
        
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
 
        xx <- as.data.frame( dataset[, sela] )
        colnames(xx) <- paste("V", sela, sep = "") 

          if ( robust == FALSE ) {
		    if ( heavy == FALSE ) {
              models[[ 1 ]] <- final <- lm( target ~., data = as.data.frame( xx ), weights = wei )
			} else { 
              models[[ 1 ]] <- final <- speedglm::speedlm( target ~., data = as.data.frame( xx ), weights = wei )
            } 
		  } else {
            models[[ 1 ]] <- final <- MASS::rlm( target ~., data = as.data.frame( xx ), maxit = 2000, weights = wei )
          }
        
      } else {
        for (i in 1:d) {
		
          if ( robust == FALSE ) {
		    if ( heavy == FALSE ) {
              models[[ i ]] <- lm( target ~., data = as.data.frame( xx[, 1:i] ), weights = wei )
			} else {  
              models[[ i ]] <- speedglm::speedlm( target ~., data = as.data.frame( xx[, 1:i] ), weights = wei )
            }
			
			} else { 
            models[[ i ]] <- MASS::rlm( target ~., data = data.frame( xx[, 1:i]), maxit = 2000, weights = wei )   
          }
		  
        }
      }

      final <- models[[ d ]]

      info <- info[1:d, ]
      if ( d == 1 )  info <- matrix(info, nrow = 1)
      info <- cbind( info, tool[ 1:d ] ) 
      colnames(info) <- c( "variables", "p-value", "stat", stopping )
      rownames(info) <- info[, 1]
    }

    result <- list(mat = t(mat), info = info, models = models, final = final, ci_test = ci_test, runtime = runtime ) 
    
  }  
  
  result
    
}
  