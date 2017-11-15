bic.glm.fsreg <- function( target, dataset, wei = NULL, tol = 0, heavy = FALSE, robust = FALSE, ncores = 1) {

  p <- dim(dataset)[2]  ## number of variables
  bico <- numeric( p )
  moda <- list()
  k <- 1   ## counter
  n <- length(target)  ## sample size
  tool <- NULL
  oiko <- NULL
  info <- matrix( 0, ncol = 2 )
  # cont = robust::lmRob.control(mxr = 2000, mxf = 2000, mxs = 2000 )
  #check for NA values in the dataset and replace them with the variable median or the mode
  if ( any( is.na(dataset) ) ) {
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    if ( is.matrix(dataset) )  {
      dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
    } else {
       poia <- unique( which( is.na(dataset), arr.ind = TRUE )[, 2] )
      for( i in poia )  {
        xi <- dataset[, i]
        if ( is.numeric(xi) ) {                    
          xi[ which( is.na(xi) ) ] <- median(xi, na.rm = TRUE) 
        } else if ( is.factor( xi ) )     xi[ which( is.na(xi) ) ] <- levels(xi)[ which.max( as.vector( table(xi) ) )]
        dataset[, i] <- xi
      }
    }
  }
  
  dataset <- as.data.frame(dataset)
  if ( heavy )  con <- log(n)
  #########
  ## if it is binomial or poisson regression
  #########
  la <- length( unique(target) ) 

  if ( la == 2  | ( sum( round(target) - target ) == 0  &  la > 2 ) ) {

   if ( la == 2 ) {
      oiko <- binomial(link = logit)  ## binomial regression
      ci_test <- "testIndLogistic"
    } else {
      oiko <- poisson(link = log)  ## poisson regression
      ci_test <- "testIndPois"
    }
    
    durat <- proc.time()
    #if ( robust == FALSE ) {
	  if ( !heavy ) {
        ini = BIC( glm( target ~ 1, family = oiko, weights = wei, y = FALSE, model = FALSE ) )   ## initial BIC
	  } else {
	    inimod <- speedglm::speedglm( target ~ 1, data = dataset, family = oiko, weights = wei )
	    ini <-  - 2 * inimod$logLik + con   ## initial BIC
	    ci_test <- "testIndSpeedglm"
	  }	
    #} else {
    #  ini = robust::glmRob( target ~ 1, family = oiko, maxit = maxit )$deviance + log(n)  ## initial BIC
    #}
    	
    if (ncores <= 1) {
	
      #if ( robust == FALSE ) {  ## Non robust
	    if ( !heavy ) {
          for (i in 1:p) {
            mi <- glm( target ~ dataset[, i], family = oiko, weights = wei, y = FALSE, model = FALSE )
            bico[i] <- BIC( mi )
          }
		  } else {
        for (i in 1:p) { 
          mi <- speedglm::speedglm( target ~ dataset[, i], data = dataset[, i], family = oiko, weights = wei )
	        bico[i] <-  - 2 * mi$logLik + length( coef(mi) ) * con   ## initial BIC
		    }	
		  }  
      # } else {  ## Robust
      #   for (i in 1:p) {
      #     mi <- robust::glmRob( target ~ dataset[, i], family = oiko, maxit = maxit )
      #     bico[i] <- mi$deviance + length( coef( mi ) ) * log(n)
      #   }
      # }
        
      mat <- cbind(1:p, bico)
      if( any( is.na(mat) ) )   mat[ which( is.na(mat) ) ] = ini

    } else {
      #if ( robust == FALSE ) {  ## Non robust
	    if ( !heavy ) {
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          bico <- numeric(p)
          mod <- foreach( i = 1:p, .combine = rbind) %dopar% {
            ww <- glm( target ~ dataset[, i], family = oiko, weights = wei )
            bico[i] <- BIC( ww )
          }
          stopCluster(cl)
		  }  else {
		     cl <- makePSOCKcluster(ncores)
         registerDoParallel(cl)
         bico <- numeric(p)
         mod <- foreach( i = 1:p, .combine = rbind, .export = "speedglm", .packages = "speedglm") %dopar% {
            ww <- speedglm( target ~ dataset[, i], data = dataset[, i], family = oiko, weights = wei )
 	          bico[i] <-  -2 * ww$logLik + length( coef(ww) ) * con   ## initial BIC
         }
         stopCluster(cl)
		
		  }  

      # } else {  ## Robust
      #   cl <- makePSOCKcluster(ncores)
      #   registerDoParallel(cl)
      #   bico <- numeric(p)
      #   mod <- foreach( i = 1:p, .combine = rbind, .export = "glmRob", .packages = "robust" ) %dopar% {
      #     ww <- robust::glmRob( target ~ dataset[, i], family = oiko, maxit = maxit )
      #     bico[i] <- ww$deviance + length( coef( ww ) ) * log(n)
      #   }
      #   stopCluster(cl)
      # }

      mat <- cbind(1:p, mod)
      if ( any( is.na(mat) ) )    mat[ which( is.na(mat) ) ] = ini
      
    }

    colnames(mat) <- c("variable", "BIC")
    rownames(mat) <- 1:p
    sel <- which.min( mat[, 2] )
    sela <- sel

    if ( ini - mat[sel, 2] > tol ) {

      info[1, ] <- mat[sel, ]
      mat <- mat[-sel, , drop = FALSE]
      #if ( robust == FALSE ) {
	    if ( !heavy ) {
          mi <- glm( target ~ dataset[, sel], family = oiko, weights = wei, y = FALSE, model = FALSE )
          tool[1] <- BIC( mi )
		   } else {
          mi <- speedglm::speedglm( target ~ dataset[, sel], data = dataset[, sel], family = oiko, weights = wei )
          tool[1] <-  - 2 * mi$logLik + length( coef(mi) ) * con   ## initial BIC
		   }  
      #} else {
      #  mi = robust::glmRob( target ~ dataset[, sel], family = oiko, maxit = maxit )
      #  tool[1] <- mi$deviance + length( coef( mi ) ) * log(n)
      #  if ( is.na(tool[1]) )  tool[1] <- ini
      #}
      moda[[ 1 ]] <- mi
    }  else  {
      info <- info  
      sela <- NULL
    }
    ######
    ###     k equals 2
    ######
    if ( length(moda) > 0  &  nrow(mat) > 0 ) {

      k <- 2
      pn <- p - k  + 1
      mod <- list()

      if ( ncores <= 1 ) {
        bico <- numeric( pn )
        #if ( robust == FALSE ) {  ## Non robust
		    if ( !heavy ) {
          for ( i in 1:pn ) {
            ma <- glm( target ~., data =  dataset[, c(sel, mat[i, 1]) ], family = oiko, weights = wei, y = FALSE, model = FALSE )
            bico[i] <- BIC( ma )
          }
        } else {
          for ( i in 1:pn ) {
            ma <- speedglm::speedglm( target ~., data = dataset[, c(sel, mat[i, 1]) ], family = oiko, weights = wei )
            bico[i] <-  - 2 * ma$logLik + length( coef(ma) ) * con
          }		  
		    }
        # } else {  ## Robust
        #   for ( i in 1:pn ) {
        #     ma <- robust::glmRob( target ~., data = as.data.frame( dataset[, c(sel, mat[i, 1]) ] ), family = oiko, maxit = maxit )
        #     bico[i] <- ma$deviance + length( coef( ma ) ) * log(n)
        #   }
        # }
        mat[, 2] <- bico

      } else {

        #if ( robust == FALSE ) {  ## Non robust
		  if ( !heavy ) {
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          bico <- numeric(pn)
          mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
            ww <- glm( target ~ dataset[, sel ] + dataset[, mat[i, 1] ], family = oiko, weights = wei )
            bico[i] <- BIC( ww )
          }
            stopCluster(cl)
		  } else {
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          bico <- numeric(pn)
          mod <- foreach( i = 1:pn, .combine = rbind, export = "speedglm", .packages = "speedglm") %dopar% {
            ww <- speedglm( target ~ dataset[, sel ] + dataset[, mat[i, 1] ], data = dataset, family = oiko, weights = wei )
            bico[i] <-  - 2 * ww$logLik + length( coef(ww) ) * con
          }
            stopCluster(cl)		  
		  }	
        # } else {  ## Robust
        #   cl <- makePSOCKcluster(ncores)
        #   registerDoParallel(cl)
        #   bico <- numeric(pn)
        #   mod <- foreach( i = 1:pn, .combine = rbind, .export = "glmRob", .packages = "robust" ) %dopar% {
        #     ww <- glmRob( target ~ dataset[, sel ] + dataset[, mat[i, 1] ], family = oiko, maxit = maxit )
        #     bico[i] <- ww$deviance + length( coef( ww ) ) * log(n)
        #   }
        # 
        #   stopCluster(cl)
        # }
        mat[, 2] <- mod
      }

      ina <- which.min( mat[, 2] )
      sel <- mat[ina, 1]
      if ( tool[1] - mat[ina, 2] <= tol ) {
        info <- info
      } else {
        tool[2] <- mat[ina, 2]
        info <- rbind(info, mat[ina, ] )
        sela <- info[, 1]
        mat <- mat[-ina, , drop = FALSE]
     }

   }
   #########
   ####      k is greater than 2
   #########
    if ( nrow(info) > 1  &  nrow(mat) > 0 ) {
      while (  k < n - 15  &  tool[ k - 1 ] - tool[ k ] > tol  & nrow(mat) > 0 ) {

        k <- k + 1
        pn <- p - k + 1

        if (ncores <= 1) {
          #if ( robust == FALSE ) {  ## Non robust
		      if ( !heavy ) {
            for ( i in 1:pn ) {
              ma <- glm( target ~., data = dataset[, c(sela, mat[i, 1]) ], family = oiko, weights = wei, y = FALSE, model = FALSE )
              mat[i, 2] <- BIC( ma )
            }
          } else {
            for ( i in 1:pn ) {
              ma <- speedglm::speedglm( target ~., data = dataset[, c(sela, mat[i, 1]) ], family = oiko, weights = wei )
              mat[i, 2] <-  - 2 * ma$logLik + length( coef(ma) ) * con
            }		  
		      }   
          # } else {  ## Robust
          #   for ( i in 1:pn ) {
          #     ma <- robust::glmRob( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), family = oiko, maxit = maxit )
          #     mat[i, 2] <- ma$deviance + length( coef( ma ) ) * log(n)
          #   }
          # }

        } else {
          #if ( robust == FALSE ) {  ## Non robust
		      if ( !heavy ) {
            cl <- makePSOCKcluster(ncores)
            registerDoParallel(cl)
            bico <- numeric(pn)
            mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
              ww <- glm( target ~., data = dataset[, c(sela, mat[i, 1]) ], family = oiko, weights = wei )
              bico[i] <- BIC( ww )
            }
            stopCluster(cl)
          } else {
            cl <- makePSOCKcluster(ncores)
            registerDoParallel(cl)
            bico <- numeric(pn)
            mod <- foreach( i = 1:pn, .combine = rbind, .export = "speedglm", .packages = "speedglm") %dopar% {
              ww <- speedglm( target ~., data = dataset[, c(sela, mat[i, 1]) ], family = oiko, weights = wei )
              bico[i] <-  - 2 * ww$logLik + length( coef(ww) ) * con
            }
            stopCluster(cl)		  
		      }
		  
          # } else {  ## Robust
          #   cl <- makePSOCKcluster(ncores)
          #   registerDoParallel(cl)
          #   bico <- numeric(pn)
          #   mod <- foreach( i = 1:pn, .combine = rbind, .export = "glmRob", .packages = "robust" ) %dopar% {
          #     ww <- glmRob( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), family = oiko, maxit = maxit )
          #     bico[i] = ww$deviance + length( coef( ww ) ) * log(n)
          #   }
          # 
          #   stopCluster(cl)
          # }
          mat[, 2] <- mod

        }

        ina <- which.min( mat[, 2] )
        sel <- mat[ina, 1]

        if ( tool[k - 1] - mat[ina, 2]  <= tol ) {
          info <- rbind( info,  c( -10, 1e300 ) )
          tool[k] <- Inf

        } else {
          tool[k] <- mat[ina, 2]
          info <- rbind(info, mat[ina, ] )
          sela <- info[, 1]
          mat <- mat[-ina, , drop = FALSE]
        }

      }

    }

    duration <- proc.time() - durat

    #####################################################
    #####################################################
    ##          ##
    ##          ##  else it is linear regression
    ##          ##
    #####################################################
    #####################################################

  } else {
  
    durat <- proc.time()
      
    if ( !robust ) {  ## Non robust
	      if ( !heavy ) {
	        ci_test <- "testIndReg"
          mi <- lm( target ~ 1, y = FALSE, model = FALSE )
          ini <- BIC( mi )
	    	} else {
	    	  ci_test <- "testIndSpeedglm"
		      mi <- speedglm::speedlm(target ~ 1, data = dataset )
		      ini <- BIC(mi)
		    }  

      } else {  ## Robust
         ci_test <- "testIndReg"
         mi = MASS::rlm( target ~ 1, maxit = 2000, weights = wei )
         ini = BIC(mi)
      }

      if (ncores <= 1) {
	
      if ( !robust ) {  ## Non robust
	      if ( !heavy ) {
          for (i in 1:p) {
            mi <- lm( target ~ dataset[, i], y = FALSE, model = FALSE )
            bico[i] <- BIC( mi )
          }
        } else {
          for (i in 1:p) {
            mi <- speedglm::speedlm( target ~ dataset[, i], data = dataset[, i] )
            bico[i] <- BIC(mi)
          }		
		    }
      } else {  ## Robust
        for (i in 1:p) {
         mi <- MASS::rlm( target ~ dataset[, i], maxit = 2000, method = "MM")
         bico[i] <- BIC(mi)
        }
      }

     mat <- cbind(1:p, bico)

    } else {
      if ( !robust ) {  ## Non robust
        if( !heavy ) {
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          bico <- numeric(p)
          mod <- foreach( i = 1:p, .combine = rbind ) %dopar% {
            ww <- lm( target ~ dataset[, i], weights = wei, y = FALSE, model = FALSE )
            bico[i] <- BIC( ww )
          }
          stopCluster(cl)
        } else {
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          bico <- numeric(p)
          mod <- foreach( i = 1:p, .combine = rbind, .export = "speedlm", .packages = "speedglm" ) %dopar% {
            ww <- speedglm::speedlm( target ~ dataset[, i], data = dataset, weights = wei )
            bico[i] <- BIC(ww)
          }
          stopCluster(cl)		
		}
	   
      } else {  ## Robust
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        bico <- numeric(p)
        mod <- foreach( i = 1:p, .combine = rbind, export = "rlm", .packages = "MASS" ) %dopar% {
          ww <- MASS::rlm( target ~ dataset[, i], maxit = 2000, method = "MM")
          bico[i] = BIC(ww)
        }
        stopCluster(cl)
      }
      mat <- cbind(1:p, mod)

    }

    colnames(mat) <- c("variable", "BIC")
    rownames(mat) <- 1:p
    sel <- which.min( mat[, 2] )
    sela <- sel

    if ( ini - mat[sel, 2] > tol ) {
	
      if ( !robust ) {
	      if ( !heavy ) {
          mi <- lm( target ~ dataset[, sel], weights = wei, y = FALSE, model = FALSE )
          tool[1] <- BIC( mi )
		    } else {
          mi <- speedglm::speedlm( target ~ dataset[, sel], data = dataset[, sel], weight = wei)
          tool[1] <- BIC(mi)	
		    }
      } else {
        mi <- MASS::rlm( target ~ dataset[, sel], maxit = 2000, method = "MM")
        tool[1] <- BIC(mi)
      }
      moda[[ 1 ]] <- mi
      info[1, ] <- mat[sel, ]
      mat <- mat[-sel, , drop = FALSE]
    } else  {
      info <- info  
      sela <- NULL
    }
    ######
    ###     k equals 2
    ######

    if ( length(moda) > 0 ) {

      k <- 2
      pn <- p - k + 1
      mod <- list()

      if ( ncores <= 1 ) {
	  
        if ( !robust ) {
		      if ( !heavy ) {
             for (i in 1:pn) {
               ma <- lm( target ~., data = dataset[, c(sel, mat[i, 1]) ], weights = wei, y = FALSE, model = FALSE )
               mat[i, 2] <- BIC( ma )
             }
		      } else {
             for (i in 1:pn) {
               ma <- speedglm::speedlm( target ~., data = dataset[, c(sel, mat[i, 1]) ], weights = wei )
               mat[i, 2] <- BIC(ma)
             }		   
		      } 	 
        } else {
          for (i in 1:pn) {
            ma <- MASS::rlm( target ~., data = dataset[, c(sel, mat[i, 1]) ], maxit = 2000, method = "MM")
            mat[i, 2] <- BIC( ma )
          }
        }

      } else {

       if ( !robust ) {  ## Non robust
	       if ( !heavy ) {
           cl <- makePSOCKcluster(ncores)
           registerDoParallel(cl)
           bico <- numeric(pn)
           mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
             ww <- lm( target ~., data = dataset[, c(sel, mat[i, 1]) ], weights = wei, y = FALSE, model = FALSE )
             bico[i] <- BIC( ww )
           }
           stopCluster(cl)
		     } else {
           cl <- makePSOCKcluster(ncores)
           registerDoParallel(cl)
           bico <- numeric(pn)
           mod <- foreach( i = 1:pn, .combine = rbind, .export = "speedlm", .packages = "speedglm") %dopar% {
             ww <- speedglm::speedlm( target ~., data = dataset[, c(sel, mat[i, 1]) ], weights = wei )
             bico[i] <- BIC(ww)
           }
           stopCluster(cl)		 
		     }  

       } else {  ## Robust
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          bico <- numeric(pn)
          mod <- foreach( i = 1:pn, .combine = rbind, export = "rlm", .packages = "MASS" ) %dopar% {
            ww <- MASS::rlm( target ~., data = dataset[, c(sel, mat[i, 1])], maxit = 2000, method = "MM")
            bico[i] <- BIC( ww )
          }
          stopCluster(cl)
       }

        mat[, 2] <- mod
      }

       ina <- which.min( mat[, 2] )
       sel <- mat[ina, 1]

       if ( tool[1] - mat[ina, 2] <= tol ) {
          info <- info

      } else {

        tool[2] <- mat[ina, 2]
        info <- rbind(info, mat[ina, ] )
        sela <- info[, 1]
        mat <- mat[-ina, , drop = FALSE]
      }

    }
   #########
   ####      k is greater than 2
   #########
    if ( nrow(info) > 1 ) {
      while ( k < n - 15 & tool[ k - 1 ] - tool[ k ] > tol & nrow(mat) > 0 ) {
        
        k <- k + 1
        pn <- p - k + 1

        if ( ncores <= 1 ) {
		
           if ( !robust ) {
		         if ( !heavy ) { 
               for ( i in 1:pn ) {
                 ma <- lm( target ~., data = dataset[, c(sela, mat[i, 1]) ], weights = wei, y = FALSE, model = FALSE )
                 mat[i, 2] <- BIC( ma )
               }
			       } else {
               for ( i in 1:pn ) {
                 ma <- speedglm::speedlm( target ~., data = dataset[, c(sela, mat[i, 1]) ], weights = wei )
                 mat[i, 2] <- BIC(ww)
               }			 
			       }  

           } else {
             for ( i in 1:pn ) {
               ma <- MASS::rlm( target ~., data = dataset[, c(sela, mat[i, 1]) ], maxit = 2000, method = "MM")
               mat[i, 2] <- BIC( ma )
             }
           }

        } else {
           if ( !robust ) {  ## Non robust
		         if ( !heavy ) {
               cl <- makePSOCKcluster(ncores)
               registerDoParallel(cl)
               bico <- numeric(pn)
               mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
                 ww <- lm( target ~., data = dataset[, c(sela, mat[i, 1]) ], weights = wei, y = FALSE, model = FALSE )
                 bico[i] <- BIC( ww )
               }
               stopCluster(cl)
			       } else {
               cl <- makePSOCKcluster(ncores)
               registerDoParallel(cl)
               bico <- numeric(pn)
               mod <- foreach( i = 1:pn, .combine = rbind, .export = "speedlm", .packages = "speedlm") %dopar% {
                 ww <- speedglm::speedlm( target ~., data = dataset[, c(sela, mat[i, 1]) ], weights = wei )
                 bico[i] <- BIC(ww)
               }
               stopCluster(cl)			
			       }  

           } else {  ## Robust
             cl <- makePSOCKcluster(ncores)
             registerDoParallel(cl)
             bico <- numeric(pn)
             mod <- foreach( i = 1:p, .combine = rbind, export = "rlm", .packages = "MASS" ) %dopar% {
               ww <- MASS::rlm( target ~., data = dataset[, c(sela, mat[i, 1]) ], maxit = 2000, method = "MM")
               bico[i] <- BIC( ww )
             }
             stopCluster(cl)
           }
          
           mat[, 2] <- mod

        }

         ina <- which.min( mat[, 2] )
         sel <- mat[ina, 1]

         if ( tool[k - 1] - mat[ina, 2]  <= tol ) {
           info <- rbind( info,  c( -10, 1e300 ) )
           tool[k] <- Inf

         } else {
 
           tool[k] <- mat[ina, 2]
           info <- rbind(info, mat[ina, ] )
           sela <- info[, 1]
           mat <- mat[-ina, , drop = FALSE]
        }

      }
    }

    duration <- proc.time() - durat
    
  }

    d <- length(sela)
    final <- NULL

    if ( d >= 1 ) {
      if ( is.null(oiko) ) {
        if ( !robust ) {
          if ( !heavy ) {
            final <- lm( target ~., data = dataset[, sela], weights = wei, y = FALSE, model = FALSE )
	        } else  final <- speedglm::speedlm( target ~., data = dataset[, sela], weights = wei )	 
        } else  final <- MASS::rlm( target ~., data = dataset[, sela], maxit = 2000, method = "MM")
      #if ( robust == FALSE ) {
     
      } else {
        if ( !heavy ) {
          final <- glm( target ~., data = dataset[, sela], family = oiko, weights = wei, y = FALSE, model = FALSE )
        } else  final <- speedglm::speedglm( target ~., data = dataset[, sela], family = oiko, weights = wei )		  
      }
      #} else {
      #  models <- final <- robust::glmRob( target ~., data = as.data.frame( xx ), family = oiko, maxit = maxit )
      #}
    }
    info <- info[1:d, , drop = FALSE]
    if ( d == 1 )  info <- matrix(info, nrow = 1)
    colnames(info) <- c( "variables", "BIC" )
    rownames(info) <- info[, 1]
    

    list( runtime = duration, mat = t(mat), info = info, ci_test = ci_test, final = final)
}
