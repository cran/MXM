bic.glm.fsreg <- function( target, dataset, wei = NULL, tol = 0, heavy = FALSE, robust = FALSE, ncores = 1) {

  p <- ncol(dataset)  ## number of variables
  bico <- numeric( p )
  moda <- list()
  k <- 1   ## counter
  n <- length(target)  ## sample size
  tool <- NULL
  oiko <- NULL
  info <- matrix( 0, ncol = 2 )
  # cont = robust::lmRob.control(mxr = 2000, mxf = 2000, mxs = 2000 )

  
    #check for NA values in the dataset and replace them with the variable median or the mode
  if(any(is.na(dataset)) == TRUE)
  {

  
    #dataset = as.matrix(dataset);
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    
    if (class(dataset) == "matrix")  {
    
       dataset = apply(dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x)}) 
              
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
  
  if ( heavy )  con <- log(n)
  
  #########
  ## if it is binomial or poisson regression
  #########
  
  la <- length( Rfast::sort_unique(target) ) 
  
  if ( la == 2  || ( sum( round(target) - target ) == 0  &  la > 2 ) ) {

   if ( la == 2 ) {
      oiko <- binomial(logit)  ## binomial regression
      ci_test <- "testIndLogistic"
    } else {
      oiko <- poisson(log)  ## poisson regression
      ci_test <- "testIndPois"
    }
    
    durat <- proc.time()
    
    #if ( robust == FALSE ) {
	  if ( heavy == FALSE ) {
        ini = BIC( glm( target ~ 1, family = oiko, weights = wei ) )   ## initial BIC
	  } else {
	    inimod <- speedglm::speedglm( target ~ 1, family = oiko, weights = wei )
	    ini <- inimod$logLik + con   ## initial BIC
	    ci_test <- "testIndSpeedglm"
	  }	
    #} else {
    #  ini = robust::glmRob( target ~ 1, family = oiko, maxit = maxit )$deviance + log(n)  ## initial BIC
    #}
    	
    if (ncores <= 1) {
	
      #if ( robust == FALSE ) {  ## Non robust
	    if ( heavy == FALSE ) {
          for (i in 1:p) {
            mi <- glm( target ~ dataset[, i], family = oiko, weights = wei )
            bico[i] <- BIC( mi )
          }
		} else {
      for (i in 1:p) { 
        mi <- speedglm::speedglm( target ~ dataset[, i], data = as.data.frame( dataset[, i] ), family = oiko, weights = wei )
	      bico[i] <- mi$logLik + length( coef(mi) ) * con   ## initial BIC
		  }	
		}  

      # } else {  ## Robust
      #   for (i in 1:p) {
      #     mi <- robust::glmRob( target ~ dataset[, i], family = oiko, maxit = maxit )
      #     bico[i] <- mi$deviance + length( coef( mi ) ) * log(n)
      #   }
      # }
        
      mat <- cbind(1:p, bico)

      if( any( is.na(mat) ) ) {
        mat[ which( is.na(mat) ) ] = ini
      }

    } else {
      #if ( robust == FALSE ) {  ## Non robust
	    if ( heavy == FALSE ) {
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          bico <- numeric(p)
          mod <- foreach( i = 1:p, .combine = rbind) %dopar% {
            ww <- glm( target ~ dataset[, i], family = oiko, weights = wei )
            bico[i] <- BIC( ww )
          }
          stopCluster(cl)
		} else {
		  cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          bico <- numeric(p)
          mod <- foreach( i = 1:p, .combine = rbind, .export = "speedglm", .packages = "speedglm") %dopar% {
            ww <- speedglm( target ~ dataset[, i], data = as.data.frame(dataset[, i]), family = oiko, weights = wei )
 	        bico[i] <- ww$logLik + length( coef(ww) ) * con   ## initial BIC

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

      if ( any( is.na(mat) ) ) {
        mat[ which( is.na(mat) ) ] = ini
      }
      
    }

    colnames(mat) <- c("variable", "BIC")
    rownames(mat) <- 1:p
    sel <- which.min( mat[, 2] )
    sela <- sel

    if ( mat[sel, 2] < ini ) {

      info[1, ] <- mat[sel, ]
      mat <- mat[-sel, ]
      if ( !is.matrix(mat) ) {
        mat <- matrix(mat, ncol = 2) 
      }
      mat <- mat[ order( mat[, 2] ), ]

      #if ( robust == FALSE ) {

	    if ( heavy == FALSE ) {
          mi <- glm( target ~ dataset[, sel], family = oiko, weights = wei )
          tool[1] <- BIC( mi )
		   } else {
          mi <- speedglm::speedglm( target ~ dataset[, sel], data = as.data.frame(dataset[, sel]), family = oiko, weights = wei )
          tool[1] <- mi$logLik + length( coef(mi) ) * con   ## initial BIC
		   }  
      #} else {
      #  mi = robust::glmRob( target ~ dataset[, sel], family = oiko, maxit = maxit )
      #  tool[1] <- mi$deviance + length( coef( mi ) ) * log(n)
      #  if ( is.na(tool[1]) )  tool[1] <- ini
      #}

      moda[[ 1 ]] <- mi
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
		    if ( heavy == FALSE ) {
          for ( i in 1:pn ) {
            ma <- glm( target ~., data = as.data.frame( dataset[, c(sel, mat[i, 1]) ] ), family = oiko, weights = wei )
            bico[i] <- BIC( ma )
          }
        } else {
          for ( i in 1:pn ) {
            ma <- speedglm::speedglm( target ~., data = as.data.frame( dataset[, c(sel, mat[i, 1]) ] ), family = oiko, weights = wei )
            bico[i] <- ma$logLik + length( coef(ma) ) * con
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
		  if ( heavy == FALSE ) {
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
              ww <- speedglm( target ~ dataset[, sel ] + dataset[, mat[i, 1] ], data = as.data.frame(dataset), family = oiko, weights = wei )
              bico[i] <- ww$logLik + length( coef(ww) ) * con
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

      if ( mat[ina, 2] >= tool[1] ) {
        info <- rbind( info,  c( -10, 1e300 ) )

      } else {
        tool[2] <- mat[ina, 2]
        info <- rbind(info, mat[ina, ] )
        sela <- info[, 1]
        mat <- mat[-ina , ]
        if ( !is.matrix(mat) ) {
          mat <- matrix(mat, ncol = 2) 
        }
        mat <- mat[ order( mat[, 2] ), ]

      }

     if ( sum( info[2, ] -  c( -10, 1e300 ) ) == 0 ) {
       info <- matrix( info[1, ], ncol = 2)
     }

   }


      #########
      ####      k is greater than 2
      #########


    if ( nrow(info) > 1  &  nrow(mat) > 0 ) {
      while ( ( k < n - 10 ) & ( tool[ k - 1 ] - tool[ k ] > tol ) & ( nrow(mat) > 0 ) ) {

        k <- k + 1
        pn <- p - k + 1

        if (ncores <= 1) {
          #if ( robust == FALSE ) {  ## Non robust
		  if ( heavy == FALSE ) {
            for ( i in 1:pn ) {
              ma <- glm( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), family = oiko, weights = wei )
              mat[i, 2] <- BIC( ma )
            }
          } else {
            for ( i in 1:pn ) {
              ma <- speedglm::speedglm( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), family = oiko, weights = wei )
              mat[i, 2] <- ma$logLik + length( coef(ma) ) * con
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
		      if ( heavy == FALSE ) {
            cl <- makePSOCKcluster(ncores)
            registerDoParallel(cl)
            bico <- numeric(pn)
            mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
              ww <- glm( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), family = oiko, weights = wei )
              bico[i] <- BIC( ww )
            }
            stopCluster(cl)
          } else {
            cl <- makePSOCKcluster(ncores)
            registerDoParallel(cl)
            bico <- numeric(pn)
            mod <- foreach( i = 1:pn, .combine = rbind, .export = "speedglm", .packages = "speedglm") %dopar% {
              ww <- speedglm( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), family = oiko, weights = wei )
              bico[i] <- ww$logLik + length( coef(ww) ) * con
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

        if ( mat[ina, 2] >= tool[k - 1] ) {
          info <- rbind( info,  c( -10, 1e300 ) )
          tool[k] <- 1e+300

        } else {

          tool[k] = mat[ina, 2]
          info <- rbind(info, mat[ina, ] )
          sela <- info[, 1]
          mat <- mat[-ina , ]
          if ( !is.matrix(mat) ) {
            mat <- matrix(mat, ncol = 2) 
          }
          mat <- mat[ order( mat[, 2] ), ]

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
      
	  if ( min(target) > 0  &  max(target) < 1 ) {
	    target <- log( target / (1 - target) )
	  
	  }
	  
      if ( robust == FALSE ) {  ## Non robust
	      if ( heavy == FALSE ) {
	        ci_test <- "testIndReg"
          mi <- lm( target ~ 1 )
          ini <- BIC( mi )
	    	} else {
	    	  ci_test <- "testIndSpeedglm"
		      mi <- speedglm::speedlm(target ~ 1)
		      ini <- mi$logLik + con
		    }  

      } else {  ## Robust
         ci_test <- "testIndReg"
         mi = MASS::rlm( target ~ 1, maxit = 2000, weights = wei )
         ini = BIC(mi)
      }

      if (ncores <= 1) {
	
      if ( robust == FALSE ) {  ## Non robust
	      if ( heavy == FALSE ) {
          for (i in 1:p) {
            mi <- lm( target ~ dataset[, i] )
            bico[i] <- BIC( mi )
          }
        } else {
          for (i in 1:p) {
            mi <- speedglm::speedlm( target ~ dataset[, i], data = as.data.frame(dataset[, i]) )
            bico[i] <- mi$logLik + length( coef(mi) ) * con
          }		
		    }
      } else {  ## Robust
        for (i in 1:p) {
         mi = MASS::rlm( target ~ dataset[, i], maxit = 2000, weights = wei )
         bico[i] = BIC(mi)
        }
      }

     mat <- cbind(1:p, bico)

    } else {
      if ( robust == FALSE ) {  ## Non robust
        if( heavy == FALSE ) {
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          bico <- numeric(p)
          mod <- foreach( i = 1:p, .combine = rbind ) %dopar% {
            ww <- lm( target ~ dataset[, i], weights = wei )
            bico[i] <- BIC( ww )
          }
          stopCluster(cl)
        } else {
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          bico <- numeric(p)
          mod <- foreach( i = 1:p, .combine = rbind, .export = "speedlm", .packages = "speedglm" ) %dopar% {
            ww <- speedglm::speedlm( target ~ dataset[, i], data = as.data.frame(dataset[, i]), weights = wei )
            bico[i] <- ww$logLik + length( coef(ww) ) * con
          }
          stopCluster(cl)		
		}
	   
      } else {  ## Robust
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        bico <- numeric(p)
        mod <- foreach( i = 1:p, .combine = rbind, export = "rlm", .packages = "MASS" ) %dopar% {
          ww <- rlm( target ~ dataset[, i], maxit = 2000, weights = wei )
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

    if ( mat[sel, 2] < ini ) {
	
      if ( robust == FALSE ) {
	      if ( heavy == FALSE ) {
          mi <- lm( target ~ dataset[, sel], weights = wei)
          tool[1] <- BIC( mi )
		    } else {
          mi <- speedglm::speedlm( target ~ dataset[, sel], data = as.data.frame(dataset[, sel]), weight = wei)
          tool[1] <- mi$logLik + length( coef(mi) ) * con	
		    }
      } else {
        mi <- MASS::rlm( target ~ dataset[, sel], maxit = 2000, weights = wei )
        tool[1] <- BIC(mi)
      }
      moda[[ 1 ]] <- mi
      info[1, ] <- mat[sel, ]
      mat <- mat[-sel, ]
      if ( !is.matrix(mat) ) {
        mat <- matrix(mat, ncol = 2) 
      }
      mat <- mat[ order( mat[, 2] ), ]    
	  
	  }

    ######
    ###     k equals 2
    ######

    if ( length(moda) > 0 ) {

      k <- 2
      pn <- p - k + 1
      mod <- list()

      if ( ncores <= 1 ) {
	  
        if ( robust == FALSE ) {
		      if ( heavy == FALSE ) {
             for (i in 1:pn) {
               ma <- lm( target ~., data = as.data.frame( dataset[, c(sel, mat[i, 1]) ] ), weights = wei )
               mat[i, 2] <- BIC( ma )
             }
		      } else {
             for (i in 1:pn) {
               ma <- speedglm::speedlm( target ~., data = as.data.frame( dataset[, c(sel, mat[i, 1]) ] ), weights = wei )
               mat[i, 2] <- ma$logLik + length( coef(ma) ) * con
             }		   
		      } 	 
        } else {
          for (i in 1:pn) {
            ma <- MASS::rlm( target ~., data = data.frame( dataset[, c(sel, mat[i, 1]) ] ), maxit = 2000, weights = wei )
            mat[i, 2] <- BIC( ma )
          }
        }

      } else {

       if ( robust == FALSE ) {  ## Non robust
	       if ( heavy == FALSE ) {
           cl <- makePSOCKcluster(ncores)
           registerDoParallel(cl)
           bico <- numeric(pn)
           mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
             ww <- lm( target ~., data = as.data.frame( dataset[, c(sel, mat[i, 1]) ] ), weights = wei )
             bico[i] <- BIC( ww )
           }
           stopCluster(cl)
		     } else {
           cl <- makePSOCKcluster(ncores)
           registerDoParallel(cl)
           bico <- numeric(pn)
           mod <- foreach( i = 1:pn, .combine = rbind, .export = "speedlm", .packages = "speedglm") %dopar% {
             ww <- speedglm::speedlm( target ~., data = as.data.frame( dataset[, c(sel, mat[i, 1]) ] ), weights = wei )
             bico[i] <- ww$logLik + length( coef(ww) ) * con
           }
           stopCluster(cl)		 
		     }  

       } else {  ## Robust
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          bico <- numeric(pn)
          mod <- foreach( i = 1:pn, .combine = rbind, export = "rlm", .packages = "MASS" ) %dopar% {
            ww <- rlm( target ~., data = data.frame( dataset[, c(sel, mat[i, 1]) ] ), maxit = 2000, weights = wei )
            bico[i] <- BIC( ww )
          }
          stopCluster(cl)
       }

        mat[, 2] <- mod
      }

       ina <- which.min( mat[, 2] )
       sel <- mat[ina, 1]

      if ( mat[ina, 2] >= tool[1] ) {
        info <- rbind( info,  c( -10, 1e300 ) )

      } else {

        tool[2] <- mat[ina, 2]
        info <- rbind(info, mat[ina, ] )
        sela <- info[, 1]
        mat <- mat[-ina , ]
        if ( !is.matrix(mat) ) mat <- matrix(mat, ncol = 2) 
        mat <- mat[ order( mat[, 2] ), ]

      }

      if ( sum( info[2, ] -  c( -10, 1e300 ) ) == 0 ) {
        info <- matrix( info[1, ], ncol = 2)
      }

    }

      #########
      ####      k is greater than 2
      #########

    if ( nrow(info) > 1 ) {
      while ( ( k < n - 10 ) & ( tool[ k - 1 ] - tool[ k ] > tol ) ) {

        k <- k + 1
        pn <- p - k + 1

        if ( ncores <= 1 ) {
		
           if ( robust == FALSE ) {
		         if ( heavy == FALSE ) { 
               for ( i in 1:pn ) {
                 ma <- lm( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), weights = wei )
                 mat[i, 2] <- BIC( ma )
               }
			       } else {
               for ( i in 1:pn ) {
                 ma <- speedglm::speedlm( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), weights = wei )
                 mat[i, 2] <- ww$logLik + length( coef(ww) ) * con
               }			 
			       }  

           } else {
             for ( i in 1:pn ) {
               ma <- MASS::rlm( target ~., data = data.frame( dataset[, c(sela, mat[i, 1]) ] ), maxit = 2000, weights = wei )
               mat[i, 2] <- BIC( ma )
             }
           }

        } else {
           if ( robust == FALSE ) {  ## Non robust
		         if ( heavy == FALSE ) {
               cl <- makePSOCKcluster(ncores)
               registerDoParallel(cl)
               bico <- numeric(pn)
               mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
                 ww <- lm( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), weights = wei )
                 bico[i] <- BIC( ww )
               }
               stopCluster(cl)
			       } else {
               cl <- makePSOCKcluster(ncores)
               registerDoParallel(cl)
               bico <- numeric(pn)
               mod <- foreach( i = 1:pn, .combine = rbind, .export = "speedlm", .packages = "speedlm") %dopar% {
                 ww <- speedglm::speedlm( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), weights = wei )
                 bico[i] <- ww$logLik + length( coef(ww) ) * con
               }
               stopCluster(cl)			
			       }  

           } else {  ## Robust
             cl <- makePSOCKcluster(ncores)
             registerDoParallel(cl)
             bico <- numeric(pn)
             mod <- foreach( i = 1:p, .combine = rbind, export = "rlm", .packages = "MASS" ) %dopar% {
               ww <- rlm( target ~., data = data.frame( dataset[, c(sela, mat[i, 1]) ] ), maxit = 2000, weights = wei )
               bico[i] = BIC( ww )
             }
             stopCluster(cl)
           }
          
           mat[, 2] <- mod

        }

         ina <- which.min( mat[, 2] )
         sel <- mat[ina, 1]

         if ( mat[ina, 2] >= tool[k - 1] ) {
           info <- rbind( info,  c( -10, 1e300 ) )
           tool[k] <- 1e+300

         } else {
 
           tool[k] <- mat[ina, 2]
           info <- rbind(info, mat[ina, ] )
           sela <- info[, 1]
           mat <- mat[-ina , ]
           if ( !is.matrix(mat) ) {
             mat <- matrix(mat, ncol = 2) 
           }
           mat <- mat[ order( mat[, 2] ), ]

        }

      }
    }

    duration <- proc.time() - durat
    
  }

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

        if ( is.null(oiko) ) {
		
           if ( robust == FALSE ) {
		     if ( heavy == FALSE ) {
               models <- final <- lm( target ~., data = as.data.frame( xx ), weights = wei )
			 } else {
               models <- final <- speedglm::speedlm( target ~., data = as.data.frame( xx ), weights = wei )	 
			 }  
           } else {
             models <- final <- MASS::rlm( target ~., data = data.frame( xx ), maxit = 2000, weights = wei )
           }

        } else {
          #if ( robust == FALSE ) {
		  if ( heavy == FALSE ) {
            models <- final <- glm( target ~., data = as.data.frame( xx ), family = oiko, weights = wei )
		  } else {
            models <- final <- speedglm::speedglm( target ~., data = as.data.frame( xx ), family = oiko, weights = wei )		  
		  }	
          #} else {
          #  models <- final <- robust::glmRob( target ~., data = as.data.frame( xx ), family = oiko, maxit = maxit )
          #}
        }

      } else {
        if ( is.null(oiko) ) {
		
          if ( robust == FALSE ) {
		         if ( heavy == FALSE ) {
               for (i in 1:d)  models[[ i ]] <- lm( target ~., data = as.data.frame( xx[, 1:i] ), weights = wei )
             } else {
               for (i in 1:d)   models[[ i ]] <- speedglm::speedlm( target ~., data = as.data.frame( xx[, 1:i] ), weights = wei )
			       }  			   
          } else {
		        for (i in 1:d)  models[[ i ]] <- MASS::rlm( target ~., data = data.frame( xx[, 1:i] ), maxit = 2000, weights = wei )
            
          }
		  
        } else {
          if ( heavy == FALSE ) {
            for (i in 1:d) {
            #if ( robust == FALSE ) {
              models[[ i ]] <- glm( target ~., data = as.data.frame( xx[, 1:i] ), family = oiko, weights = wei )
		       	}  
		    } else {
		       for (i in 1:d)  models[[ i ]] <- speedglm::speedglm( target ~., data = as.data.frame( xx[, 1:i] ), family = oiko, weights = wei )		
		  }  
            #} else {
            #  models[[ i ]] <- robust::glmRob( target ~., data = as.data.frame( xx[, 1:i] ), family = oiko, maxit = maxit )
            #}
        }

        final <- summary( models[[ d ]] )
        
      }

      info <- info[1:d, ]
      if ( d == 1 )  info <- matrix(info, nrow = 1)
      colnames(info) <- c( "variables", "BIC" )
      rownames(info) <- info[, 1]
    
    }

    list(mat = t(mat), info = info, models = models, final = final, ci_test = ci_test, runtime = duration)

}
