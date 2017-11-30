glm.bsreg <- function(target, dataset, threshold = 0.05, wei = NULL, heavy = FALSE, robust = FALSE) {
  
  threshold <- log(threshold)
  dm <- dim(dataset)
  if ( is.null(dm) ) {
    n <- length(target)
    p <- 1
  } else {
    n <- dm[1]  ## sample size 
    p <- dm[2]  ## number of variables
  }  
  if ( p > n ) {
    res <- paste("The number of variables is hiher than the sample size. 
    No backward procedure was attempted")
  } else {
  #check for NA values in the dataset and replace them with the variable median or the mode
    if ( any(is.na(dataset)) ) {
      #dataset = as.matrix(dataset);
      warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
      if ( is.matrix(dataset) )  {
        dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
      } else {
        poia <- unique( which( is.na(dataset), arr.ind = TRUE )[, 2] )
        for( i in poia )  {
          xi <- dataset[, i]
          if( is.numeric(xi) )   {                    
            xi[ which( is.na(xi) ) ] <- median(xi, na.rm = TRUE) 
          } else if ( is.factor( xi ) ) {
            xi[ which( is.na(xi) ) ] <- levels(xi)[ which.max( as.vector( table(xi) ) )]
          }
          dataset[, i] <- xi
        }
      }
    }
    dataset <- as.data.frame(dataset)
    ##################################
    # target checking and initialize #
    ################################## 
    if ( is.matrix(target)  &  NCOL(target) == 2 )  {
      
      ci_test <- "testIndBinom"
      wei <- target[, 2]
      ywei <- target[, 1] / wei
      
      tic <- proc.time()
      if ( !heavy ) {
        ini <- glm( ywei ~., data = dataset, family = binomial(logit), weights = wei, y = FALSE, model = FALSE )
        tab <- drop1( ini, test = "Chisq" )
        dof <- tab[-1, 1]
        stat <- tab[-1, 4]
        
      }  else {
        ci_test <- "testIndSpeedglm"
        ini <- speedglm::speedglm( ywei ~., data = dataset, family = binomial(logit), weights = wei )
        dofini <- length( coef(ini) )
        stat <- dof <- numeric(p)
		    if (p == 1) {
          mod <- speedglm::speedglm( ywei ~ 1, data = dataset, family = binomial(logit), weights = wei )
          stat <- mod$deviance - ini$deviance
          dof <- dofini - length( coef(mod) ) 
		    } else {
          for (i in 1:p) {
            mod <- speedglm::speedglm( ywei ~.,  data = dataset[, -i, drop = FALSE], family = binomial(logit), weights = wei )
            stat[i] <- mod$deviance - ini$deviance
            dof[i] <- dofini - length( coef(mod) ) 
          } 
		    } 
      }	
      
      mat <- cbind(1:p, pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
      colnames(mat) <- c("variable", "log.p-values", "statistic" )
      rownames(mat) <- 1:p 
      
      sel <- which.max( mat[, 2] )
      info <- matrix( c(0, -10, -10) , ncol = 3 )
      sela <- sel 
      
      if ( mat[sel, 2] < threshold ) {
        runtime <- proc.time() - tic
        res <- list(runtime = runtime, info = matrix(0, 0, 3), mat = mat, ci_test = ci_test, final = ini ) 
        
      } else {
        info[1, ] <- mat[sel, ]
        mat <- mat[-sel, , drop = FALSE] 
        dat <- dataset[, -sel, drop = FALSE]

      i <- 1  
      if ( info[1, 2] > threshold & dim(mat)[1] > 0 ) {
        
        if ( !heavy ) {
          
          while ( info[i, 2] > threshold  &  dim(dat)[2] > 0 )  {   
            
            ini <- glm( ywei ~., data = dat, family = binomial(logit), weights = wei, y = FALSE, model = FALSE )
            i <- i + 1
            k <- p - i + 1
            if ( k == 1 ) {
              mod <- glm(ywei ~ 1, family = binomial(logit), weights = wei)
              stat <- 2 * ( logLik(ini) - logLik(mod) )
              dof <- length( coef(ini) ) - length( coef(mod) ) 
              pval <- pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE)
              
              if (pval > threshold ) {
                final <- "No variables were selected"
                info <- rbind(info, c(mat[, 1], pval, stat) )
                dat <- dataset[, -info[, 1], drop = FALSE ] 
                mat <- NULL
              } else {
                info <- rbind(info, c(0, -10, -10)) 
                final <- ini
                mat[, 2:3] <- c(pval, stat)
              }  
            } else {
              tab <- drop1( ini, test = "Chisq" )
              dof <- tab[-1, 1]
              stat <- tab[-1, 4]
              mat[, 2:3] <- cbind( pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
              sel <- which.max( mat[, 2] )
            
             if ( mat[sel, 2] < threshold ) {
               final <- ini
               info <- rbind(info, c(0, -10, -10) )
              
             } else {
               info <- rbind(info, mat[sel, ] )
               mat <- mat[-sel, , drop = FALSE] 
               dat <- dataset[, -info[, 1], drop = FALSE ] 
            }
           }  
          }
          
        } else {
          
          while ( info[i, 2] > threshold  &  dim(dat)[2] > 0 )  {   
            
            i <- i + 1       
            k <- p - i + 1
            ini <- speedglm::speedglm( ywei ~., data = dat, family = binomial(logit), weights = wei )
            dofini <- length( coef(ini) ) 
            
            if ( k == 1 ) {
              mod <- speedglm::speedglm(ywei ~ 1, data = dat, family = binomial(logit), weights = wei)
              stat <- 2 * ( logLik(ini) - logLik(mod) )
              dof <- dofini - length( coef(mod) ) 
              pval <- pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE)
              
              if (pval > threshold ) {
                final <- "No variables were selected"
                info <- rbind(info, c(mat[, 1], pval, stat) )
                dat <- dataset[, -info[, 1], drop = FALSE ]
                mat <- NULL
              } else {
                info <- rbind(info, c(0, -10, -10)) 
                final <- ini
                mat[, 2:3] <- c(pval, stat)
              }
              
            } else {
              
             stat <- dof <- numeric(k)
             for (j in 1:k) {
              mod <- speedglm::speedglm( ywei ~., data = dat[, -j, drop = FALSE], family = binomial(logit), weights = wei )
              stat[j] <- mod$deviance - ini$deviance
              dof[j] <- dofini - length( coef(mod) ) 		   
            } 
            
            mat[, 2:3] <- cbind( pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
            sel <- which.max( mat[, 2] )

            if ( mat[sel, 2] < threshold ) {
              final <- ini
              info <- rbind(info, c(0, -10, -10) )
              
            } else {
              
              info <- rbind(info, mat[sel, ] )
              mat <- mat[-sel, , drop = FALSE] 
              dat <- dataset[, -info[, 1], drop = FALSE ]
            }
            
           }	
         }    
          
        }
        runtime <- proc.time() - tic		
        info <- info[ info[, 1] > 0, , drop = FALSE]
        res <- list(runtime = runtime, info = info, mat = mat, ci_test = ci_test, final = final ) 
        
      }  else {
        runtime <- proc.time() - tic
        if (ci_test == "testIndLogistic")   mod <- glm(ywei ~ 1, binomial)
        res <- list(runtime = runtime, info = info, mat = NULL, ci_test = ci_test, final = mod ) 
      }
      
      }   
      
    } else {  

   ############ 
   ###  Poisson or logistic regression
   ############
   	la <- length( unique(target) ) 
	
    if ( la == 2  || ( la > 2  &  sum( round(target) - target ) == 0 ) ) {

	    tic <- proc.time()

      if ( la == 2 ) {
        oiko <- binomial(logit)
        ci_test <- "testIndLogistic"
        
      } else if ( la > 2  &  sum( round(target) - target ) == 0 ) {
        oiko <- poisson(log)
        ci_test <- "testIndPois"
      }
   
      if ( !heavy ) {
        ini <- glm( target ~.,  data = dataset, family = oiko, weights = wei, y = FALSE, model = FALSE )
        tab <- drop1( ini, test = "Chisq" )
        dof <- tab[-1, 1]
        stat <- tab[-1, 4]
		
	  }  else {
	     ci_test <- "testIndSpeedglm"
       ini <- speedglm::speedglm( target ~.,  data = dataset, family = oiko, weights = wei )
		   dofini <- length( coef(ini) )
	     stat <- dof <- numeric(p)
		   if (p == 1) {
		     mod <- speedglm::speedglm( target ~ 1, data = dataset, family = oiko, weights = wei )
		     stat <- mod$deviance - ini$deviance
		     dof <- dofini - length( coef(mod) ) 		 
		   } else {
	       for (i in 1:p) {
		       mod <- speedglm::speedglm( target ~.,  data = dataset[, -i, drop = FALSE], family = oiko, weights = wei )
		       stat[i] <- mod$deviance - ini$deviance
		       dof[i] <- dofini - length( coef(mod) ) 
		     }
		   }  
	  }	
	  
      mat <- cbind(1:p, pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
      colnames(mat) <- c("variable", "log.p-values", "statistic" )
      rownames(mat) <- 1:p 

      sel <- which.max( mat[, 2] )
      info <- matrix( c(0, -10, -10) , ncol = 3 )
      sela <- sel 

      if ( mat[sel, 2] < threshold ) {
        runtime <- proc.time() - tic
        res <- list(runtime = runtime, info = matrix(0, 0, 3), mat = mat, ci_test = ci_test, final = ini ) 
        
      } else {
       
        info[1, ] <- mat[sel, ]
        mat <- mat[-sel, , drop = FALSE] 
        dat <- dataset[, -sel, drop = FALSE] 

      i <- 1  

      if ( info[1, 2] > threshold  &  dim(mat)[1] > 0 ) {
        
       if ( !heavy ) {
	   
         while ( info[i, 2] > threshold  &  dim(dat)[2] > 0 )  {   
           i <- i + 1
           k <- p - i + 1
           ini <- glm( target ~., data = dat, family = oiko, weights = wei, y = FALSE, model = FALSE )
           
           if ( k == 1 ) {
             mod <- glm(target ~ 1, data = dat, family = oiko, weights = wei)
             stat <- 2 * ( logLik(ini) - logLik(mod) )
             dof <- length( coef(ini) ) - length( coef(mod) ) 
             pval <- pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE)
             
             if (pval > threshold ) {
               final <- "No variables were selected"
               info <- rbind(info, c(mat[, 1], pval, stat) )
               dat <- dataset[, -info[, 1], drop = FALSE ]
               mat <- NULL
             } else {
               info <- rbind(info, c(0, -10, -10)) 
               final <- ini
               mat[, 2:3] <- c(pval, stat)
             }
             
           } else {
           
             tab <- drop1( ini, test = "Chisq" )
             dof <- tab[-1, 1]
             stat <- tab[-1, 4]
             mat[, 2:3] <- cbind( pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
             sel <- which.max( mat[, 2] )
      
             if ( mat[sel, 2] < threshold ) {
               final <- ini
               info <- rbind(info, c(0, -10, -10) )

             } else {
               info <- rbind(info, mat[sel, ] )
               mat <- mat[-sel, , drop = FALSE] 
               dat <- dataset[, -info[, 1], drop = FALSE ]
             }
           }  
         }
         
	   } else {

         while ( info[i, 2] > threshold  &  dim(dat)[2] > 0 )  {   

           i <- i + 1       
           k <- p - i + 1
           ini <- speedglm::speedglm( target ~., data = dat, family = oiko, weights = wei )
		       dofini <- length( coef(ini) ) 
		       if ( k == 1 ) {
		         mod <- speedglm::speedglm(target ~ 1, data = dat, family = oiko, weights = wei)
		         stat <- 2 * ( logLik(ini) - logLik(mod) )
		         dof <- dofini - length( coef(mod) ) 
		         pval <- pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE)
		         
		         if (pval > threshold ) {
		           final <- "No variables were selected"
		           info <- rbind(info, c(mat[, 1], pval, stat) )
		           dat <- dataset[, -info[, 1], drop = FALSE ]
		           mat <- NULL
		         } else {
		           info <- rbind(info, c(0, -10, -10)) 
		           final <- ini
		           mat[, 2:3] <- c(pval, stat)
		         }  
		       } else {		       
		       
             stat <- dof <- numeric(k)
		   
		         for (j in 1:k) {
               mod <- speedglm::speedglm( target ~., data = dat[, -j, drop = FALSE], family = oiko, weights = wei )
		           stat[j] <- mod$deviance - ini$deviance
		           dof[j] <- dofini - length( coef(mod) ) 		   
		        } 

             mat[, 2:3] <- cbind( pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
             sel <- which.max( mat[, 2] )
      
             if ( mat[sel, 2] < threshold ) {
               final <- ini
               info <- rbind(info, c(0, -10, -10) )

             } else {
               info <- rbind(info, mat[sel, ] )
               mat <- mat[-sel, , drop = FALSE] 
               dat <- dataset[, -info[, 1], drop = FALSE]
             }

          }	   
        }
	   }
      runtime <- proc.time() - tic		
      info <- info[ info[, 1] > 0, , drop = FALSE]
      res <- list(runtime = runtime, info = info, mat = mat, ci_test = ci_test, final = final ) 
        
    } else {
      runtime <- proc.time() - runtime
      if (ci_test == "testIndPois")       mod <- glm(target ~ 1, poisson)
      if (ci_test == "testIndLogistic")    mod <- glm(target ~ 1, binomial)
      res <- list(runtime = runtime, info = info, mat = NULL, ci_test = ci_test, final = mod ) 
    }
      
    }
    ############ 
    ###  Linear regression
    ############
    } else { 
      
	    tic <- proc.time()

	  if ( !heavy ) {
	  
	    ci_test <- "testIndReg"  
      if ( !robust ) {
        ini <- lm( target ~., data = dataset, weights = wei, y = FALSE, model = FALSE )
        df2 <- n - length( coef(ini) )
        tab <- drop1( ini, test = "F" )
        df1 <- tab[-1, 1]
        stat <- tab[-1, 5]
        mat <- cbind(1:p, pf( stat, df1, df2, lower.tail = FALSE, log.p = TRUE), stat )
      } else {
	    	mat <- matrix(ncol = 3, nrow = p)
	      mat[, 1] <- 1:p
        ini <- MASS::rlm( target ~., data = dataset, weights = wei, maxit = 2000)
		    lik1 <- as.numeric( logLik(ini) )
		    dofini <- length( coef(ini) ) 
			  if (p == 1) {
			    fit2 <- MASS::rlm( target ~ 1, data = dataset, weights = wei, maxit = 2000 )
          mat[1, 3] <- 2 * abs( lik1 - as.numeric( logLik(fit2) ) ) 
		      dof <- dofini - length( coef(fit2) )
          mat[1, 2] = pchisq( mat[1, 3], dof, lower.tail = FALSE, log.p = TRUE )  
			  } else {
		      for (j in 1:p) {
            fit2 <- MASS::rlm( target ~., data = dataset[, -j, drop = FALSE], weights = wei, maxit = 2000 )
            mat[j, 3] <- 2 * abs( lik1 - as.numeric( logLik(fit2) ) ) 
		        dof <- dofini - length( coef(fit2) )
            mat[j, 2] = pchisq( mat[j, 3], dof, lower.tail = FALSE, log.p = TRUE )  
		      }
        }			
	    }
	  
    } else {
      ci_test <- "testIndSpeedglm"
	    mat <- matrix(ncol = 3, nrow = p)
	    mat[, 1] <- 1:p
      ini <- speedglm::speedlm( target ~., data = dataset, weights = wei )
      d1 <- length( coef(ini) )
	  if (p == 1) {
	      fit2 = speedglm::speedlm( target ~ 1, data = dataset, weights = wei )
          df1 = d1 - length( coef(fit2) )
          df2 = n - d1
          mat[1, 3] = (fit2$RSS - ini$RSS)/df1 / ( ini$RSS /df2 )
          mat[1, 2] = pf( mat[1, 3], df1, df2, lower.tail = FALSE, log.p = TRUE )
	  } else {
	    for (j in 1:p) {
          fit2 = speedglm::speedlm( target ~., data = dataset[, -j, drop = FALSE], weights = wei )
          df1 = d1 - length( coef(fit2) )
          df2 = n - d1
          mat[j, 3] = (fit2$RSS - ini$RSS)/df1 / ( ini$RSS /df2 )
          mat[j, 2] = pf( mat[j, 3], df1, df2, lower.tail = FALSE, log.p = TRUE )
        }  
       }		
	  }
	  
      colnames(mat) <- c("variable", "log.p-values", "statistic" )
      rownames(mat) <- 1:p 

      sel <- which.max( mat[, 2] )
      info <- matrix( c(0, -10, -10) , ncol = 3 )
      sela <- sel 

      if ( mat[sel, 2] < threshold ) {
        runtime <- proc.time() - tic
        res <- list(runtime = runtime, info = matrix(0, 0, 3), mat = mat, ci_test = ci_test, final = ini ) 
        
      } else {
        info[1, ] <- mat[sel, ]
        mat <- mat[-sel, , drop = FALSE] 
        dat <- dataset[, -sel, drop = FALSE] 

      i <- 1  

      if ( info[1, 2] > threshold & dim(mat)[1] > 0 ) {

        if ( !heavy ) {
      
          while ( info[i, 2] > threshold  &  dim(dat)[2] > 0 )  {   
            i <- i + 1
            k <- p - i + 1
            
	        if ( !robust ) {
            ini <- lm( target ~.,  data = dat, weights = wei, y = FALSE, model = FALSE )
            df2 <- n - length( coef(ini) )
            tab <- drop1( ini, test = "F" )
            df1 <- tab[-1, 1]
            stat <- tab[-1, 5]
            pval <- pf( stat, df1, df2, lower.tail = FALSE, log.p = TRUE)
            
            if ( k == 1 ) {
              if (pval > threshold ) {
                final <- "No variables were selected"
                info <- rbind(info, c(mat[, 1], pval, stat) )
                dat <- dataset[, -info[, 1], drop = FALSE ]
                mat <- NULL
              } else {
                info <- rbind(info, c(0, -10, -10)) 
                final <- ini
                mat[, 2:3] <- c(pval, stat)
              }  
            } else {
              
              mat[, 2:3] <- cbind( pval, stat )
              
              sel <- which.max( mat[, 2] )
              if ( mat[sel, 2] < threshold ) {
                final <- ini
                info <- rbind(info, c(0, -10, -10) )
                
              } else {
                
                info <- rbind(info, mat[sel, ] )
                mat <- mat[-sel, , drop = FALSE] 
                dat <- dataset[, -info[, 1], drop = FALSE] 
              }
            }
            
			    } else {
            ini <- MASS::rlm( target ~., data = dat, weights = wei, maxit = 2000 )
			      lik1 <- as.numeric( logLik(ini) )
            dofini <- length( coef(ini) )		
            
            if ( k == 1 ) {
              mod <- MASS::rlm(target ~ 1, data = dat, weights = wei, maxit = 2000)
              stat <- 2 * abs( lik1 - logLik(mod) )
              dof <- dofini - length( coef(mod) ) 
              pval <- pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE)
              
              if (pval > threshold ) {
                final <- "No variables were selected"
                info <- rbind(info, c(mat[, 1], pval, stat) )
                dat <- dataset[, -info[, 1], drop = FALSE ]
                mat <- NULL
              } else {
                info <- rbind(info, c(0, -10, -10)) 
                final <- ini
                mat[, 2:3] <- c(pval, stat)
              }  
            } else { 
            
			        for (j in 1:k) {
                fit2 <- MASS::rlm( target ~., data = dat[, -j, drop = FALSE], weights = wei, maxit = 2000 )
                mat[j, 3] <- 2 * abs( lik1 - as.numeric( logLik(fit2) ) ) 
		            dof <- dofini - length( coef(fit2) )
                mat[j, 2] <- pchisq( mat[j, 3], dof, lower.tail = FALSE, log.p = TRUE )  
		          }  
            
              sel <- which.max( mat[, 2] )
              if ( mat[sel, 2] < threshold ) {
			          final <- ini
                info <- rbind(info, c(0, -10, -10) )

              } else {
      
                info <- rbind(info, mat[sel, ] )
                mat <- mat[-sel, , drop = FALSE] 
                dat <- dataset[, -info[, 1] ,drop = FALSE]
              }
			      }
          } 
         }
          
        } else {

          while ( info[i, 2] > threshold  &  dim(dat)[2] > 0 )  {   
            
            i <- i + 1
            k <- p - i + 1
            ini <- speedglm::speedlm( target ~.,  data = dat, weights = wei )	
            d1 <- length( coef(ini) )
            
            if ( k == 1 ) {
              mod <- speedglm::speedlm(target ~ 1, data = dat, weights = wei)
              stat <- (mod$RSS - ini$RSS)/d1 / ( ini$RSS /(n - d1) )
              pval <- pf( stat, d1 - 1, n - d1, lower.tail = FALSE, log.p = TRUE )
              
              if (pval > threshold ) {
                final <- "No variables were selected"
                info <- rbind(info, c(mat[, 1], pval, stat) )
                dat <- as.data.frame( dataset[, -info[, 1] ] )
                mat <- NULL
              } else {
                info <- rbind(info, c(0, -10, -10)) 
                final <- ini
                mat[, 2:3] <- c(pval, stat)
              }  
            } else {
              
   	          for (j in 1:k) {
                fit2 = speedglm::speedlm( target ~., data = dat[, -j, drop = FALSE], weights = wei )
                df1 = d1 - length( coef(fit2) )
                df2 = n - d1
                mat[j, 3] = (fit2$RSS - ini$RSS)/df1 / ( ini$RSS /df2 )
                mat[j, 2] = pf( mat[j, 3], df1, df2, lower.tail = FALSE, log.p = TRUE )
              }    		  

              sel <- which.max( mat[, 2] )

              if ( mat[sel, 2] < threshold ) {
                final <- ini
                info <- rbind(info, c(0, -10, -10) )

              } else {
                info <- rbind(info, mat[sel, ] )
                mat <- mat[-sel, , drop = FALSE] 
                dat <- dataset[, -info[, 1], drop = FALSE]
              }

            } 		
          }
          
        }
        
        runtime <- proc.time() - tic		
        info <- info[ info[, 1] > 0, , drop = FALSE]
        res <- list(runtime = runtime, info = info, mat = mat, ci_test = ci_test, final = final ) 
        
      }	else {
        runtime <- proc.time() - tic
        if (ci_test == "testIndReg")     mod <- lm(target ~ 1, weights = wei)
        res <- list(runtime = runtime, info = info, mat = NULL, ci_test = ci_test, final = mod ) 
      }

      }
     }  
    
    }

    
  }  
  res
}    

     


 
  
    