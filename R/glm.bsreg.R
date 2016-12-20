glm.bsreg <- function(target, dataset, threshold = 0.05, wei = NULL, heavy = FALSE, robust = FALSE) {
  
  threshold <- log(threshold)
  
  dm <- dim(dataset)
  n <- dm[1]  ## sample size 
  p <- dm[2]  ## number of variables

  if ( p > n ) {
    res <- paste("The number of variables is hiher than the sample size. 
    No backward procedure was attempted")
  
  } else {


  #check for NA values in the dataset and replace them with the variable median or the mode
    if( any(is.na(dataset)) ) {
      #dataset = as.matrix(dataset);
      warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
      if (class(dataset) == "matrix")  {
        dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
      }else{
        poia <- which( is.na(dataset), arr.ind = TRUE )[2]
        for( i in poia )  {
          xi <- dataset[, i]
          if(class(xi) == "numeric")
          {                    
            xi[ which( is.na(xi) ) ] <- median(xi, na.rm = TRUE) 
          } else if ( is.factor( xi ) ) {
            xi[ which( is.na(xi) ) ] <- levels(xi)[ which.max( as.vector( table(xi) ) )]
          }
          dataset[, i] <- xi
        }
      }
    }
    
  
    ##################################
    # target checking and initialize #
    ################################## 
   
    if ( is.null( colnames(dataset) ) )  colnames(dataset) <- paste("X", 1:p, sep = "")
  
    if ( is.matrix(target)  &  NCOL(target) == 2 )  {
      
      ci_test <- "testIndBinom"
      runtime <- proc.time()
      
      wei <- target[, 2]
      y <- target[, 1] / wei
      
      tic <- proc.time()
      if ( !heavy ) {
        ini <- glm( y ~.,  data = as.data.frame(dataset), family = binomial(logit), weights = wei, y = FALSE, model = FALSE )
        tab <- drop1( ini, test = "Chisq" )
        dof <- tab[-1, 1]
        stat <- tab[-1, 4]
        
      }  else {
        ci_test <- "testIndSpeedglm"
        ini <- speedglm::speedglm( y ~.,  data = data.frame(dataset), family = binomial(logit), weights = wei )
        dofini <- length( coef(ini) )
        stat <- dof <- numeric(p)
        for (i in 1:p) {
          mod <- speedglm::speedglm( y ~.,  data = data.frame(dataset[, -i]), family = binomial(logit), weights = wei )
          stat[i] <- mod$deviance - ini$deviance
          dof[i] <- dofini - length( coef(mod) ) 
        }
      }	
      
      mat <- cbind(1:p, pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
      colnames(mat) <- c("variable", "log.p-values", "statistic" )
      rownames(mat) <- 1:p 
      
      sel <- which.max( mat[, 2] )
      info <- matrix( c(0, -10, -10) , ncol = 3 )
      sela <- sel 
      
      if ( mat[sel, 2] < threshold ) {
        res <- list(mat = mat, final = ini  ) 
        
      } else {
        
        info[1, ] <- mat[sel, ]
        mat <- mat[-sel, ] 
        if ( !is.matrix(mat) )   mat <- matrix(mat, ncol = 3) 
        dat <- as.data.frame( dataset[, -sel] ) 
        
      } 
      
      i <- 1  
      
      if ( info[1, 2] > threshold ) {
        
        if ( !heavy ) {
          
          while ( info[i, 2] > threshold  &  NCOL(dat) > 0 )  {   
            
            ini <- glm( y ~., data = data.frame(dat), family = binomial(logit), weights = wei, y = FALSE, model = FALSE )
            i <- i + 1
            k <- p - i + 1
            
            if ( k == 1 ) {
              mod <- glm(target ~ 1, family = binomial(logit), weights = wei)
              stat <- 2 * abs( logLik(ini) - logLik(mod) )
              dof <- length( coef(ini) ) - length( coef(mod) ) 
              pval <- pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE)
              
              if (pval > threshold ) {
                final <- "No variables were selected"
                info <- rbind(info, c(mat[, 1], stat, pval) )
                dat <- as.data.frame( dataset[, -info[, 1] ] )
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
               mat <- mat[-sel, ] 
               if ( !is.matrix(mat) ) {
                 mat <- matrix(mat, ncol = 3) 
               }
               dat <- as.data.frame( dataset[, -info[, 1] ] )
            }
            
           }  
            
          }
          
        } else {
          
          while ( info[i, 2] > threshold  &  NCOL(dat) > 0 )  {   
            
            i <- i + 1       
            k <- p - i + 1
            ini <- speedglm::speedglm( y ~., data = data.frame(dat), family = binomial(logit), weights = wei )
            dofini <- length( coef(ini) ) 
            
            if ( k == 1 ) {
              mod <- speedglm::speedglm(target ~ 1, data = data.frame(dat), family = binomial(logit), weights = wei)
              stat <- 2 * abs( logLik(ini) - logLik(mod) )
              dof <- dofini - length( coef(mod) ) 
              pval <- pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE)
              
              if (pval > threshold ) {
                final <- "No variables were selected"
                info <- rbind(info, c(mat[, 1], stat, pval) )
                dat <- as.data.frame( dataset[, -info[, 1] ] )
                mat <- NULL
              } else {
                info <- rbind(info, c(0, -10, -10)) 
                final <- ini
                mat[, 2:3] <- c(pval, stat)
              }
              
            } else {
              
             stat <- dof <- numeric(k)
            
             for (j in 1:k) {
              mod <- speedglm::speedglm( y ~., data = data.frame(dat[, -j]), family = binomial(logit), weights = wei )
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
              mat <- mat[-sel, ] 
              if ( !is.matrix(mat) )  mat <- matrix(mat, ncol = 3) 
              dat <- as.data.frame( dataset[, -info[, 1] ] )
            }
            
           }	
         }    
          
        }
        
      }
      
      runtime <- proc.time() - tic		
      
    } else {  

   ############ 
   ###  Poisson or logistic regression
   ############
     
   	la <- length( Rfast::sort_unique(target) ) 
	
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
        ini <- glm( target ~.,  data = as.data.frame(dataset), family = oiko, weights = wei, y = FALSE, model = FALSE )
        tab <- drop1( ini, test = "Chisq" )
        dof <- tab[-1, 1]
        stat <- tab[-1, 4]
		
	  }  else {
	     ci_test <- "testIndSpeedglm"
       ini <- speedglm::speedglm( target ~.,  data = as.data.frame(dataset), family = oiko, weights = wei )
		   dofini <- length( coef(ini) )
	     stat <- dof <- numeric(p)
	     for (i in 1:p) {
		    mod <- speedglm::speedglm( target ~.,  data = as.data.frame(dataset[, -i]), family = oiko, weights = wei )
		    stat[i] <- mod$deviance - ini$deviance
		    dof[i] <- dofini - length( coef(mod) ) 
		   }
	  }	
	  
      mat <- cbind(1:p, pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
      colnames(mat) <- c("variable", "log.p-values", "statistic" )
      rownames(mat) <- 1:p 

      sel <- which.max( mat[, 2] )
      info <- matrix( c(0, -10, -10) , ncol = 3 )
      sela <- sel 

      if ( mat[sel, 2] < threshold ) {
        res <- list(mat = mat, final = ini  ) 
    
      } else {
       
        info[1, ] <- mat[sel, ]
        mat <- mat[-sel, ] 
        if ( !is.matrix(mat) )   mat <- matrix(mat, ncol = 3) 
        dat <- as.data.frame( dataset[, -sel] ) 

      } 
    
      i <- 1  

      if ( info[1, 2] > threshold ) {
        
       if ( !heavy ) {
	   
         while ( info[i, 2] > threshold  &  NCOL(dat) > 0 )  {   
           i <- i + 1
           k <- p - i + 1

           ini <- glm( target ~., data = data.frame(dat), family = oiko, weights = wei, y = FALSE, model = FALSE )
           
           if ( k == 1 ) {
             mod <- glm(target ~ 1, data = data.frame(dat), family = oiko, weights = wei)
             stat <- 2 * abs( logLik(ini) - logLik(mod) )
             dof <- length( coef(ini) ) - length( coef(mod) ) 
             pval <- pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE)
             
             if (pval > threshold ) {
               final <- "No variables were selected"
               info <- rbind(info, c(mat[, 1], stat, pval) )
               dat <- as.data.frame( dataset[, -info[, 1] ] )
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
               mat <- mat[-sel, ] 
               if ( !is.matrix(mat) )   mat <- matrix(mat, ncol = 3) 
               dat <- as.data.frame( dataset[, -info[, 1] ] )
             }
            
           }  

         }
         
	   } else {

         while ( info[i, 2] > threshold  &  NCOL(dat) > 0 )  {   

           i <- i + 1       
           k <- p - i + 1
           ini <- speedglm::speedglm( target ~., data = data.frame(dat), family = oiko, weights = wei )
		       dofini <- length( coef(ini) ) 
		       if ( k == 1 ) {
		         mod <- speedglm::speedglm(target ~ 1, data = data.frame(dat), family = oiko, weights = wei)
		         stat <- 2 * abs( logLik(ini) - logLik(mod) )
		         dof <- dofini - length( coef(mod) ) 
		         pval <- pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE)
		         
		         if (pval > threshold ) {
		           final <- "No variables were selected"
		           info <- rbind(info, c(mat[, 1], stat, pval) )
		           dat <- as.data.frame( dataset[, -info[, 1] ] )
		           mat <- NULL
		         } else {
		           info <- rbind(info, c(0, -10, -10)) 
		           final <- ini
		           mat[, 2:3] <- c(pval, stat)
		         }  
		       } else {		       
		       
             stat <- dof <- numeric(k)
		   
		         for (j in 1:k) {
               mod <- speedglm::speedglm( target ~., data = data.frame(dat[, -j]), family = oiko, weights = wei )
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
               mat <- mat[-sel, ] 
               if ( !is.matrix(mat) )  mat <- matrix(mat, ncol = 3) 
               dat <- as.data.frame( dataset[, -info[, 1] ] )
             }

          }	   
        }
	   }
	   
    }
      
    runtime <- proc.time() - tic		


   ############ 
   ###  Linear regression
   ############
      

    } else { 
      
	    tic <- proc.time()

      ## percentages
      if ( min( target ) > 0 & max( target ) < 1 )  {  ## are they percentages?
        target <- log( target / (1 - target) )       
      }
      
	if ( !heavy ) {
	  
	    ci_test <- "testIndReg"  
      if ( !robust ) {
        ini <- lm( target ~., data = data.frame(dataset), weights = wei, y = FALSE, model = FALSE )
        df2 <- n - length( coef(ini) )
        tab <- drop1( ini, test = "F" )
        df1 <- tab[-1, 1]
        stat <- tab[-1, 5]
        mat <- cbind(1:p, pf( stat, df1, df2, lower.tail = FALSE, log.p = TRUE), stat )
      } else {
	    	mat <- matrix(ncol = 3, nrow = p)
	      mat[, 1] <- 1:p
        ini <- MASS::rlm( target ~., data = as.data.frame(dataset), weights = wei, maxit = 2000)
		    lik1 <- as.numeric( logLik(ini) )
		    dofini <- length( coef(ini) ) 
		    for (j in 1:p) {
          fit2 <- MASS::rlm( target ~., data = as.data.frame(dataset[, -j]), weights = wei, maxit = 2000 )
          mat[j, 3] <- 2 * abs( lik1 - as.numeric( logLik(fit2) ) ) 
		      dof <- dofini - length( coef(fit2) )
          mat[j, 2] = pchisq( mat[j, 3], dof, lower.tail = FALSE, log.p = TRUE )  
		    }  
	    }
	  
    } else {
      ci_test <- "testIndSpeedglm"
	    mat <- matrix(ncol = 3, nrow = p)
	    mat[, 1] <- 1:p
      ini <- speedglm::speedlm( target ~., data = as.data.frame(dataset), weights = wei )
      d1 <- length( coef(ini) )
	    for (j in 1:p) {
          fit2 = speedglm::speedlm( target ~., data = as.data.frame(dataset[, -j]), weights = wei )
          df1 = d1 - length( coef(fit2) )
          df2 = n - d1
          mat[j, 3] = (fit2$RSS - ini$RSS)/df1 / ( ini$RSS /df2 )
          mat[j, 2] = pf( mat[j, 3], df1, df2, lower.tail = FALSE, log.p = TRUE )
        }  		  
	  }
	  
      colnames(mat) <- c("variable", "log.p-values", "statistic" )
      rownames(mat) <- 1:p 

      sel <- which.max( mat[, 2] )
      info <- matrix( c(0, -10, -10) , ncol = 3 )
      sela <- sel 

      if ( mat[sel, 2] < threshold ) {
        res <- list(mat = mat, final = ini ) 
        
      } else {
        info[1, ] <- mat[sel, ]
        mat <- mat[-sel, ] 
        if ( !is.matrix(mat) )   mat <- matrix(mat, ncol = 3) 
        dat <- as.data.frame( dataset[, -sel] ) 
      } 
    
      i <- 1  

      if ( info[1, 2] > threshold ) {

        if ( !heavy ) {
      
          while ( info[i, 2] > threshold  &  NCOL(dat) > 0 )  {   
            i <- i + 1
            k <- p - i + 1
            
	        if ( !robust ) {
            ini <- lm( target ~.,  data = as.data.frame(dat), weights = wei, y = FALSE, model = FALSE )
            df2 <- n - length( coef(ini) )
            tab <- drop1( ini, test = "F" )
            df1 <- tab[-1, 1]
            stat <- tab[-1, 5]
            pval <- pf( stat, df1, df2, lower.tail = FALSE, log.p = TRUE)
            
            if ( k == 1 ) {
              if (pval > threshold ) {
                final <- "No variables were selected"
                info <- rbind(info, c(mat[, 1], stat, pval) )
                dat <- as.data.frame( dataset[, -info[, 1] ] )
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
                mat <- mat[-sel, ] 
                if ( !is.matrix(mat) )   mat <- matrix(mat, ncol = 3) 
                dat <- as.data.frame( dataset[, -info[, 1]] )
              }
            }
            
			    } else {
            ini <- MASS::rlm( target ~., data = data.frame(dat), weights = wei, maxit = 2000 )
			      lik1 <- as.numeric( logLik(ini) )
            dofini <- length( coef(ini) )		
            
            if ( k == 1 ) {
              mod <- MASS::rlm(target ~ 1, data = data.frame(dat), weights = wei, maxit = 2000)
              stat <- 2 * abs( lik1 - logLik(mod) )
              dof <- dofini - length( coef(mod) ) 
              pval <- pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE)
              
              if (pval > threshold ) {
                final <- "No variables were selected"
                info <- rbind(info, c(mat[, 1], stat, pval) )
                dat <- as.data.frame( dataset[, -info[, 1] ] )
                mat <- NULL
              } else {
                info <- rbind(info, c(0, -10, -10)) 
                final <- ini
                mat[, 2:3] <- c(pval, stat)
              }  
            } else { 
            
			        for (j in 1:k) {
                fit2 <- MASS::rlm( target ~., data = data.frame(dat[, -j]), weights = wei, maxit = 2000 )
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
                mat <- mat[-sel, ] 
                if ( !is.matrix(mat) )   mat <- matrix(mat, ncol = 3) 
                dat <- as.data.frame( dataset[, -info[, 1]] )
              }
			      }
          } 
         }
          
        } else {

          while ( info[i, 2] > threshold  &  ncol(dat) > 0 )  {   
            
            i <- i + 1
            k <- p - i + 1
            ini <- speedglm::speedlm( target ~.,  data = data.frame(dat), weights = wei )	
            d1 <- length( coef(ini) )
            
            if ( k == 1 ) {
              mod <- speedglm::speedlm(target ~ 1, data = data.frame(dat), weights = wei)
              stat <- (mod$RSS - ini$RSS)/d1 / ( ini$RSS /(n - d1) )
              pval <- pf( stat, d1 - 1, n - d1, lower.tail = FALSE, log.p = TRUE )
              
              if (pval > threshold ) {
                final <- "No variables were selected"
                info <- rbind(info, c(mat[, 1], stat, pval) )
                dat <- as.data.frame( dataset[, -info[, 1] ] )
                mat <- NULL
              } else {
                info <- rbind(info, c(0, -10, -10)) 
                final <- ini
                mat[, 2:3] <- c(pval, stat)
              }  
            } else {
              
   	          for (j in 1:k) {
                fit2 = speedglm::speedlm( target ~., data = data.frame(dat[, -j]), weights = wei )
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
                mat <- mat[-sel, ] 
                if ( !is.matrix(mat) )   mat <- matrix(mat, ncol = 3) 
                dat <- as.data.frame( dataset[, -info[, 1]] )
              }

            } 		
          }
          
        }    
      }	

      runtime <- proc.time() - tic 
	  
     }  
    
    }
    
    info <- info[ info[, 1] > 0, ]
    res <- list(runtime = runtime, info = info, mat = mat, ci_test = ci_test, final = final ) 
    
  }  
	
  res

}    

     


 
  
    