glm.bsreg = function(target, dataset, threshold = 0.05, wei = NULL, heavy = FALSE, robust = FALSE) {
  
  threshold <- log(threshold)
  
  n <- dim(dataset)[1]  ## sample size 
  p <- dim(dataset)[2]  ## number of variables

  if ( p > n ) {
    res <- paste("The number of variables is hiher than the sample size. 
    No backward procedure was attempted")
  
  } else {


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
  
    ##################################
    # target checking and initialize #
    ################################## 
   
    if ( is.null( colnames(dataset) ) )  {  ## checks for column names
      colnames(dataset) <- paste("X", 1:p, sep = "")
    }	
  

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
   
      if ( heavy == FALSE ) {
        ini <- glm( target ~.,  data = as.data.frame(dataset), family = oiko, weights = wei )
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
      colnames(mat) <- c("variable", "p-value", "statistic" )
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
        
       if ( heavy == FALSE ) {
	   
         while ( info[i, 2] > threshold  &  NCOL(dat) > 0 )  {   
         
           ini <- glm( target ~., data = as.data.frame(dat), family = oiko, weights = wei )
           tab <- drop1( ini, test = "Chisq" )
           dof <- tab[-1, 1]
           stat <- tab[-1, 4]
           mat[, 2:3] <- cbind( pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
           sel <- which.max( mat[, 2] )
           i <- i + 1
      
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
         
         if ( dim(info)[1] == p )  final <- glm( target ~ 1, family = oiko, weights = wei )
           
	   
	   } else {

         while ( info[i, 2] > threshold  &  NCOL(dat) > 0 )  {   

           i <- i + 1       
           ini <- speedglm::speedglm( target ~., data = as.data.frame(dat), family = oiko, weights = wei )
		   dofini <- length( coef(ini) ) 
		   k <- p - i 
           stat <- dof <- numeric(k)
		   
		   for (j in 1:k) {
             mod <- speedglm::speedglm( target ~., data = as.data.frame(dat[, -i]), family = oiko, weights = wei )
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
             if ( !is.matrix(mat) ) {
               mat <- matrix(mat, ncol = 3) 
             }
             dat <- as.data.frame( dataset[, -info[, 1] ] )
           }

         }	   
	   
	     if ( dim(info)[1] == p )  final <- speedglm::speedglm( target ~ 1, data = data.frame(dataset[, 1]), family = oiko, weights = wei )
	     
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
      
	if ( heavy == FALSE ) {
	  
	    ci_test <- "testIndReg"  
      if ( robust == FALSE ) {
        ini <- lm( target ~., data = as.data.frame(dataset), weights = wei )
        df2 <- n - length( coef(ini) )
        tab <- drop1( ini, test = "F" )
        df1 <- tab[-1, 1]
        stat <- tab[-1, 5]
        mat <- cbind(1:p, pf( stat, df1, df2, lower.tail = F, log.p = TRUE), stat )
      } else {
	    	mat <- matrix(ncol = 3, nrow = p)
	      mat[, 1] <- 1:p
        ini <- MASS::rlm( target ~., data = as.data.frame(dataset), weights = wei )
		    lik1 <- as.numeric( logLik(ini) )
		    dofini <- length( coef(ini) ) 
		    for (j in 1:p) {
          fit2 <- MASS::rlm( target ~., data = as.data.frame(dataset[, -i]), weights = wei, maxit = 2000 )
          mat[j, 3] <- 2 * abs( lik1 - as.numeric( logLik(fit2) ) ) 
		      dof <- dofini - length( coef(fit2) )
          mat[j, 2] = pchisq( mat[i, 3], dof, lower.tail = FALSE, log.p = TRUE )  
		    }  
	    }
	  
    } else {
      ci_test <- "testIndSpeedglm"
	    mat <- matrix(ncol = 3, nrow = p)
	    mat[, 1] <- 1:p
      ini <- speedglm::speedlm( target ~., data = as.data.frame(dataset), weights = wei )
      d1 <- length( coef(ini) )
	    for (j in 1:p) {
          fit2 = speedglm::speedlm( target ~., data = as.data.frame(dataset[, -i]), weights = wei )
          d2 = length( coef(fit2) )
          df1 = d2 - d1
          df2 = n - d2
          mat[j, 3] = (ini$RSS - fit2$RSS)/df1 / ( fit2$RSS /df2 )
          mat[j, 2] = pf( mat[i, 3], df1, df2, lower.tail = FALSE, log.p = TRUE )
        }  		  
	  
	}
	  
      colnames(mat) <- c("variable", "p-value", "statistic" )
      rownames(mat) <- 1:p 

      sel <- which.max( mat[, 2] )
      info <- matrix( c(0, -10, -10) , ncol = 3 )
      sela <- sel 

      if ( mat[sel, 2] < threshold ) {
        res <- list(mat = mat, final = ini ) 
        
      } else {
       
        info[1, ] <- mat[sel, ]
        mat <- mat[-sel, ] 
        if ( !is.matrix(mat) ) {
          mat <- matrix(mat, ncol = 3) 
        }
        dat <- as.data.frame( dataset[, -sel] ) 

      } 
    
      i <- 1  

      if ( info[1, 2] > threshold ) {

        if ( heavy == FALSE ) {
      
          while ( info[i, 2] > threshold  &  ncol(dat) > 0 )  {   
      
	        if ( robust == FALSE ) {
            ini <- lm( target ~.,  data = as.data.frame(dat), weights = wei )
            df2 <- n - length( coef(ini) )
            tab <- drop1( ini, test = "F" )
            df1 <- tab[-1, 1]
            stat <- tab[-1, 5]
            mat[, 2:3] <- cbind( pf( stat, df1, df2, lower.tail = FALSE, log.p = TRUE), stat )
			  
			    } else {
            ini <- MASS::rlm( target ~., data = as.data.frame(dat), weights = wei, maxit = 2000 )
			      lik1 <- as.numeric( logLik(ini) )
            dofini <- length( coef(ini) )			   
			      for (j in 1:p) {
              fit2 <- MASS::rlm( target ~., data = as.data.frame(dat[, -i]), weights = wei, maxit = 2000 )
              mat[j, 3] <- 2 * abs( lik1 - as.numeric( logLik(fit2) ) ) 
		          dof <- dofini - length( coef(fit2) )
              mat[j, 2] = pchisq( mat[i, 3], dof, lower.tail = FALSE, log.p = TRUE )  
		        }  
			
			    }
			
            sel <- which.max( mat[, 2] )
            i <- i + 1

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
	   
          if ( dim(info)[1] == p )  {
            if ( robust ) {
              final <- MASS::rlm( target ~., data = as.data.frame(dat), weights = wei, maxit = 2000 )
            } else  final <- lm( target ~ 1, weights = wei )
          
          } 
          
        } else {

          while ( info[i, 2] > threshold  &  ncol(dat) > 0 )  {   
      
            ini <- speedglm::speedlm( target ~.,  data = as.data.frame(dat), weights = wei )	
            d1 <- length( coef(ini) )
		        for (j in 1:p) {
              fit2 = speedglm::speedlm( target ~., data = as.data.frame(dat[, -i]), weights = wei )
              d2 = length( coef(fit2) )
              df1 = d2 - d1
              df2 = n - d2
              mat[j, 3] = (ini$RSS - fit2$RSS)/df1 / ( fit2$RSS /df2 )
              mat[j, 2] = pf( mat[i, 3], df1, df2, lower.tail = FALSE, log.p = TRUE )
            }  		  

            sel <- which.max( mat[, 2] )
            i <- i + 1

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
	  
          if ( dim(info)[1] == p )  final <- speedglm::speedlm( target ~ 1, data = data.frame(dataset[, 1]), weights = wei )
          
	     }
	   
      }	

      runtime <- proc.time() - tic 
	  
     }  
    
     info <- info[ info[, 1] > 0, ]
     res <- list(runtime = runtime, info = info, ci_test <- ci_test, final = final ) 
    
  } 
	
  res

}    

     


 
  
    