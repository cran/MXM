glm.bsreg = function(target, dataset, threshold = 0.05) {
  
  threshold <- log(threshold)
  n <- nrow(dataset)  ## sample size 
  p <- ncol(dataset)  ## number of variables

  if ( p > n ) {
    res <- paste("The number of variables is hiher than the sample size. 
    No backward procedure was attempted")
  
  } else {


    #check for NA values in the dataset and replace them with the variable median or the mode
    if(any(is.na(dataset)) == TRUE)
    {
      #dataset = as.matrix(dataset);
      warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    
      if(class(dataset) == "matrix")  {
        dataset = apply( dataset, 2, function(x){ x[ which(is.na(x)) ] = median(x,na.rm = TRUE) } );
      } else {
         for(i in 1:ncol(dataset)) {
          if ( any( is.na(dataset[, i]) ) )  {
            xi = dataset[,i]
            if(class(xi) == "numeric")  {                    
              xi[ which(is.na(xi)) ] = median(xi,na.rm = TRUE) 
            } else if( class(xi) == "factor" ) {
              xi[ which(is.na(xi)) ] = levels(xi)[ which.max(xi) ]
            }
            dataset[, i] = xi
          }
        }
      }
    }
  
    ##################################
    # target checking and initialize #
    ################################## 

    dataset <- as.data.frame(dataset)  
    
    if ( is.null( colnames(dataset) ) )  {  ## checks for column names
      colnames(x) <- paste("X", 1:p, sep = "")
    }	
  

   ############ 
   ###  Poisson or logistic regression
   ############
   
    if ( length( unique(target) ) == 2  || ( length( unique(target) ) > 2  &  sum( round(target) - target ) == 0 ) ) {

	    tic <- proc.time()

      if ( length( unique(target) ) == 2 ) {
        oiko = "binomial"

      } else if ( length( unique(target) ) > 2  &  sum( round(target) - target ) == 0 ) {
        oiko = "poisson"
      }
   
      mod <- glm( target ~., data = dataset, family = oiko )
      tab <- drop1( mod, test = "Chisq" )
      dof <- tab[-1, 1]
      stat <- tab[-1, 4]
      mat <- cbind(1:p, pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
      colnames(mat) <- c("variable", "p-value", "statistic" )
      rownames(mat) <- 1:p 

      sel <- which.max( mat[, 2] )
      info <- matrix( c(0, 10, 10) , ncol = 3 )
      sela <- sel 

      if ( mat[sel, 2] > threshold ) {
        final = paste("No significant variable found")
        res <- list(mat = mat, final = final  ) 
    
      } else {
       
        info[1, ] <- mat[sel, ]
        mat <- mat[-sel, ] 
        if ( !is.matrix(mat) ) {
          mat <- matrix(mat, ncol = 3) 
        }
        dat <- as.data.frame( dataset[, -sel] ) 

      } 
    
      i <- 1  

      if ( info[1, 2] < threshold ) {
        
      
       while ( info[i, 2] < threshold )  {   
      
        mod <- glm( target ~., data = dat, family = oiko )
        tab <- drop1( mod, test = "Chisq" )
        dof <- tab[-1, 1]
        stat <- tab[-1, 4]
        mat[, 2:3] <- cbind( pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
        sel <- which.max( mat[, 2] )
        i <- i + 1

        if ( mat[sel, 2] > threshold ) {
          final = summary( glm( target ~., data = as.data.frame( dataset[, info[, 1] ] ), family = oiko ) )
          res <- list(mat = mat, final = summary(final) ) 
          info <- rbind(info, c(0, 10, 10) )

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
      
      runtime <- proc.time() - tic		
      res <- list(mat = mat, final = final, runtime = runtime ) 


   ############ 
   ###  Linear regression
   ############
      

    } else { 
      
	    tic <- proc.time()
	   
      ## multivariate data
      if ( sum( class(target) == "matrix" ) > 0 ) {
        a <- rowSums(target)
        if ( min( target ) > 0 & round( sd(a), 16 ) == 0 ) { ## are they compositional data?
          target <- log( target[, -1] / target[, 1] )
        }
      }
    
      ## percentages
      if ( min( target ) > 0 & max( target ) < 1 )  {  ## are they percentages?
        target <- log( target / (1 - target) )       
      }
      
      mod <- lm( target ~., dataset )
      df2 <- n - length( coef(mod) )
      tab <- drop1( mod, test = "F" )
      df1 <- tab[-1, 1]
      stat <- tab[-1, 5]
      mat <- cbind(1:p, pf( stat, df1, df2, lower.tail = F, log.p = T), stat )
      colnames(mat) <- c("variable", "p-value", "statistic" )
      rownames(mat) <- 1:p 

      sel <- which.max( mat[, 2] )
      info <- matrix( c(0, 10, 10) , ncol = 3 )
      sela <- sel 

      if ( mat[sel, 2] > threshold ) {
        final = paste("No significant variable found")
        res <- list(mat = mat, final = final ) 
        
      } else {
       
        info[1, ] <- mat[sel, ]
        mat <- mat[-sel, ] 
        if ( !is.matrix(mat) ) {
          mat <- matrix(mat, ncol = 3) 
        }
        dat <- as.data.frame( dataset[, -sel] ) 

      } 
    
      i <- 1  

      if ( info[1, 2] < threshold ) {
        
       while ( info[i, 2] < threshold )  {   
      
        mod <- lm( target ~., dat )
        df2 <- n - length( coef(mod) )
        tab <- drop1( mod, test = "F" )
        df1 <- tab[-1, 1]
        stat <- tab[-1, 5]
        mat[, 2:3] <- cbind( pf( stat, df1, df2, lower.tail = F, log.p = T), stat )
        sel <- which.max( mat[, 2] )
        i <- i + 1

        if ( mat[sel, 2] > threshold ) {
          final = summary( lm( target ~., data = as.data.frame( dataset[, info[, 1] ] ) ) )
          res <- list(mat = mat, final = final )  
          info <- rbind(info, c(0, 10, 10) )

        } else {
      
          info <- rbind(info, mat[sel, ] )
          mat <- mat[-sel, ] 
          if ( !is.matrix(mat) ) {
            mat <- matrix(mat, ncol = 3) 
          }
          dat <- as.data.frame( dataset[, -info[, 1]] )
        }

       } 
	   
      }
      
      runtime <- proc.time() - tic 
      res <- list(mat = mat, final = final, runtime = runtime ) 
      
    } 
	
  }
  
  res  

}    

     


 
  
    