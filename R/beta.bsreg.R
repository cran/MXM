beta.bsreg <- function(target, dataset, threshold = 0.05, wei = NULL) {
  
  threshold <- log(threshold)
  
  dm <- dim(dataset)
  n <- dm[1]  ## sample size 
  p <- dm[2]  ## number of variables
  
  if ( p > n ) {
    res <- paste("The number of variables is hiher than the sample size. No backward procedure was attempted")
    
  } else {
    
    tic <- proc.time()
    
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
    
    if ( is.null( colnames(dataset) ) )   colnames(dataset) <- paste("X", 1:p, sep = "")

      ###################
      ###################
      
      ini <- beta.mod( target,  dataset, wei = wei )
      dofini <- length( ini$be )
      likini <- ini$loglik 
      stat <- dof <- numeric(p)
      
      for (j in 1:p) {
        mod <- beta.reg( target, dataset[, -j], wei = wei )
        stat[j] <- 2 * abs( likini - mod$loglik )
        dof[j] <- dofini - length( mod$be )
      }
      
      mat <- cbind(1:p, pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
      colnames(mat) <- c("variable", "log.p-value", "statistic" )
      rownames(mat) <- 1:p 
      
      sel <- which.max( mat[, 2] )
      info <- matrix( c(0, -10, -10) , ncol = 3 )
      sela <- sel 
      
      if ( mat[sel, 2] < threshold ) {
        final <- ini 
        
      } else {
        info[1, ] <- mat[sel, ]
        mat <- mat[-sel, ] 
        if ( !is.matrix(mat) )   mat <- matrix(mat, ncol = 3) 
        dat <- as.data.frame( dataset[, -sel] ) 
      } 
      
      i <- 1  
      
      if ( info[1, 2] > threshold ) {
        
        while ( info[i, 2] > threshold  &  NCOL(dat) > 0 )  {   
          
          ini <- beta.mod( target, dat, wei = wei )
          likini <- ini$loglik
          dofini <- length(ini$be)
          
          i <- i + 1        
          k <- p - i + 1
          
          if ( k == 1 ) {
            if ( is.null(wei) ) {
              mod <- Rfast::beta.mle(target)
            } else mod <- betamle.wei(target, wei)
            
            stat <- 2 * abs( likini - mod$loglik )
            dof <- dofini - length( mod$be ) 
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
              mod <- beta.reg( target,  dat[, -j], wei = wei )
              stat[j] <- 2 * abs( likini - mod$loglik )
              dof[j] <- dofini - length( mod$be ) 
            }
            
            mat[, 2:3] <- cbind( pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
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
      
      runtime <- proc.time() - tic		
      info <- info[ info[, 1] > 0, ]
      res <- list(runtime = runtime, info = info, ci_test = "testIndBeta", final = final ) 
    }
    
  res
  
} 
  






