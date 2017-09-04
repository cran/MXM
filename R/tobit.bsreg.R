tobit.bsreg <- function(target, dataset, threshold = 0.05, wei = NULL, heavy = FALSE, robust = FALSE) {
  
  threshold <- log(threshold)
  dm <- dim(dataset)
  if ( is.null(dm) ) {
    n <- length(dataset)
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
    runtime <- proc.time()

      ini <- survival::survreg( target ~.,  data = as.data.frame(dataset), weights = wei, dist = "gaussian" )
      dofini <- length( coef(ini) )
      stat <- dof <- numeric(p)
      for (i in 1:p) {
        mod <- survival::survreg( target ~.,  data = as.data.frame(dataset[, -i]), weights = wei, dist = "gaussian" )
        stat[i] <- 2 * abs(logLik(mod) - logLik(ini) )
        dof[i] <- dofini - length( coef(mod) ) 
      }

    mat <- cbind(1:p, pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
    colnames(mat) <- c("variable", "log.p-values", "statistic" )
    rownames(mat) <- 1:p 
    
    sel <- which.max( mat[, 2] )
    info <- matrix( c(0, -10, -10) , ncol = 3 )

    if ( mat[sel, 2] < threshold ) {
      res <- list(mat = mat, final = ini  ) 
      
    } else {
      
      info[1, ] <- mat[sel, ]
      mat <- mat[-sel, , drop = FALSE] 
      dat <- as.data.frame( dataset[, -sel] ) 
    } 
    
    i <- 1  
    
    if ( info[1, 2] > threshold ) {
      
      if ( !heavy ) {
        
        while ( info[i, 2] > threshold  &  NCOL(dat) > 0 )  {   
          i <- i + 1
          k <- p - i + 1
          ini <- survival::survreg( target ~., data = data.frame(dat), weights = wei, dist = "gaussian" )
          
          if ( k == 1 ) {
            mod <- survival::survreg(target ~ 1, data = data.frame(dat), weights = wei, dist = "gaussian")
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
              mat[, 2:3] <- c(stat, pval)
            }
            
          } else {
            
            stat <- dof <- numeric(k)
            
            for (j in 1:k) {
              mod <- survival::survreg( target ~.,  data = as.data.frame(dat[, -j]), weights = wei, dist = "gaussian" )
              stat[j] <- 2 * abs( logLik(mod) - logLik(ini) )
              dof[j] <- length( coef(ini) ) - length( coef(mod) ) 
            }
            mat[, 2:3] <- cbind( pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
            sel <- which.max( mat[, 2] )
            
            if ( mat[sel, 2] < threshold ) {
              final <- ini
              info <- rbind(info, c(0, -10, -10) )
              
            } else {
              info <- rbind(info, mat[sel, ] )
              mat <- mat[-sel, ,drop = FALSE] 
              dat <- as.data.frame( dataset[, -info[, 1] ] )
            }
          }  
        }  ## end while
        
      } else {
        
        while ( info[i, 2] > threshold  &  NCOL(dat) > 0 )  {   
          
          i <- i + 1       
          k <- p - i + 1
          ini <- survival::survreg( target ~., data = data.frame(dat), weights = wei )
          dofini <- length( coef(ini) ) 
          if ( k == 1 ) {
            mod <- survival::survreg(target ~ 1, data = data.frame(dat), weights = wei)
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
              mat[, 2:3] <- c(stat, pval)
            }  
          } else {		       
            
            stat <- dof <- numeric(k)
            
            for (j in 1:k) {
              mod <- survival::survreg( target ~., data = data.frame(dat[, -j]), weights = wei )
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
              dat <- as.data.frame( dataset[, -info[, 1] ] )
            }
            
          }	   
        }  ## end while
      }  ## end else 
    }  else final <- ini
    
  }  
  
  info <- info[ info[, 1] > 0, , drop = FALSE ]
  res <- list(runtime = runtime, info = info, mat = mat, ci_test = "testIndTobit", final = final ) 
  
}    






