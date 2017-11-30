normlog.bsreg <- function(target, dataset, threshold = 0.05, wei = NULL, heavy = FALSE, robust = FALSE) {
  
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
    runtime <- proc.time()
    dataset <- as.data.frame(dataset)
    
    if ( !heavy ) {
      ini <- glm( target ~.,  data = dataset, family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE )
      dofini <- length( ini$coefficients )
      tab <- drop1( ini, test = "Chisq" )
      dof <- tab[-1, 1]
      stat <- tab[-1, 4]
      
    }  else {
      ini <- speedglm::speedglm( target ~.,  data = dataset, family = gaussian(link = log), weights = wei )
      dofini <- length( ini$coefficients )
      stat <- dof <- phi <- numeric(p)
      if (p == 1) {
        mod <- speedglm::speedglm( target ~ 1,  data = dataset, family = gaussian(link = log), weights = wei )
        phi <- mod$RSS/mod$df
        stat <- 2 * ( logLik(ini) - logLik(mod) ) / phi
        dof <- length( mod$coefficients ) 		  
      } else {
        for (i in 1:p) {
          mod <- speedglm::speedglm( target ~.,  data = dataset[, -i, drop = FALSE], family = gaussian(link = log), weights = wei )
          stat[i] <- 2 * ( logLik(ini) - logLik(mod) ) / ( mod$RSS/mod$df )
          dof[i] <- length( mod$coefficients ) 
        }
      }
    }	
    
    mat <- cbind(1:p, pf( stat, dofini - dof, n - dof, lower.tail = FALSE, log.p = TRUE), stat )
    colnames(mat) <- c("variable", "log.p-values", "statistic" )
    rownames(mat) <- 1:p 
    
    sel <- which.max( mat[, 2] )
    info <- matrix( c(0, -10, -10) , ncol = 3 )
    sela <- sel 
    
    if ( mat[sel, 2] < threshold ) {
      runtime <- proc.time() - runtime
      res <- list(runtime = runtime, info = matrix(0, 0, 3), mat = mat, ci_test = "testIndNormLog", final = ini ) 
      
    } else {
      
      info[1, ] <- mat[sel, ]
      mat <- mat[-sel, , drop = FALSE] 
      dat <- dataset[, -sel, drop = FALSE] 
      
      i <- 1  
      
      if ( info[1, 2] > threshold & dim(mat)[1]>0 ) {
        
        if ( !heavy ) {
          
          while ( info[i, 2] > threshold  &  dim(dat)[2] > 0 )  {   
            i <- i + 1
            k <- p - i + 1
            ini <- glm( target ~., data = dat, family = gaussian(link = log), weights = wei, y = FALSE, model = FALSE )
            dofini <- length(ini$coefficients)
            if ( k == 1 ) {
              mod <- glm(target ~ 1, data = dat, family = gaussian(link = log), weights = wei)
              stat <- 2 * ( logLik(ini) - logLik(mod) ) / summary(mod)[[ 14 ]]
              dof <- length( mod$coefficients ) 
              pval <- pf( stat, dofini - dof, n - dof, lower.tail = FALSE, log.p = TRUE)
              
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
              mat[, 2:3] <- cbind( pf(stat, dof, n - (dofini - dof), lower.tail = FALSE, log.p = TRUE), stat )
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
          }  ## end while
          
        } else {
          
          while ( info[i, 2] > threshold  &  dim(dat)[2] > 0 )  {   
            
            i <- i + 1       
            k <- p - i + 1
            ini <- speedglm::speedglm( target ~., data = dat, family = gaussian(link = log), weights = wei )
            dofini <- length( ini$coefficients ) 
            if ( k == 1 ) {
              mod <- speedglm::speedglm(target ~ 1, data = dat, family = gaussian(link = log), weights = wei)
              stat <- 2 * ( logLik(ini) - logLik(mod) ) /  (mod$RSS/mod$df)
              dof <- length( mod$coefficients ) 
              pval <- pf( stat, dofini - dof, n - dof, lower.tail = FALSE, log.p = TRUE)
              
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
                mod <- speedglm::speedglm( target ~., data = dat[, -j, drop = FALSE], family = gaussian(link = log), weights = wei )
                stat[j] <- 2 * ( logLik(ini) - logLik(mod) ) / (mod$RSS/mod$df)
                dof[j] <- length( mod$coefficients ) 		   
              } 
              
              mat[, 2:3] <- cbind( pf(stat, dofini - dof, n - dof, lower.tail = FALSE, log.p = TRUE), stat )
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
          }  ## end while
        }  ## end else 
        runtime <- proc.time() - runtime		
        info <- info[ info[, 1] > 0, , drop = FALSE]
        res <- list(runtime = runtime, info = info, mat = mat, ci_test = "testIndNormLog", final = final ) 
        
      } else {
        runtime <- proc.time() - runtime
        res <- list(runtime = runtime, info = info, mat = NULL, ci_test = "testIndNormLog", final = mod ) 
      }
    }  
  }  
  res
  
}    






