gammabsreg <- function(target, dataset, threshold = 0.05, wei = NULL, heavy = FALSE, robust = FALSE) {
  
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
    
       if ( !heavy ) {
          ini <- glm( target ~.,  data = as.data.frame(dataset), family = Gamma(link = log), weights = wei, y = FALSE, model = FALSE )
          tab <- drop1( ini, test = "Chisq" )
          dof <- tab[-1, 1]
          stat <- tab[-1, 4]
          
        }  else {
          ini <- speedglm::speedglm( target ~.,  data = as.data.frame(dataset), family = Gamma(link = log), weights = wei )
          dofini <- length( coef(ini) )
          stat <- dof <- numeric(p)
          for (i in 1:p) {
            mod <- speedglm::speedglm( target ~.,  data = as.data.frame(dataset[, -i]), family = Gamma(link = log), weights = wei )
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
          mat <- mat[-sel, , drop = FALSE] 
          dat <- as.data.frame( dataset[, -sel] ) 
        } 
        
        i <- 1  
        
        if ( info[1, 2] > threshold ) {
          
          if ( !heavy ) {
            
            while ( info[i, 2] > threshold  &  NCOL(dat) > 0 )  {   
              i <- i + 1
              k <- p - i + 1
              ini <- glm( target ~., data = data.frame(dat), family = Gamma(link = log), weights = wei, y = FALSE, model = FALSE )
              
              if ( k == 1 ) {
                mod <- glm(target ~ 1, data = data.frame(dat), family = Gamma(link = log), weights = wei)
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
                  dat <- as.data.frame( dataset[, -info[, 1] ] )
                }
              }  
            }  ## end while
            
          } else {
            
            while ( info[i, 2] > threshold  &  NCOL(dat) > 0 )  {   
              
              i <- i + 1       
              k <- p - i + 1
              ini <- speedglm::speedglm( target ~., data = data.frame(dat), family = Gamma(link = log), weights = wei )
              dofini <- length( coef(ini) ) 
              if ( k == 1 ) {
                mod <- speedglm::speedglm(target ~ 1, data = data.frame(dat), family = Gamma(link = log), weights = wei)
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
                  mod <- speedglm::speedglm( target ~., data = data.frame(dat[, -j]), family = Gamma(link = log), weights = wei )
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
        
        info <- info[ info[, 1] > 0, , drop = FALSE ]
        res <- list(runtime = runtime, info = info, mat = mat, ci_test = "testIndGamma", final = final )       
    }  
  res

}    






