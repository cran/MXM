######## Internal for linear regression
internaliamb.lmbs <- function(target, dataset, threshold, wei, p, heavy = FALSE, robust = FALSE) {
  if ( !is.null(dataset) |  p > 0 ) {
    n <- length(target)
    if ( p > 1 ) {
      if ( !heavy ) {
        if ( !robust ) {
          ini <- lm( target ~., data = data.frame(dataset), weights = wei, y = FALSE, model = FALSE )
          df2 <- n - length( coef(ini) )
          tab <- drop1( ini, test = "F" )
          df1 <- tab[-1, 1]
          stat <- tab[-1, 5]
          mat <- cbind(1:p, pf( stat, df1, df2, lower.tail = FALSE, log.p = TRUE), stat )
        } else {
          ini <- MASS::rlm( target ~., data = as.data.frame(dataset), weights = wei, maxit = 2000, method = "MM" )
          lik1 <- as.numeric( logLik(ini) )
          dofini <- length( coef(ini) ) 
          for (i in 1:p) {
            fit2 <- MASS::rlm( target ~., data = as.data.frame(dataset[, -i]), weights = wei, maxit = 2000, method = "MM" )
            stat[i] <- 2 * abs( lik1 - as.numeric( logLik(fit2) ) ) 
            dof[i] <- dofini - length( coef(fit2) )
          }
          mat <- cbind(1:p, pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
        }  
        
      } else {
        ini <- speedglm::speedlm( target ~.,  data = data.frame(dataset), weights = wei )
        dofini <- length( coef(ini) )
        stat <- pval <- numeric(p)   ;   d1 <- length( coef(ini) )
        for (i in 1:p) {
          fit2 = speedglm::speedlm( target ~., data = as.data.frame(dataset[, -i]), weights = wei )
          df1 = d1 - length( coef(fit2) )
          df2 = n - d1
          stat[i] = (fit2$RSS - ini$RSS)/df1 / ( ini$RSS /df2 )
          pval[i] = pf( stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE )
        }
        mat <- cbind(1:p, pval, stat)
      }	
    } else {
      if ( !heavy ) {
        if ( !robust ) {
          ini <- lm( target ~., data = data.frame(dataset), weights = wei, y = FALSE, model = FALSE )
          df2 <- n - length( coef(ini) )
          tab <- drop1( ini, test = "F" )
          df1 <- tab[-1, 1]
          stat <- tab[-1, 5]
          mat <- cbind(1, pf( stat, df1, df2, lower.tail = FALSE, log.p = TRUE), stat )
        } else  {
          ini <- MASS::rlm( target ~., data = as.data.frame(dataset), weights = wei, maxit = 2000, method = "MM" )
          mod0 <- MASS::rlm(target ~ 1, weights = wei, maxit = 2000, method = "MM" ) 
          stat <- 2 * abs( as.numeric( logLik(ini) ) - as.numeric( logLik(mod0) ) )
          mat <- cbind(1, pchisq( stat, length( coef(ini) ) - 1, lower.tail = FALSE, log.p = TRUE), stat )
        }  
      } else {
        ini <- speedglm::speedlm( target ~.,  data = data.frame(dataset), weights = wei )
        fit2 <- speedglm::speedlm( target ~ 1, data = as.data.frame(dataset), weights = wei )
        df1 <- length( coef(ini) ) - length( coef(fit2) )
        df2 <- n - length( coef(ini) )
        stat <- (fit2$RSS - ini$RSS)/df1 / ( ini$RSS /df2 )
        mat <- cbind(1, pf( stat, df1, df2, lower.tail = FALSE, log.p = TRUE), stat )
      }  
    }
    
    colnames(mat) <- c("variable", "log.p-values", "statistic" )
    rownames(mat) <- 1:p  
    
    info <- matrix( c(0, -10, -10) , ncol = 3 )
    sel <- which( mat[, 2] > threshold ) 
    
    if ( length(sel) == 0 ) {
      final <- ini 
      
    } else {
      
      info <- mat[sel, , drop = FALSE]
      mat <- mat[-sel, , drop = FALSE] 
      dat <- as.data.frame( dataset[, -sel] ) 
      
      if ( p - length(sel) == 0 ) {
        final <- "No variables were selected"
        mat <- NULL
      } else if ( p - length(sel) == 1 ) {
        if ( !heavy ) {
          if ( !robust ) {
            ini <- lm( target ~., data = data.frame(dataset), weights = wei, y = FALSE, model = FALSE )
            df2 <- n - length( coef(ini) )
            tab <- drop1( ini, test = "F" )
            df1 <- tab[-1, 1]
            stat <- tab[-1, 5]
            pval <- pf(stat, df1, df2, lower.tail = FALSE, log.p = TRUE)
          } else {
            ini <- MASS::rlm( target ~., data = as.data.frame(dataset), weights = wei, maxit = 2000, method = "MM" )
            mod0 <- MASS::rlm( target ~ 1, weights = wei, maxit = 2000, method = "MM" ) 
            stat <- 2 * abs( as.numeric( logLik(ini) ) - as.numeric( logLik(mod0) ) )
            pval <- pchisq(stat, length( coef(ini) ) - 1, lower.tail = FALSE, log.p = TRUE)
          }
        } else  {
          mod1 <- speedglm::speedlm( target ~ .,  data = data.frame(dat), weights = wei )
          mod0 <- speedglm::speedlm( target ~ 1,  data = data.frame(dat), weights = wei )
          df1 <- length( coef(mod1) ) - length( coef(mod0) )
          df2 <- n - length( coef(mod1) )
          stat <- (fit2$RSS - ini$RSS)/df1 / ( ini$RSS /df2 )
          pval <- pf( stat, df1, df2, lower.tail = FALSE, log.p = TRUE)
        }  
        if (pval > threshold ) {
          final <- "No variables were selected"
          mat <- NULL
        } else final <- mod1
      } else {
        if ( !heavy ) {
          if ( !robust ) {
            final <- lm( target ~., data = data.frame(dataset), weights = wei, y = FALSE, model = FALSE )
          } else    final <- MASS::rlm( target ~., data = as.data.frame(dataset), weights = wei, maxit = 2000, method = "MM" )
        } else  final <- speedglm::speedlm( target ~ .,  data = data.frame(dat), weights = wei )
      }
    }
    info <- info[ info[, 1] > 0, ]
    
  } else { 
    info <- NULL  
    mat <- NULL 
    final <- "No variables were selected"
  } 
  
  list(info = info, mat = mat, final = final) 
} 