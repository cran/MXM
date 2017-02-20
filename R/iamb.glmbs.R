iamb.glmbs <- function(target, dataset, threshold = 0.05, wei = NULL, test = NULL, heavy = FALSE, robust = FALSE) {
  
  threshold <- log(threshold)
  dm <- dim(dataset)
  n <- dm[1]  ## sample size 
  p <- dm[2]  ## number of variables
  
  if ( p > n ) {
    res <- paste("The number of variables is hiher than the sample size. No backward procedure was attempted")
    
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
    
    if ( !is.null(test) ) {
      ci_test <- test 
      if ( test == "testIndBinom" ) {
        wei <- target[, 2]
        target <- target[, 1] / wei
        funa <- internaliamb.binombs
      } else if (test == "testIndLogistic" ) {
        funa <- internaliamb.binombs
      } else if ( test == "tesIndPois" ) {
        funa <- internaliamb.poisbs
      } else   funa <- internaliamb.lmbs
      
    } else {
      la <- length( unique(target) )
      if ( la == 2 ) {
        ci_test <- "testIndLogistic"
        funa <- internaliamb.binombs
      } else if ( is.matrix(target) ) {
        ci_test <- "testIndBinom"
      } else if (la > 2  &  sum( target - round(target) ) == 0 ) {
        ci_test <- "testIndPois"
        funa <- internaliamb.poisbs
      } else {
        ci_test <- "testIndReg"
        funa <- internaliamb.lmbs
      }  
    }
    
    a1 <- funa( target = target, dataset = dataset, threshold = threshold, wei = wei, p = p, heavy = heavy, robust = robust ) 
    ind <- 1:p
    a2 <- list()
    poies <- a1$mat[, 1]
    if ( length(poies) > 0 ) {
      ind[-poies] <- 0
      ind <- ind[ind > 0]
      dat <- dataset[, poies ]
      a2 <- funa(target = target, dataset = dat, threshold = threshold, wei = wei, p = length(ind), heavy = heavy, robust = robust )  
      poies <- a2$mat[, 1]
      ind[-poies] <- 0
      ind <- ind[ind > 0]
      if ( length(poies) == 1 )   dat <- dat  else   dat <- dat[, poies]
      i <- 2
    } else {
      ind <- NULL
      a2$mat <- NULL  
    }
    while ( length(a1$mat[, 1]) - length(a2$mat[, 1]) != 0 ) {
      i <- i + 1
      a1 <- a2
      a2 <- funa( target = target, dataset = dat, threshold = threshold, wei = wei, p = length(ind), heavy = heavy, robust = robust ) 
      poies <- a2$mat[, 1]
      if ( length(poies) > 0 ) {
        ind[-poies] <- 0
        ind <- ind[ind > 0]
        if ( length(poies) == 1 )   dat <- dat  else   dat <- dat[, poies]
      } else {
        ind <- NULL
        dat <- NULL  
      }  
    }
    
    res <- list(info = ind, mat = a2$mat, ci_test = ci_test, final = a2$final ) 
  }
  
  res
}  




######## Internal for logistic regression
internaliamb.binombs <- function(target, dataset, threshold, wei, p, heavy = FALSE, robust = FALSE) {
  if ( !is.null(dataset) |  p > 0 ) {
    if ( p > 1 ) {
      if ( !heavy ) {
        ini <- glm( target ~.,  data = data.frame(dataset), family = binomial(logit), weights = wei, y = FALSE, model = FALSE )
        tab <- drop1( ini, test = "Chisq" )
        dof <- tab[-1, 1]
        stat <- tab[-1, 4]
      
      } else {
        ci_test <- "testIndSpeedglm"
        ini <- speedglm::speedglm( target ~.,  data = data.frame(dataset), family = binomial(logit), weights = wei )
        dofini <- length( coef(ini) )
        stat <- dof <- numeric(p)
        for (i in 1:p) {
          mod <- speedglm::speedglm( target ~.,  data = data.frame(dataset[, -i]), family = binomial(logit), weights = wei )
          stat[i] <- mod$deviance - ini$deviance
          dof[i] <- dofini - length( coef(mod) ) 
        }
      }	
    } else {
      if ( !heavy ) {
        ini <- glm( target ~.,  data = data.frame(dataset), family = binomial(logit), weights = wei, y = FALSE, model = FALSE )
        mod0 <- glm(target ~ 1, family = binomial(logit), weights = wei, y = FALSE, model = FALSE)
      } else {
        ini <- speedglm::speedglm( target ~.,  data = data.frame(dataset), family = binomial(logit), weights = wei )
        mod0 <- speedglm::speedglm( target ~ 1,  data = data.frame(dataset), family = binomial(logit), weights = wei )
      } 
      stat <- mod0$deviance - ini$deviance
      dof <- length( coef(ini) ) - 1
    }
    
    mat <- cbind(1:p, pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
    colnames(mat) <- c("variable", "log.p-values", "statistic" )
    rownames(mat) <- 1:p  
    
    info <- matrix( c(0, -10, -10) , ncol = 3 )
    sel <- which( mat[, 2] > threshold ) 
    
    if ( length(sel) == 0 ) {
      final <- ini 
      
    } else {
      
      info <- mat[sel, ]
      mat <- mat[-sel, ] 
      if ( !is.matrix(info) )   info <- matrix(info, ncol = 3) 
      if ( !is.matrix(mat) )   mat <- matrix(mat, ncol = 3) 
      dat <- as.data.frame( dataset[, -sel] ) 
      
      if ( p - length(sel) == 0 ) {
        final <- "No variables were selected"
        mat <- NULL
      } else if ( p - length(sel) == 1 ) {
        if ( !heavy ) {
          mod1 <- glm(target ~., data = data.frame(dat), family = binomial(logit), weights = wei, y = FALSE, model = FALSE)
          mod0 <- glm(target ~ 1, family = binomial(logit), weights = wei, y = FALSE, model = FALSE)
        } else  {
          mod1 <- speedglm::speedglm( target ~ .,  data = data.frame(dat), family = binomial(logit), weights = wei )
          mod0 <- speedglm::speedglm( target ~ 1,  data = data.frame(dat), family = binomial(logit), weights = wei )
        }  
        stat <- abs( mod1$deviance - mod0$deviance )
        pval <- pchisq( stat, length( coef(mod1) ) - 1, lower.tail = FALSE, log.p = TRUE)
        if (pval > threshold ) {
          final <- "No variables were selected"
          mat <- NULL
        } else final <- mod1
      } else {
        if ( !heavy ) {
          final <- glm(target ~., data = data.frame(dat), family = poisson(log), weights = wei, y = FALSE, model = FALSE)
        } else  final <- speedglm::speedglm( target ~ .,  data = data.frame(dat), family = poisson(log), weights = wei )
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



######## Internal for Poisson regression
internaliamb.poisbs <- function(target, dataset, threshold, wei, p, heavy = FALSE, robust = FALSE) {
  if ( !is.null(dataset) |  p > 0 ) {
    if ( p > 1 ) {
      if ( !heavy ) {
        ini <- glm( target ~.,  data = data.frame(dataset), family = poisson(log), weights = wei, y = FALSE, model = FALSE )
        tab <- drop1( ini, test = "Chisq" )
        dof <- tab[-1, 1]
        stat <- tab[-1, 4]
        
      } else {
        ci_test <- "testIndSpeedglm"
        ini <- speedglm::speedglm( target ~.,  data = data.frame(dataset), family = poisson(log), weights = wei )
        dofini <- length( coef(ini) )
        stat <- dof <- numeric(p)
        for (i in 1:p) {
          mod <- speedglm::speedglm( target ~.,  data = data.frame(dataset[, -i]), family = poisson(log), weights = wei )
          stat[i] <- mod$deviance - ini$deviance
          dof[i] <- dofini - length( coef(mod) ) 
        }
      }	
    } else {
      if ( !heavy ) {
        ini <- glm( target ~.,  data = data.frame(dataset), family = poisson(log), weights = wei, y = FALSE, model = FALSE )
        mod0 <- glm(target ~ 1, family = poisson(log), weights = wei, y = FALSE, model = FALSE)
      } else {
        ini <- speedglm::speedglm( target ~.,  data = data.frame(dataset), family = poisson(log), weights = wei )
        mod0 <- speedglm::speedglm( target ~ 1,  data = data.frame(dataset), family = poisson(log), weights = wei )
      }  
      stat <- mod0$deviance - ini$deviance
      dof <- length( coef(ini) ) - 1
    }
    
    mat <- cbind(1:p, pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
    colnames(mat) <- c("variable", "log.p-values", "statistic" )
    rownames(mat) <- 1:p  
    
    info <- matrix( c(0, -10, -10) , ncol = 3 )
    sel <- which( mat[, 2] > threshold ) 
    
    if ( length(sel) == 0 ) {
      final <- ini 
      
    } else {
      
      info <- mat[sel, ]
      mat <- mat[-sel, ] 
      if ( !is.matrix(info) )   info <- matrix(info, ncol = 3) 
      if ( !is.matrix(mat) )   mat <- matrix(mat, ncol = 3) 
      dat <- as.data.frame( dataset[, -sel] ) 
      
      if ( p - length(sel) == 0 ) {
        final <- "No variables were selected"
        mat <- NULL
      } else if ( p - length(sel) == 1 ) {
        if ( !heavy ) {
          mod1 <- glm(target ~., data = data.frame(dat), family = poisson(log), weights = wei, y = FALSE, model = FALSE)
          mod0 <- glm(target ~ 1, family = poisson(log), weights = wei, y = FALSE, model = FALSE)
        } else  {
          mod1 <- speedglm::speedglm( target ~ .,  data = data.frame(dat), family = poisson(log), weights = wei )
          mod0 <- speedglm::speedglm( target ~ 1,  data = data.frame(dat), family = poisson(log), weights = wei )
        }  
        stat <- abs( mod1$deviance - mod0$deviance )
        pval <- pchisq( stat, length( coef(mod1) ) - 1, lower.tail = FALSE, log.p = TRUE)
        if (pval > threshold ) {
          final <- "No variables were selected"
          mat <- NULL
        } else final <- mod1
      } else {
        if ( !heavy ) {
          final <- glm(target ~., data = data.frame(dat), family = poisson(log), weights = wei, y = FALSE, model = FALSE)
        } else  final <- speedglm::speedglm( target ~ .,  data = data.frame(dat), family = poisson(log), weights = wei )
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
          ini <- MASS::rlm( target ~., data = as.data.frame(dataset), weights = wei, maxit = 2000)
          lik1 <- as.numeric( logLik(ini) )
          dofini <- length( coef(ini) ) 
          for (i in 1:p) {
            fit2 <- MASS::rlm( target ~., data = as.data.frame(dataset[, -i]), weights = wei, maxit = 2000 )
            stat[i] <- 2 * abs( lik1 - as.numeric( logLik(fit2) ) ) 
            dof[i] <- dofini - length( coef(fit2) )
          }
          mat <- cbind(1:p, pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
        }  

      } else {
        ci_test <- "testIndSpeedglm"
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
          ini <- MASS::rlm( target ~., data = as.data.frame(dataset), weights = wei, maxit = 2000)
          mod0 <- MASS::rlm(target ~ 1, weights = wei, maxit = 2000) 
          stat <- 2 * abs( lik1 - as.numeric( logLik(mod0) ) )
          mat <- cbind(1, pchisq( stat, length( coef(ini) ) - 1, lower.tail = FALSE, log.p = TRUE), stat )
        }  
      } else {
        ini <- speedglm::speedlm( target ~.,  data = data.frame(dataset), weights = wei )
        fit2 <- speedglm::speedlm( target ~1, data = as.data.frame(dataset), weights = wei )
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
      
      info <- mat[sel, ]
      mat <- mat[-sel, ] 
      if ( !is.matrix(info) )   info <- matrix(info, ncol = 3) 
      if ( !is.matrix(mat) )   mat <- matrix(mat, ncol = 3) 
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
            ini <- MASS::rlm( target ~., data = as.data.frame(dataset), weights = wei, maxit = 2000)
            mod0 <- MASS::rlm(target ~ 1, weights = wei, maxit = 2000) 
            stat <- 2 * abs( lik1 - as.numeric( logLik(mod0) ) )
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
          } else    final <- MASS::rlm( target ~., data = as.data.frame(dataset), weights = wei, maxit = 2000)
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

