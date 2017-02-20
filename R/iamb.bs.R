iamb.bs <- function(target, dataset, threshold = 0.05, wei = NULL, test = NULL, user_test = NULL, robust = FALSE) {
  
  threshold <- log(threshold)
  dm <- dim(dataset)
  n <- dm[1]  ## sample size 
  p <- dm[2]  ## number of variables
  
  if ( p > n ) {
    res <- paste("The number of variables is hiher than the sample size. No backward procedure was attempted")
    
  } else {
    
    runtime <- proc.time()
    
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
    
    la <- length( unique(target) )

    if ( is.null(test)  &  is.null(user_test) ) {
      
      ## surival data
      if ( sum( class(target) == "Surv" ) == 1 ) {
        ci_test <- test <- "censIndCR"
        
        ## ordinal, multinomial or perhaps binary data
      } else if ( is.factor(target) ||  is.ordered(target) || la== 2 ) {
        ci_test <- test <- "testIndLogistic"
        
        ## count data
      } else if ( la > 2  &  !is.factor(target) ) {
        if ( sum( round(target) - target ) == 0 ) {
          ci_test <- test <- "testIndPois"
        } else  ci_test <- test <- "testIndReg"  
      }
    }
    
    av_models <- c("testIndReg", "testIndBeta", "censIndCR", "testIndRQ", "censIndWR", "testIndLogistic", "testIndPois", "testIndNB", "testIndZIP", "testIndSpeedglm");
    
    ci_test <- test
    test <- match.arg(test, av_models, TRUE);
    
    ############ 
    ###  GLMs 
    ############
    
    heavy <- FALSE
    if (test == "testIndSpeedglm")   heavy <- TRUE
    
    if ( test == "testIndPois"  ||  test == "testIndReg"  ||  ( test == "testIndLogistic"  &  la == 2 )  ||  test == "testIndSpeedglm" ) {
      
      res <- iamb.glmbs(target = target, dataset = dataset, wei = wei, threshold = exp(threshold), heavy = heavy, robust = robust ) 
      
    } else if ( test == "testIndBeta" ) {
      
      res <- iamb.betabs(target = target, dataset = dataset, wei = wei, threshold = exp(threshold) ) 
      
    } else if ( test == "testIndZIP" ) {
      
      res <- iamb.zipbs(target = target, dataset = dataset, wei = wei, threshold = exp(threshold) ) 	
      
    } else {
      
      if ( test == "censIndCR" ) {
        test <- survival::coxph 
        robust <- FALSE
        
      } else if ( test == "censIndWR" ) {
        test <- survival::survreg 
        robust <- FALSE
        
      } else if ( test == "testIndLogistic" ) {
        
        if ( is.ordered(target) )  {
          test <- ordinal::clm
          robust <- FALSE
          
        } else {
          test <- nnet::multinom
          robust <- FALSE
        }
        
      } else if ( test == "testIndNB" ) {
        test <- MASS::glm.nb
        robust <- FALSE
        
      } else if ( test == "testIndRQ" ) {
        test <- quantreg::rq
        robust <- FALSE

      } else if ( !is.null(user_test) ) {
        test <- user_test
      }
      
      ###################
      ###################
      
      a1 <- internaliamb.bs( target = target, dataset = dataset, threshold = threshold, test = test, wei = wei, p = p, robust = robust ) 
      ind <- 1:p
      a2 <- list()
      poies <- a1$mat[, 1]
      if ( length(poies) > 0 ) {
        ind[-poies] <- 0
        ind <- ind[ind > 0]
        dat <- dataset[, poies ]
        a2 <- internaliamb.bs(target = target, dataset = dat, threshold = threshold, test = test, wei = wei, p = length(ind), robust = robust )  
        poies <- a2$mat[, 1]
        ind[-poies] <- 0
        ind <- ind[ind > 0]
        if ( length(poies) == 1 )   dat <- dat  else   dat <- dat[, poies]
        i <- 2
      } else {
        ind <- NULL
        a2$mat <- NULL  
        a2$final <- paste("No variables were selected")
      } 
      while ( length(a1$mat[, 1]) - length(a2$mat[, 1]) != 0 ) {
        i <- i + 1
        a1 <- a2
        a2 <- internaliamb.bs( target = target, dataset = dat, threshold = threshold, test = test, wei = wei, p = length(ind), robust = robust ) 
        poies <- a2$mat[, 1]
        if ( length(poies) > 0 ) {
          ind[-poies] <- 0
          ind <- ind[ind > 0]
          if ( length(poies) == 1 )   dat <- dat  else   dat <- dat[, poies]
        } else  {
          dat <- NULL  
          ind <- NULL
        }  
      }
      
      runtime <- proc.time() - runtime

    # if ( !is.null(a2$mat) ) {
    #   final <- test( target ~.,  data = as.data.frame(dataset[, ind]), weights = wei )
    #   a2$mat[, 1] <- ind 
    # } else final <- "No variables were selected"
    res <- list(runtime = runtime, ci_test = ci_test, vars = ind, mat = a2$mat, final = a2$final ) 
  
    }
  }  
  res
}
      


      
      
         
internaliamb.bs <- function(target, dataset, threshold, test, wei, p, robust) {
  if ( !is.null(dataset) |  p > 0 ) {
    ini <- test( target ~.,  data = as.data.frame(dataset), weights = wei )
    dofini <- length( coef(ini) )
    likini <- logLik(ini) 
    stat <- dof <- numeric(p)
    
    if (p > 1) {
      for (j in 1:p) {
        mod <- test( target ~.,  data = as.data.frame(dataset[, -j]), weights = wei )
        stat[j] <- 2 * abs( likini - logLik(mod) )
        dof[j] <- dofini - length( coef(mod) ) 
      }
    } else {
      mod <- test( target ~ 1,  data = as.data.frame(dataset), weights = wei )
      stat <- 2 * abs( likini - logLik(mod) )
      dof <- dofini - length( coef(mod) ) 
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
        mod1 <- test(target ~ ., data = data.frame(dat), weights= wei )
        mod0 <- test(target ~ 1, weights = wei)
        stat <- 2 * abs( logLik(mod1) - logLik(mod0) )
        dof <- length( coef(mod1) ) - length( coef(mod0) ) 
        pval <- pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE)
        if (pval > threshold ) {
          final <- "No variables were selected"
          mat <- NULL
        } else final <- mod1
      } else  final <- test(target ~ ., data = data.frame(dat), weights= wei )
    }
    info <- info[ info[, 1] > 0, ]
    
  } else { 
    info <- NULL  
    mat <- NULL 
    final <- "No variables were selected"
  } 
  
  list(info = info, mat = mat, final = final) 
}  
