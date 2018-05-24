fbed.reg <- function(target, dataset, ini = NULL, test = NULL, threshold = 0.05, wei = NULL, K = 0, method = "LR", gam = NULL, backward = TRUE) {
  
  #check for NA values in the dataset and replace them with the variable median or the mode
  if ( any( is.na(dataset) ) ) {
    #dataset = as.matrix(dataset);
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    if ( is.matrix(dataset) )  {
      dataset <- apply( dataset, 2, function(za){ za[which(is.na(za))] = median(za, na.rm = TRUE) ; return(za) } ) 
    } else {
       poia <- unique( which( is.na(dataset), arr.ind = TRUE )[, 2] )
      for ( i in poia )  {
        xi <- dataset[, i]
        if ( is.numeric(xi) ) {                    
          xi[ which( is.na(xi) ) ] <- median(xi, na.rm = TRUE) 
        } else if ( is.factor( xi ) ) {
          xi[ which( is.na(xi) ) ] <- levels(xi)[ which.max( as.vector( table(xi) ) )]
        }
        dataset[, i] <- xi
      }
    }
  }
  
  id <- Rfast::check_data(dataset)
  if ( sum(id>0) > 0 )  dataset[, id] <- rnorm( dim(dataset)[1] * length(id) )
  dataset <- as.data.frame(dataset)
  
  if ( length(K) > 1 ) {
    
    result <- kfbed.reg(y = target, x = dataset, univ = ini, test = test, alpha = threshold, wei = NULL, K = K, method = method, gam = gam, backward = backward)
  
  } else {
   
  if (test == "gSquare") {
    
    result <- fbed.g2(y = target, x = dataset, alpha = threshold, univ = ini, K = K, backward = backward)
    
  } else {
    
  if (method == "LR") {
    
    if (test == "testIndReg") {
      result <- fbed.lm(target, dataset, alpha = threshold, univ = ini, wei = wei, K = K)
      
    } else if (test == "testIndFisher") {
      result <- Rfast::cor.fbed(target, as.matrix(dataset), alpha = threshold, K = K)
      result$univ <- NULL
      
    } else if (test == "testIndPois") {
      result <- fbed.glm(target, dataset, alpha = threshold, univ = ini, wei = wei, K = K, type = "poisson")
      
    } else if (test == "testIndNB") {
      result <- fbed.nb(target, dataset, alpha = threshold, univ = ini, wei = wei, K = K)
    
    } else if (test == "testIndLogistic") {
      result <- fbed.glm(target, dataset,alpha = threshold, univ = ini, wei = wei, K = K, type = "logistic")
      
    } else if ( test == "testIndMultinom" ) {
      result <- fbed.multinom(target, dataset, alpha = threshold, univ = ini, wei = wei, K = K)
      
    } else  if (test == "testIndOrdinal") {
      result <- fbed.ordinal(target, dataset, alpha = threshold, univ = ini, wei = wei, K = K)
    
    } else if (test == "testIndMMReg") {
      result <- fbed.mmreg(target, dataset, alpha = threshold, univ = ini, wei = wei, K = K)
      
    } else if (test == "testIndBinom") {
      result <- fbed.glm(target, dataset, alpha = threshold, univ = ini, wei = wei, K = K, type = "binomial")
      
    } else if (test == "censIndCR") {
      result <- fbed.cr(target, dataset, alpha = threshold, univ = ini, wei = wei, K = K)
      
    } else if (test == "censIndWR") {
      result <- fbed.wr(target, dataset, alpha = threshold, univ = ini, wei = wei, K = K)

    } else if (test == "testIndBeta") {
      result <- fbed.beta(target, dataset, alpha = threshold, univ = ini, wei = wei, K = K)
      
    } else if (test == "testIndZIP") {
      result <- fbed.zip(target, dataset, alpha = threshold, univ = ini, wei = wei, K = K)
      
    } else if (test == "testIndGamma") {
      result <- fbed.glm2(target, dataset, alpha = threshold, univ = ini, wei = wei, K = K, type = "gamma")
      
    } else if (test == "testIndNormLog") {
      result <- fbed.glm2(target, dataset, alpha = threshold, univ = ini, wei = wei, K = K, type = "normlog")
      
    } else if (test == "testIndTobit") {
      result <- fbed.tobit(target, dataset, alpha = threshold, univ = ini, wei = wei, K = K)
      
    } else if (test == "testIndClogit") {
      result <- fbed.clogit(target, dataset, alpha = threshold, univ = ini, wei = wei, K = K)
      
    } else if (test == "testIndQBinom") {
      result <- fbed.glm2(target, dataset, alpha = threshold, univ = ini, wei = wei, K = K, type = "quasibinomial")
      
    } else if (test == "testIndQPois") {
      result <- fbed.glm2(target, dataset, alpha = threshold, univ = ini, wei = wei, K = K, type = "quasipoisson")
    }
    
    result$back.rem <- 0
    result$back.n.tests <- 0
    
    if ( backward ) {
      
      if (result$info[1, 1] > 0) {
        a <- bs.reg(target, dataset[, result$res[, 1], drop = FALSE], threshold = threshold, wei = wei, test = test)
        
        if ( typeof(a) == "list" ) {
          result$back.rem <- result$res[a$info[, 1], 1]
          back.n.tests <- sum( dim(result$res)[1] : dim(a$mat)[1] )
          sel <- result$res[a$mat[, 1], 1] 
          stat <- a$mat[, 3]
          pval <- a$mat[, 2]
          result$res <- cbind(sel, stat, pval)
          result$back.n.tests <- back.n.tests
          result$runtime <- result$runtime + a$runtime
        } else {
          back.rem <- 0
          back.n.tests <- 0
          result$back.rem <- back.rem
          result$back.n.tests <- back.n.tests
          result$runtime <- result$runtime 
        }  
      }   ## end if (result$info[1, 1] > 0)
    }  ## end if ( backward )
    
  } else {  ## end of method =="LR"
    
    if (test == "testIndReg"  |  test == "testIndFisher") {
      test <- "testIndReg"
      result <- ebic.fbed.lm(target, dataset, univ = ini, gam = gam, wei = wei, K = K) 
      
    } else if (test == "testIndPois") {
      result <- ebic.fbed.glm(target, dataset, univ = ini, gam = gam, wei = wei, K = K, type = "poisson")
      
    } else if (test == "testIndNB") {
      result <- ebic.fbed.nb(target, dataset, univ = ini, gam = gam, wei = wei, K = K)
      
    } else if (test == "testIndLogistic") {
      result <- ebic.fbed.glm(target, dataset, univ = ini, gam = gam, wei = wei, K = K, type = "logistic")
      
    } else if (test == "testIndMultinom") {  
      result <- ebic.fbed.multinom(target, dataset, univ = ini, gam = gam, wei = wei, K = K)
        
    } else if (test == "testIndOrdinal") {
      result <- ebic.fbed.ordinal(target, dataset, univ = ini, gam = gam, wei = wei, K = K)
      
    } else if (test == "testIndMMReg") {
      result <- ebic.fbed.mmreg(target, dataset, univ = ini, gam = gam, wei = wei, K = K)
      
    } else if (test == "testIndBinom") {
      result <- ebic.fbed.glm(target, dataset, univ = ini, gam = gam, wei = wei, K = K, type = "binomial")
      
    } else if (test == "censIndCR") {
      result <- ebic.fbed.cr(target, dataset, univ = ini, gam = gam, wei = wei, K = K)
      
    } else if (test == "censIndWR") {
      result <- ebic.fbed.wr(target, dataset, univ = ini, gam = gam, wei = wei, K = K)

    } else if (test == "testIndBeta") {
      result <- ebic.fbed.beta(target, dataset, univ = ini, gam = gam, wei = wei, K = K)
      
    } else if (test == "testIndZIP") {
      result <- ebic.fbed.zip(target, dataset, univ = ini, gam = gam, wei = wei, K = K)
      
    } else if (test == "testIndGamma") {
      result <- ebic.fbed.glm(target, dataset, univ = ini, gam = gam, wei = wei, K = K, type = "gamma")
      
    } else if (test == "testIndNormLog") {
      result <- ebic.fbed.glm(target, dataset, univ = ini, gam = gam, wei = wei, K = K, type = "gaussian")
      
    } else if (test == "testIndTobit") {
      result <- ebic.fbed.tobit(target, dataset, univ = ini, gam = gam, wei = wei, K = K)
      
    } else if (test == "testIndClogit") {
      result <- ebic.fbed.clogit(target, dataset, univ = ini, gam = gam, wei = wei, K = K)
    }
    
    result$back.rem <- 0
    result$back.n.tests <- 0
    
    if ( backward ) {
      
      if ( result$info[1, 1] > 0 ) {
        a <- ebic.bsreg(target, dataset[, result$res[, 1], drop = FALSE], test = test, wei = wei, gam = gam) 
       
        if ( typeof(a) == "list" ) {
          back.n.tests <- sum( dim(result$res)[1] : length(a$mat[, 1]) )
          result$back.rem <- result$res[a$info[, 1], 1]
          sel <- result$res[ a$mat[, 1], 1]
          val <- a$mat[, 2]
          result$res <- cbind(sel, val)
          colnames(result$res) <- c("Vars", "eBIC difference")
          result$back.n.tests <- back.n.tests
          result$runtime <- result$runtime + a$runtime
        } else {
          back.rem <- 0
          back.n.tests <- 0
          result$back.rem <- back.rem
          result$back.n.tests <- back.n.tests
          result$runtime <- result$runtime 
        }  
      }   ## end if (result$info[1, 1] > 0) 
      
    }  ## end if ( backward )
    
  }  ## end of method == "eBIC"
  
  }  ## end if (test == "gSquare")
  
  }  ## end if ( length(K) > 1 )
  result
}