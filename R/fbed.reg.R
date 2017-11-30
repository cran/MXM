fbed.reg <- function(y, x, test = NULL, alpha = 0.05, wei = NULL, K = 0, method = "LR", gam = NULL, backward = TRUE) {
  
  #check for NA values in the dataset and replace them with the variable median or the mode
  if ( any( is.na(x) ) ) {
    #dataset = as.matrix(dataset);
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    if ( is.matrix(x) )  {
      x <- apply( x, 2, function(za){ za[which(is.na(za))] = median(za, na.rm = TRUE) ; return(za) } ) 
    } else {
       poia <- unique( which( is.na(x), arr.ind = TRUE )[, 2] )
      for ( i in poia )  {
        xi <- x[, i]
        if ( is.numeric(xi) ) {                    
          xi[ which( is.na(xi) ) ] <- median(xi, na.rm = TRUE) 
        } else if ( is.factor( xi ) ) {
          xi[ which( is.na(xi) ) ] <- levels(xi)[ which.max( as.vector( table(xi) ) )]
        }
        x[, i] <- xi
      }
    }
  }
  
  id <- Rfast::check_data(x)
  if ( sum(id>0) > 0 )  x[, id] <- rnorm( dim(x)[1] * length(id) )
  x <- as.data.frame(x)
  
  if (method == "LR") {
    
    if (test == "testIndReg") {
      result <- fbed.lm(y, x, alpha = alpha, wei = wei, K = K)
      
    } else if (test == "testIndFisher") {
      result <- Rfast::cor.fbed(y, as.matrix(x), alpha = alpha, K = K)
      
    } else if (test == "testIndPois") {
      result <- fbed.glm(y, x, alpha = alpha, wei = wei, K = K, type = "poisson")

    } else if (test == "testIndNB") {
      result <- fbed.nb(y, x, alpha = alpha, wei = wei, K = K)
    
    } else if (test == "testIndLogistic") {
      
      if ( length(unique(y) ) == 2 ) {
        result <- fbed.glm(y, x, alpha = alpha, wei = wei, K = K, type = "logistic")
      } else if ( !is.ordered(y) ) {
        result <- fbed.multinom(y, x, alpha = alpha, wei = wei, K = K)
      } else  result <- fbed.ordinal(y, x, alpha = alpha, wei = wei, K = K)
    
    } else if (test == "testIndMMreg") {
      result <- fbed.mmreg(y, x, alpha = alpha, wei = wei, K = K)
      
    } else if (test == "testIndBinom") {
      result <- fbed.glm(y, x, alpha = alpha, wei = wei, K = K, type = "binomial")
      
    } else if (test == "censIndCR") {
      result <- fbed.cr(y, x, alpha = alpha, wei = wei, K = K)
      
    } else if (test == "censIndWR") {
      result <- fbed.wr(y, x, alpha = alpha, wei = wei, K = K)

    } else if (test == "testIndBeta") {
      result <- fbed.beta(y, x, alpha = alpha, wei = wei, K = K)
      
    } else if (test == "testIndZIP") {
      result <- fbed.zip(y, x, alpha = alpha, wei = wei, K = K)
      
    } else if (test == "testIndGamma") {
      result <- fbed.glm(y, x, alpha = alpha, wei = wei, K = K, type = "gamma")
      
    } else if (test == "testIndNormLog") {
      result <- fbed.glm(y, x, alpha = alpha, wei = wei, K = K, type = "normlog")
      
    } else if (test == "testIndTobit") {
      result <- fbed.tobit(y, x, alpha = alpha, wei = wei, K = K)
      
    } else if (test == "testIndClogit") {
      result <- fbed.clogit(y, x, alpha = alpha, wei = wei, K = K)
    }
    
    result$back.rem <- 0
    result$back.n.tests <- 0
    
    if ( backward ) {
      
      if (test == "testIndMMreg") {
        test <- "testIndReg"
        robust <- TRUE
      } else robust =  FALSE
      
      if (result$info[1, 1] > 0) {
        a <- bs.reg(y, x[, result$res[, 1], drop = FALSE], threshold = alpha, wei = wei, test = test, robust = robust)
        
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
    
  } else {
    
    if (test == "testIndReg") {
      result <- ebic.fbed.lm(y, x, gam = gam, wei = wei, K = K) 
      
    } else if (test == "testIndPois") {
      result <- ebic.fbed.glm(y, x, gam = gam, wei = wei, K = K, type = "poisson")
      
    } else if (test == "testIndNB") {
      result <- ebic.fbed.nb(y, x, gam = gam, wei = wei, K = K)
      
    } else if (test == "testIndLogistic") {
      
      if ( length(unique(y) ) == 2 ) {
        result <- ebic.fbed.glm(y, x, gam = gam, wei = wei, K = K, type = "logistic")
      } else if ( !is.ordered(y) ) {
        result <- ebic.fbed.multinom(y, x, gam = gam, wei = wei, K = K)
      } else  result <- ebic.fbed.ordinal(y, x, gam = gam, wei = wei, K = K)
      
    } else if (test == "testIndMMreg") {
      result <- ebic.fbed.mmreg(y, x, gam = gam, wei = wei, K = K)
      
    } else if (test == "testIndBinom") {
      result <- ebic.fbed.glm(y, x, gam = gam, wei = wei, K = K, type = "binomial")
      
    } else if (test == "censIndCR") {
      result <- ebic.fbed.cr(y, x, gam = gam, wei = wei, K = K)
      
    } else if (test == "censIndWR") {
      result <- ebic.fbed.wr(y, x, gam = gam, wei = wei, K = K)

    } else if (test == "testIndBeta") {
      result <- ebic.fbed.beta(y, x, gam = gam, wei = wei, K = K)
      
    } else if (test == "testIndZIP") {
      result <- ebic.fbed.zip(y, x, gam = gam, wei = wei, K = K)
      
    } else if (test == "testIndGamma") {
      result <- ebic.fbed.glm(y, x, gam = gam, wei = wei, K = K, type = "gamma")
      
    } else if (test == "testIndNormLog") {
      result <- ebic.fbed.glm(y, x, gam = gam, wei = wei, K = K, type = "gaussian")
      
    } else if (test == "testIndTobit") {
      result <- ebic.fbed.tobit(y, x, gam = gam, wei = wei, K = K)
      
    } else if (test == "testIndClogit") {
      result <- ebic.fbed.clogit(y, x, gam = gam, wei = wei, K = K)
    }
    
    result$back.rem <- 0
    result$back.n.tests <- 0
    
    if ( backward ) {
      
      if (result$info[1, 1] > 0) {
        a <- ebic.bsreg(y, x[, result$res[, 1], drop = FALSE], test = test, wei = wei, gam = gam) 
       
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
    
  } 
  
  result
}