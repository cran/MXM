ebic.bsreg <- function(y, x, test = NULL, wei = NULL, gam = 0) {
  
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
  
  if (test == "testIndReg") {
    result <- ebic.lm.bsreg(y, x, gam = gam, wei = wei)
  
  } else if (test == "testIndPois") {
    result <- ebic.glm.bsreg(y, x, gam = gam, wei = wei, type = "poisson")
  
  } else if (test == "testIndNB") {
    result <- ebic.nb.bsreg(y, x, gam = gam, wei = wei)
  
  } else if (test == "testIndLogistic") {
    
    if ( length(unique(y) == 2) ) {
      result <- ebic.glm.bsreg(y, x, gam = gam, wei = wei, type = "logistic")
    } else if ( length(unique(y) > 2) &  !is.ordered(y) ) {
      result <- ebic.multinom.bsreg(y, x, gam = gam, wei = wei)  
    } else  result <- ebic.ordinal.bsreg(y, x, gam = gam, wei = wei)
  
  } else if (test == "testIndMMreg") {
    result <- ebic.mm.bsreg(y, x, gam = gam, wei = wei)
    
  } else if (test == "testIndBinom") {
    result <- ebic.glm.bsreg(y, x, gam = gam, wei = wei, type = "binomial")
  
  } else if (test == "censIndCR") {
    result <- ebic.cr.bsreg(y, x, gam = gam, wei = wei)
  
  } else if (test == "censIndWR") {
    result <- ebic.wr.bsreg(y, x, gam = gam, wei = wei)
  
  } else if (test == "testIndBeta") {
    result <- ebic.beta.bsreg(y, x, gam = gam, wei = wei)
  
  } else if (test == "testIndZIP") {
    result <- ebic.zip.bsreg(y, x, gam = gam, wei = wei)
  
  } else if (test == "testIndGamma") {
    result <- ebic.glm.bsreg(y, x, gam = gam, wei = wei, type = "gamma")
  
  } else if (test == "testIndNormLog") {
    result <- ebic.glm.bsreg(y, x, gam = gam, wei = wei, type = "gaussian")
  
  } else if (test == "testIndTobit") {
    result <- ebic.tobit.bsreg(y, x, gam = gam, wei = wei)
  }

  back.rem <- result$info[, 1]
  back.n.tests <- sum( dim(result$mat)[1] : length(result$mat[, 1]) ) 
  result$res <- result$mat[, 1]
  result$back.rem <- back.rem
  result$back.n.tests. <- back.n.tests
  result$runtime <- result$runtime 
  result
}   
