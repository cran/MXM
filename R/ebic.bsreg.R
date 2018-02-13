ebic.bsreg <- function(target, dataset, test = NULL, wei = NULL, gam = NULL) {
  
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
  
  if (test == "testIndReg") {
    result <- ebic.lm.bsreg(target, dataset, gam = gam, wei = wei)
  
  } else if (test == "testIndPois") {
    result <- ebic.glm.bsreg(target, dataset, gam = gam, wei = wei, type = "poisson")
  
  } else if (test == "testIndNB") {
    result <- ebic.nb.bsreg(target, dataset, gam = gam, wei = wei)
  
  } else if (test == "testIndLogistic") {
    
    if ( length(unique(target) ) == 2 ) {
      result <- ebic.glm.bsreg(target, dataset, wei = wei, gam = gam, type = "logistic")
    } else if ( length(unique(target) ) > 2 &  !is.ordered(target) ) {
      result <- ebic.multinom.bsreg(target, dataset, gam = gam, wei = wei)  
    } else  result <- ebic.ordinal.bsreg(target, dataset, gam = gam, wei = wei)
  
  } else if (test == "testIndMMReg") {
    result <- ebic.mm.bsreg(target, dataset, gam = gam, wei = wei)
    
  } else if (test == "testIndBinom") {
    result <- ebic.glm.bsreg(target, dataset, gam = gam, wei = wei, type = "binomial")
  
  } else if (test == "censIndCR") {
    result <- ebic.cr.bsreg(target, dataset, gam = gam, wei = wei)
  
  } else if (test == "censIndWR") {
    result <- ebic.wr.bsreg(target, dataset, gam = gam, wei = wei)
  
  } else if (test == "testIndBeta") {
    result <- ebic.beta.bsreg(target, dataset, gam = gam, wei = wei)
  
  } else if (test == "testIndZIP") {
    result <- ebic.zip.bsreg(target, dataset, gam = gam, wei = wei)
  
  } else if (test == "testIndGamma") {
    result <- ebic.glm.bsreg(target, dataset, gam = gam, wei = wei, type = "gamma")
  
  } else if (test == "testIndNormLog") {
    result <- ebic.glm.bsreg(target, dataset, gam = gam, wei = wei, type = "gaussian")
  
  } else if (test == "testIndTobit") {
    result <- ebic.tobit.bsreg(target, dataset, gam = gam, wei = wei)
    
  } else if (test == "testIndClogit") {
    result <- ebic.clogit.bsreg(target, dataset, gam = gam, wei = wei)
  }
  
  back.rem <- result$info[, 1]
  back.n.tests <- dim(dataset)[2]:dim(result$mat)[1] 
  result$mat <- result$mat
  result$back.rem <- back.rem
  result$back.n.tests <- sum(back.n.tests)
  result$runtime <- result$runtime 
  result
}   
