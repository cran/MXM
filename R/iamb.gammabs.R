iamb.gammabs <- function(target, dataset, threshold = 0.05, wei = NULL, heavy = FALSE) {
  
  threshold <- log(threshold)
  dm <- dim(dataset)
  n <- dm[1]  ## sample size 
  p <- dm[2]  ## number of variables
  if ( p > n ) {
    res <- paste("The number of variables is hiher than the sample size. No backward procedure was attempted")
  } else {
    #check for NA values in the dataset and replace them with the variable median or the mode
    if ( any(is.na(dataset)) ) {
      #dataset = as.matrix(dataset);
      warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
      if ( is.matrix(dataset) )  {
        dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
      } else {
        poia <- unique( which( is.na(dataset), arr.ind = TRUE )[, 2] )
        for ( i in poia )  {
          xi <- dataset[, i]
          if ( is.numeric(xi) )  {                    
            xi[ which( is.na(xi) ) ] <- median(xi, na.rm = TRUE) 
          } else if ( is.factor( xi ) ) {
            xi[ which( is.na(xi) ) ] <- levels(xi)[ which.max( as.vector( table(xi) ) )]
          }
          dataset[, i] <- xi
        }
      }
    }

    a1 <- internaliamb.gammabs( target = target, dataset = dataset, threshold = threshold, wei = wei, p = p, heavy = heavy ) 
    ind <- 1:p
    a2 <- list()
    poies <- a1$mat[, 1]
    if ( length(poies) > 0 ) {
      ind[-poies] <- 0
      ind <- ind[ind > 0]
      dat <- dataset[, poies ]
      a2 <- internaliamb.gammabs( target = target, dataset = dat, threshold = threshold, wei = wei, p = length(ind), heavy = heavy ) 
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
      a2 <- internaliamb.gammabs( target = target, dataset = dat, threshold = threshold, wei = wei, p = length(ind), heavy = heavy ) 
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
    
    res <- list(info = ind, mat = a2$mat, ci_test = "testIndGamma", final = a2$final ) 
  }
  
  res
}  









