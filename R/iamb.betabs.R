iamb.betabs <- function(target, dataset, threshold = 0.05, wei = NULL) {
  
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
    
    ci_test <- "testIndBeta"
    a1 <- internaliamb.betabs( target = target, dataset = dataset, threshold = threshold, wei = wei, p = p ) 
    ind <- 1:p
    a2 <- list()
    poies <- a1$mat[, 1]
    if ( length(poies) > 0 ) {
      ind[-poies] <- 0
      ind <- ind[ind > 0]
      dat <- dataset[, poies ]
      a2 <- internaliamb.betabs(target = target, dataset = dat, threshold = threshold, wei = wei, p = length(ind) )  
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
      a2 <- internaliamb.betabs( target = target, dataset = dat, threshold = threshold, wei = wei, p = length(ind) ) 
      poies <- a2$mat[, 1]
      if ( length(poies) > 0 ) {
        ind[-poies] <- 0
        ind <- ind[ind > 0]
        if ( length(poies) == 1 )   dat <- dat  else   dat <- dat[, poies]
      } else  dat <- NULL  
    }
    
    runtime <- proc.time() - runtime
    # if ( !is.null(a2$mat) ) {
    #   final <- beta.mod( target, dataset[, ind], wei = wei )
    #   a2$mat[, 1] <- ind 
    # } else  final <- "No variables were selected"
    # 
    res <- list(info = ind, mat = a2$mat, ci_test = ci_test, final = a2$final ) 
  }
  
  res
}  



internaliamb.betabs <- function(target, dataset, threshold, wei, p) {
  if ( !is.null(dataset) |  p > 0 ) {
     ini <- beta.reg( target, dataset, wei = wei )
     dofini <- length( ini$be )
     likini <- ini$loglik 
     stat <- dof <- numeric(p)
     if (p > 1) {
       for (j in 1:p) {
        mod <- beta.reg( target, dataset[, -j], wei = wei )
        stat[j] <- 2 * abs( likini - mod$loglik )
        dof[j] <- dofini - length( mod$be ) 
       }
     } else {
       if ( is.null(wei) ) {
         mod0 <- Rfast::beta.mle(target)
       } else betamle.wei(target, wei = wei)
      stat <- 2 * abs( likini - mod0$loglik )
      dof <- dofini - 1
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
         mod1 <- beta.reg(target, dat, wei = wei )
         if ( is.null(wei) ) {
           mod0 <- Rfast::beta.mle(target)
         } else betamle.wei(target, wei = wei)
         stat <- 2 * abs( mod1$loglik - mod0$loglik )
         pval <- pchisq( stat, length( mod1$be ) - 1, lower.tail = FALSE, log.p = TRUE)
         if (pval > threshold ) {
           final <- "No variables were selected"
           mat <- NULL
         } else final <- mod1
       } else  final <- beta.reg(target, dat, wei = wei)
     }
     info <- info[ info[, 1] > 0, ]
     
  } else { 
    info <- NULL  
    mat <- NULL 
    final <- "No variables were selected"
  } 
  
  list(info = info, mat = mat, final = final) 
}  



