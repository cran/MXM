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
    info <- info[ info[, 1] > 0, , drop = FALSE]
    
  } else { 
    info <- NULL  
    mat <- NULL 
    final <- "No variables were selected"
  } 
  
  list(info = info, mat = mat, final = final) 
}  



