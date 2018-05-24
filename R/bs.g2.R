bs.g2 <- function(target, dataset, threshold = 0.05) {
  
  tic <- proc.time() 
  dataset <- as.matrix(dataset)
  target <- as.numeric(target)
  z <- cbind(target, dataset)
  all.dc <- Rfast::colrange(z, cont = FALSE)
  p <- dim(z)[2]
  ind <- 2:p
  stat <- rep( Inf, p - 1 )
  pval <- rep( -100, p - 1 )
  dof <- rep( 100, p - 1 )
  sig <- log(threshold)
  
  for ( i in 2:p )  {
    dc <- c(all.dc[1], all.dc[i], all.dc[-c(1, i)])
    mod <- Rfast::g2Test(z, x = 1, y = i, cs = ind[ ind != i ], dc = dc) 
    dof[i] <- mod$df
    stat[i] <- mod$statistic
  }  ## end for (i in 2:p)  
  pval <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
  no <- which(pval == 0)
  stat[ no ] <- stat[no]/dof[no]
  if ( length(no) < 1 ) { 
    sel <- which.max( pval[2] )
  } else sel <- which.min(stat)
  
  mat <- cbind(1:p, pval, stat)
  colnames(mat) <- c("variable", "log.p-values", "statistic" )
  info <- matrix( c(0, -10, -10), ncol = 3 )
  colnames(info) <- c("Variables", "log.p-values", "statistic")
  
  if ( pval[sel] < sig ) {
    runtime <- proc.time() - tic
    res <- list(info = matrix(0, 0, 3), mat = mat, ci_test = "gSquare") 
    
  } else {
    
    info[1, ] <- mat[sel, ]
    mat <- mat[-sel, , drop = FALSE] 
    dat <- dataset[, -sel, drop = FALSE] 
    dc2 <- all.dc[-sel]
    p <- p - 1
    ind <- 2:p
    stat <- rep( Inf, p - 1 )
    pval <- rep( -100, p - 1 )
    dof <- rep( 100, p - 1 )
    for ( i in 2:p )  {
      dc <- c(dc2[1], dc2[i], dc2[-c(1, i)])
      mod <- Rfast::g2Test(z, x = 1, y = i, cs = ind[ ind != i ], dc = dc)           
      stat[i] <- mod$statistic
      dof[i] <- mod$df
    }  ## end for (i in 2:p)  
    pval <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    no <- which(pval == 0)
    stat[ no ] <- stat[no]/dof[no]
    if ( length(no) < 1 ) { 
      sel <- which.max( pval )
    } else sel <- which.min(stat)
    
    while ( pval[sel] > sig  &  p > 1 ) {
      info <- rbind(info, mat[sel, ])
      mat <- mat[-sel, , drop = FALSE] 
      dat <- dataset[, -sel, drop = FALSE] 
      dc2 <- all.dc[-sel]
      p <- p - 1
      
      if ( p == 1 ) {
        mat <- NULL
        runtime <- proc.time() - tic
        
      } else {
      
        ind <- 2:p
        stat <- rep( Inf, p - 1 )
        pval <- rep( -100, p - 1 )
        dof <- rep( 100, p - 1 )

        for ( i in 2:p )  {
          dc <- c(dc2[1], dc2[i], dc2[-c(1, i)])
          mod <- Rfast::g2Test(z, x = 1, y = i, cs = ind[ ind != i ], dc = dc)           
          stat[i] <- mod$statistic
          dof[i] <- mod$df
        }  ## end for (i in 2:p)     
        pval <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
        no <- which(pval == 0)
        stat[ no ] <- stat[no]/dof[no]
        if ( length(no) < 1 ) { 
          sel <- which.max( pval )
        } else sel <- which.min(stat)

      }  ## end else of if (p == 1)
      
    }  ## end  while ( pval[sel] > sig  &  p > 1 )
    runtime <- proc.time() - tic
    info <- info[ info[, 1] > 0, , drop = FALSE]
    res <- list(runtime = runtime, info = info, mat = mat, ci_test = "gSquare" ) 
  }  ## end else of if ( pval[sel] < sig ) 
  
  res
}  
  
  
  