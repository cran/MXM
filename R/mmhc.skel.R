mmhc.skel <- function(dataset, max_k = 3, threshold = 0.05, test = "testIndFisher", rob = FALSE, fast = FALSE, symmetry = TRUE, nc = 1, graph = FALSE) {
  ## dataset is either conitnuous or categorical data  
  ## max_k is the maximum number of variables upon which to condition
  ## threshold is the level of significance to reject the independence
  ## test can be either testIndFisher (default) or testIndSpearman for continuous data
  ## OR gSquare (default) for categorical data
  ## rob is for robust correlation
  ## nc is the number of cores to use, set to 1 by default
  
  dataset <- as.matrix(dataset)
  n <- dim(dataset)[2]
  G <- matrix(0, n, n)
  
   if ( fast ) {  
      
     pa <- proc.time()
     a <- MMPC(1, dataset, max_k = max_k, test = test, threshold = threshold, robust = rob)
     sel <- a@selectedVars
     G[1, sel] <- 1 
    
     for ( i in 2:n ) {
       ina <- 1:n
       che <- which( G[ 1:c(i - 1), i ] == 0 ) 
       ina[c(i, che)] <- 0   ;  ina <- ina[ina>0]
       a <- MMPC(dataset[, i], as.matrix(dataset[, -c(i, che)]), max_k = max_k, test = test, threshold = threshold, robust = rob)
       if ( !is.null(a) ) {
         sel <- a@selectedVars
         sela <- ina[sel]
         G[i, sela] <- 1 
       }  else {
          G[i, ] <- 0
       }
     }
     runtime <- proc.time() - pa
      
   } else {
      
     if (nc <= 1  ||  is.null(nc) ) {
        
       pa <- proc.time()
       for (i in 1:n) {
         a <- MMPC(i, dataset, max_k = max_k, test = test, threshold = threshold, robust = rob)
         sel <- a@selectedVars
         G[i, sel] <- 1 
       } 
       runtime <- proc.time() - pa
    
     }  else {
   
        pa <- proc.time() 
        cl <- makePSOCKcluster(nc)
        registerDoParallel(cl)
        sel <- numeric(n)
        mod <- foreach(i = 1:n, .combine = rbind, .export = c("MMPC") ) %dopar% {
          ## arguments order for any CI test are fixed
          sel <- numeric(n)
          a <- MMPC(i, dataset, max_k = max_k, test = test, threshold = threshold, robust = rob)
          sel[a@selectedVars] <- 1
          return(sel)
        }
        
       stopCluster(cl)
       G <- as.matrix(mod)
       runtime <- proc.time() - pa
    }
   
  }
  
  diag(G) <- 0
  
  if (symmetry ) {
    a <- which( G == 1  &  t(G) == 1 ) 
    G[ -a ] <- 0
  } else {
    G <- G + t(G)
    G[ G > 0 ] <- 1
  }
  
  info <- summary( Rfast::rowsums(G) )
  density <- sum(G) / n / ( n - 1 ) 
  
  if (is.null( colnames(dataset) ) ) {
    colnames(G) <- rownames(G) <- paste("X", 1:n, sep = "")
  } else  colnames(G) <- rownames(G) <- colnames(dataset)
   
  if( graph )  plotnetwork(G, paste("Skeleton of the MMHC algorithm for", deparse( substitute(dataset) ) ) )
  
  list(runtime = runtime, density = density, info = info, G = G)
}