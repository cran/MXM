mmmb = function(target, dataset , max_k = 3 , threshold = 0.05 , test = "testIndFisher" , user_test = NULL, robust = FALSE, ncores = 1) {
  
  durat <- proc.time()  
  
  mmpcobject <- MMPC(target, dataset, max_k = max_k, threshold = threshold, test = test, user_test = user_test, robust = robust, ncores = ncores, backward = TRUE)
  varsToIterate <- mmpcobject@selectedVars;
  pct <- varsToIterate
  met <- 1:length(pct)
  d <- ncol(dataset) 
  ci_test <- test <- mmpcobject@test
  lista <- list()
  aa <- NULL
  
  if ( length(pct) > 0 ) {
    
    for ( i in met) {
      tar <- dataset[, varsToIterate[i] ];
      datas <- cbind( dataset[, -varsToIterate[i] ], target)
      res <- MMPC(tar, datas, max_k = 3, threshold = threshold, test = test, user_test = user_test, robust = robust, ncores = ncores, backward = TRUE) 
      poies <- sort( res@selectedVars )
      poies <- poies[poies != d ]
      poies[ poies >= varsToIterate[i] ] = poies[ poies >= varsToIterate[i] ] + 1
      lista[[ i ]] <- poies     
    }

    if ( length( unlist(lista) ) > 0 ) {
      
      lista <- unlist( lista[met] )
      lista <- unique(lista)
      lista <- setdiff(lista, pct )  ## remove the target and the PC set from the candidate MB set
      ina <- 1:length(lista)
      ci_test <- test <- mmpcobject@test
      av_tests <- c("testIndFisher", "testIndSpearman", "gSquare", NULL);
      
      if ( length(test) == 1 ) {  #avoid vectors, matrices etc
        test <- match.arg(test, av_tests ,TRUE);
        if (test == "testIndFisher")   {
          #an einai posostiaio target
          if ( min(target) > 0 & max(target) < 1 )  target = log( target/(1 - target) ) 
          test <- testIndFisher;
        }  else if (test == "testIndSpearman")  {
          #an einai posostiaio target
          if ( min(target) > 0 & max(target) < 1 )   target = log( target / (1 - target) ) ## logistic normal 
          target <- rank(target)
          dataset <- apply(dataset, 2, rank)  
          test <- testIndSpearman;  ## Spearman is Pearson on the ranks of the data
          robust <- FALSE
        }
        else if (test == "gSquare")  {
          test <- gSquare;
          robust <- FALSE
        }
      } else  stop('invalid test option');
      
      a <- log(threshold)
      mat <- list()
      for (i in 1:max_k)      mat[[ i ]] <- Rfast::comb_n(pct, i)
      
      for ( l in 1:length(lista) ) {
        
        for ( i in 1:max_k ) {
          condset <- mat[[ i ]]
          pval <- test(target, dataset, xIndex = lista[l], csIndex = condset[, 1], robust = robust)$pvalue
          k <- 1
          dm2 <- dim(condset)[2]
          while ( pval > a  &  k < dm2 ) {
            k <- k + 1
            pval <- test(target, dataset, xIndex = lista[l], csIndex = condset[, k], robust = robust)$pvalue
          }  ## end  while  
        }   ##  end  for ( i in 1:max_k )
        if ( pval > a )   ina[l] = 0
      }   ##  end  for ( l in 1:length(lista) )
    }   ##  end if ( length( unlist(lista) ) > 0 )
	   aa <- lista[ina]
  }   ##  end if ( length(varsToIterate) > 0 ) 
  
  runtime <- proc.time() - durat   
  
  list( mb = sort( c(mmpcobject@selectedVars[met], aa) ), ci_test = ci_test, runtime = runtime )
}