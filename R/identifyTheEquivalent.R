identifyTheEquivalent = function(equal_case, queues, target, dataset, cvar, z, test, wei, threshold, univariateModels, pvalues, hash, dataInfo, stat_hash, pvalue_hash, robust, ncores)
{
  z = t(z);
  #case 3
  #if we have more than one equivalent vars in z , we select the first one
  #for loop if we dont use the below lapply function
  #equalsY = 0;
  for ( i in 1:ncol(z) ) {
    w = z[,i];
    w = t( t(w) );
    zPrime = c(setdiff(z , w) , cvar);
    cur_results = test(target = target, dataset = dataset, xIndex = w, csIndex = zPrime, wei = wei, dataInfo = dataInfo, univariateModels = univariateModels, hash = hash, stat_hash = stat_hash, pvalue_hash = pvalue_hash, robust=robust);
    if ( cur_results$flag & (cur_results$pvalue > threshold) )  {  
      queues[[w]] = as.matrix( c(queues[[w]] , queues[[cvar]]) );
      break;
      #equalsY = equalsY+1;
    }
  }
  #cat("\nEquals Ys in the current z subset = %d",equalsY);
  return(queues);
}