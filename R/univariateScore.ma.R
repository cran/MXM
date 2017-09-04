univariateScore.ma = function(targ, data , test, statistic, dataInfo, hash, stat_hash, pvalue_hash, targetID, robust, ncores)
{
  #how many tests
  nTests = ncol(data[[ 1 ]]);
  #data structure to be returned
  univariateModels = NULL;
  univariateModels$pvalue = numeric(nTests) 
  univariateModels$stat = numeric(nTests)
  univariateModels$flag = numeric(nTests) 
  #univariateModels$uniModelFit = rep(NA,nTests);
  test_results = NULL;
  #for way to initialize the univariateModel
  if ( ncores == 1 | is.null(ncores) | ncores <= 0 ) {
    
    for(i in 1:nTests)  {
      #arguments order for any CI test are fixed
      if (i != targetID){
        test_results = test(targ, data, i, 0, statistic = statistic, dataInfo = dataInfo, hash = hash, stat_hash = stat_hash, pvalue_hash = pvalue_hash, robust = robust)
        univariateModels$pvalue[[i]] = test_results$pvalue;
        univariateModels$stat[[i]] = test_results$stat;
        univariateModels$flag[[i]] = test_results$flag;
        univariateModels$stat_hash = test_results$stat_hash
        univariateModels$pvalue_hash = test_results$pvalue_hash      
      }else{
        univariateModels$pvalue[[i]] = log(1);
        univariateModels$stat[[i]] = 0;
        univariateModels$flag[[i]] = 1;
      }
    }
    
  } else {
    # require(doParallel, quietly = TRUE, warn.conflicts = FALSE)  
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    test = test
    mod <- foreach(i = 1:nTests, .combine = rbind) %dopar% {
      ## arguments order for any CI test are fixed
      if (i != targetID) {
        test_results = test(targ, data, i, 0, statistic = statistic, dataInfo = dataInfo, hash = FALSE, stat_hash = NULL, pvalue_hash = NULL, robust = robust)
        return( c(test_results$pvalue, test_results$stat, test_results$flag) )
      } else{
        return( c(log(1), 0, 1) )
      }
    }
    stopCluster(cl)
    univariateModels$pvalue = as.vector( mod[, 1] )
    univariateModels$stat = as.vector( mod[, 2] )
    univariateModels$flag = as.vector( mod[, 3] )
    univariateModels$stat_hash = NULL
    univariateModels$pvalue_hash = NULL   
  }
  return(univariateModels);
}
