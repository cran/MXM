wald.Internalses = function(target, dataset, max_k, threshold, test = NULL, ini, wei=NULL, equal_case=3, user_test=NULL, dataInfo=NULL, hash=FALSE, varsize, stat_hash, pvalue_hash, targetID, robust, ncores)
{
  #get the current time
  runtime = proc.time();

  if ( is.null(ini) ) { 
    univariateModels = wald.univregs(target = target, dataset = dataset, targetID = targetID, test = test, user_test = user_test, wei = wei, dataInfo = dataInfo, ncores = ncores) 
  } else  univariateModels = ini
  
  pvalues = univariateModels$pvalue;      
  stats = univariateModels$stat;
  flags = univariateModels$flag;
  #   stat_hash = univariateModels$stat_hash;
  #   pvalue_hash = univariateModels$pvalue_hash;
  #if we dont have any associations , return
  if ( min(pvalues, na.rm = TRUE) > threshold )  {
    #cat('No associations!');
    results = NULL;
    results$selectedVars = c();
    class(results$selectedVars) = "numeric";
    results$selectedVarsOrder = c();
    class(results$selectedVarsOrder) = "numeric";
    results$queues = c();
    class(results$queues) = 'list';
    results$signatures = matrix(nrow=1,ncol=1);
    class(results$signatures) = 'matrix';
    results$hashObject = NULL;
    class(results$hashObject) = 'list';
    class(results$univ) = 'list';
    results$pvalues = exp(pvalues);
    results$stats = stats;
    results$univ = univariateModels
    results$max_k = max_k;
    results$threshold = exp(threshold);
    runtime = proc.time() - runtime;
    results$runtime = runtime;
    results$rob = robust
    results$n.tests <- length(flags)
    
    return(results);
  }
  
  #Initialize the data structs
  selectedVars = numeric(varsize);
  selectedVarsOrder = numeric(varsize);
  queues = vector('list' , varsize);
  queues <- lapply(1:varsize , function(i){queues[[i]] = i;})
  #select the variable with the highest association
  #selectedVar = which(flags == 1 & stats == stats[[which.max(stats)]]);
  selectedVar = which(flags == 1 & pvalues == pvalues[[which.min(pvalues)]]);
  selectedVars[selectedVar] = 1;
  selectedVarsOrder[selectedVar] = 1; #CHANGE
  #print(paste("rep: ",0,", selected var: ",selectedVar,", pvalue = ",exp(pvalues[selectedVar])))
  #lets check the first selected var
  #cat('First selected var: %d, p-value: %.6f\n', selectedVar, pvalues[selectedVar]);
  #remaining variables to be considered
  remainingVars = numeric(varsize) + 1;
  remainingVars[selectedVar] = 0;
  remainingVars[pvalues > threshold] = 0;
  if (targetID > 0)   remainingVars[targetID] = 0;
  #main SES loop
  #loop until there are not remaining vars
  loop = any( as.logical(remainingVars) );
  #rep = 1;
  while (loop)  {
    #lets find the equivalences
    IdEq_results <- IdentifyEquivalence(equal_case, queues, target, dataset, test, wei, threshold , max_k , selectedVars , pvalues , stats , remainingVars , univariateModels, selectedVarsOrder, hash=hash, dataInfo, stat_hash=stat_hash, pvalue_hash=pvalue_hash, robust = robust, ncores = ncores);
    queues = IdEq_results$queues;
    selectedVars = IdEq_results$selectedVars;
    remainingVars = IdEq_results$remainingVars;
    pvalues = IdEq_results$pvalues;
    stats = IdEq_results$stats;
    stat_hash=IdEq_results$stat_hash;
    pvalue_hash=IdEq_results$pvalue_hash;
    #lets find the variable with the max min association
    max_min_results = max_min_assoc(target, dataset, test, wei, threshold, max_k, selectedVars, pvalues, stats, remainingVars, univariateModels, selectedVarsOrder, hash = hash, dataInfo, stat_hash = stat_hash, pvalue_hash = pvalue_hash, robust = robust, ncores = ncores);
    selectedVar = max_min_results$selected_var;
    selectedPvalue = max_min_results$selected_pvalue;
    remainingVars = max_min_results$remainingVars;
    pvalues = max_min_results$pvalues;
    stats = max_min_results$stats;
    stat_hash=max_min_results$stat_hash;
    pvalue_hash=max_min_results$pvalue_hash;
    #if the selected variable is associated with target , add it to the selected variables
    if (selectedPvalue <= threshold)  {
      #print(paste("rep: ",rep,", selected var: ",selectedVar,", pvalue = ",exp(selectedPvalue)))
      #rep = rep + 1;
      selectedVars[selectedVar] = 1;
      selectedVarsOrder[selectedVar] = max(selectedVarsOrder) + 1;
      remainingVars[selectedVar] = 0;
    }
    
    loop = any(as.logical(remainingVars));
  }
  #lets find the variables to be discarded
  IdEq_results <- IdentifyEquivalence(equal_case, queues, target, dataset, test, wei, threshold, max_k, selectedVars, pvalues, stats, remainingVars, univariateModels, selectedVarsOrder, hash=hash, dataInfo, stat_hash = stat_hash, pvalue_hash = pvalue_hash, robust = robust, ncores = ncores);
  queues = IdEq_results$queues;
  selectedVars = IdEq_results$selectedVars;
  pvalues = IdEq_results$pvalues;
  stats = IdEq_results$stats;
  remainingVars = IdEq_results$remainingVars;
  stat_hash=IdEq_results$stat_hash;
  pvalue_hash=IdEq_results$pvalue_hash;
  selectedVarsOrder[which(!selectedVars)] = varsize;#
  numberofSelectedVars = sum(selectedVars);#
  selectedVarsOrder = sort(selectedVarsOrder);#
  #   selectedVars = selectedVarsOrder[1:numberofSelectedVars];
  #queues correctness
  all_queues = queues
  queues = queues[which(selectedVars==1)];
  queues <- lapply(1:length(queues) , function(i){queues[[i]] = unique(queues[[i]]);});
  #adjusting the results
  if (targetID > 0) {
    toAdjust <- which(selectedVars > targetID);
    selectedVars[toAdjust] = selectedVars[toAdjust] + 1;
  }
  results = NULL;
  results$selectedVars = which(selectedVars == 1);
  svorder = sort(pvalues[results$selectedVars] , index.return = TRUE);
  svorder = results$selectedVars[svorder$ix];
  results$selectedVarsOrder = svorder;
  results$queues = queues;
  results$signatures = as.matrix(do.call(expand.grid, results$queues))
  hashObject = NULL;
  hashObject$stat_hash = stat_hash;
  hashObject$pvalue_hash = pvalue_hash;
  results$hashObject = hashObject;
  class(results$hashObject) = 'list';
  results$pvalues = exp(pvalues);
  results$stats = stats;
  results$univ = univariateModels
  #   results$all_queues = all_queues;
  #   already known
  #   results$data = dataset;
  #   results$target = target;
  #   results$test = test;
  results$max_k = max_k;
  results$threshold = exp(threshold);
  runtime = proc.time() - runtime;
  results$runtime = runtime;
  results$rob = robust
  results$n.tests <- length(flags) + length( hashObject$stat_hash )
  
  return(results);
}