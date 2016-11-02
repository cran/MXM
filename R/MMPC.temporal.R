MMPC.temporal = function(target , reps = NULL, group, dataset , max_k = 3 , threshold = 0.05 , test = NULL, ini = NULL, wei = NULL, user_test = NULL, hash=FALSE, hashObject=NULL, slopes = FALSE, ncores = 1)
{
 #get the log threshold
  threshold = log(threshold)
      
  ##############################
  # initialization part of MMPC #
  ##############################


  options(warn=0);

  equal_case = 3;
  stat_hash = NULL;
  pvalue_hash = NULL;
  
  if(hash == TRUE)
  {
    if(requireNamespace("hash"))
    {
      if(is.null(hashObject))
      {
        stat_hash = hash();
        pvalue_hash = hash();
      }else if(class(hashObject) == "list"){
        stat_hash = hashObject$stat_hash;
        pvalue_hash = hashObject$pvalue_hash;
      }else{
        stop('hashObject must be a list of two hash objects (stat_hash, pvalue_hash)')
      }
    }else{
      cat('The hash version of SES requires the hash package');
      return(NULL);
    }
  }
  
  dataInfo = NULL;
  
  ###################################
  # dataset checking and initialize #
  ###################################
  
  if(!is.null(dataset))
  {
    #check if dataset is an ExpressionSet object of Biobase package
    #if(class(dataset) == "ExpressionSet")
    #{
    #  #get the elements (numeric matrix) of the current ExpressionSet object.
    #  dataset = Biobase::exprs(dataset);
    #  dataset = t(dataset); #take the features as columns and the samples as rows
    #} else if ( class(dataset) == "matrix" | class(dataset) == "data.frame" ){
    #  target = target 
    #} else {
    #  stop('Invalid dataset class. It must be either a matrix, a dataframe or an ExpressionSet');
    #}
	
  }
    if(is.null(dataset) || is.null(target) )
    {
      stop('invalid dataset or target (class feature) arguments.');
    }else{
      target = target;
    }
  
  #check for NA values in the dataset and replace them with the variable mean
  if(any(is.na(dataset)) == TRUE)
  {
    dataset = as.matrix(dataset);
    warning("The dataset contains missing values and they were replaced automatically by the variable (column) mean.")
    dataset = apply(dataset, 2, function(x){ x[which(is.na(x))] = mean(x,na.rm = TRUE) ; return(x) });
  }
  
  ##################################
  # target checking and initialize #
  ##################################
  
  targetID = -1;
  
  #check if the target is a string
  if (is.character(target) & length(target) == 1){
    findingTarget <- target == colnames(dataset);#findingTarget <- target %in% colnames(dataset);
    if(!sum(findingTarget)==1){
      warning('Target name not in colnames or it appears multiple times');
      return(NULL);
    }
    targetID <- which(findingTarget);
    target <- dataset[ , targetID];
  }
  
  #checking if target is a single number
  if (is.numeric(target) & length(target) == 1){
    if(target > dim(dataset)[2]){
      warning('Target index larger than the number of variables');
      return(NULL);
    }
    targetID <- target;
    target <- dataset[ , targetID];
  }
  
  ################################
  # test checking and initialize #
  ################################
  
  if(typeof(user_test) == "closure")
  {
    test = user_test;
  }else{
    #auto detect independence test in case of not defined by the user and the test is null or auto
    if(is.null(test) || test == "auto")
    {
      #if target is a factor then use the Logistic test
      if("factor" %in% class(target))
      {
        test = "testIndGLMM";
        if(is.ordered(target) == TRUE)
        {
          dataInfo$target_type = "binary";
          
          cat('\nTarget variable type: Binary')
        }else{
          if(length(unique(target)) == 2)
          {
            dataInfo$target_type = "binary"
            cat('\nTarget variable type: Binomial')
          }else{
            dataInfo$target_type = "nominal"
            cat('\nTarget variable type: Nominal')
          }
        }
      }else if(class(target) == "numeric" || class(target) == "matrix"){
        if(class(target) == "matrix")
        {
          if(dim(target)[2]!=1)
          {
            stop('Target can not be a matrix')
          }
        }
        
        if(identical(floor(target),target) == TRUE)
        {
          if(length(unique(target)) == 2)
          {
            dataInfo$target_type = "binary";
            cat('\nTarget variable type: Binary')
            test = "testIndGLMM";
          }else{
            test = "testIndGLMM";
            dataInfo$target_type = "discrete";
            cat('\nTarget variable type: Discrete')
          }
        }else{
          dataInfo$target_type = "normal";
          cat('\nTarget variable type: Normal')
          test = "testIndGLMM";  
        }
      }else{
        stop('Target must be a factor, vector, matrix with one column or a Surv object');
      }
    }
    
    if(test == "testIndGLMM")
    {

       if(identical(floor(target),target) == TRUE)
        {
          if(length(unique(target)) == 2)
          {
            dataInfo$target_type = "binary";
            cat('\nTarget variable type: Binary')
            test = "testIndGLMM";
          }else{
            test = "testIndGLMM";
            dataInfo$target_type = "discrete";
            cat('\nTarget variable type: Discrete')
          }
        }else{
          dataInfo$target_type = "normal";
          cat('\nTarget variable type: Normal')
          test = "testIndGLMM";  
        }
      
      if(is.ordered(target) == TRUE)
      {
        dataInfo$target_type = "binary";
        #cat('\nTarget variable type: Binary')
        
        if(requireNamespace("lme4", quietly = TRUE, warn.conflicts = FALSE)==FALSE)
        {
          cat("The testIndGLMM test requires the lme4 package for the glmm. Please install it.");
          return(NULL);
        }
        
      }else{
        if(identical(floor(target),target) == TRUE)
        {
          if(length(unique(target)) == 2)
          {
            dataInfo$target_type = "binary";
            cat('\nTarget variable type: Binary')
            test = "testIndGLMM";
          }else{
            test = "testIndGLMM";
            dataInfo$target_type = "discrete";
            cat('\nTarget variable type: Discrete')
          }
        }else{
          dataInfo$target_type = "normal";
          cat('\nTarget variable type: Normal')
          test = "testIndGLMM";  
        }
        
        if(requireNamespace("lme4", quietly = TRUE, warn.conflicts = FALSE)==FALSE)
        {
          cat("The testIndGLMM test requires the lme4 package for the glmm. Please install it.");
          return(NULL);
        }
        
      }
    }
    
    cat("\nConditional independence test used: ");cat(test);cat("\n");
    
    #available conditional independence tests
    av_tests = c("testIndGLMM", "auto" ,  NULL);
    
    ci_test = test
    if(length(test) == 1) #avoid vectors, matrices etc
    {
      test = match.arg(test , av_tests ,TRUE);
      #convert to closure type
      if(test == "testIndGLMM")
      {
        #an einai posostiaio target
         if ( all(target > 0 & target < 1) ){
         target = log( target/(1 - target) ) ## logistic normal 
         }
        
        test = testIndGLMM;
      }
      
    }else{
      stop('invalid test option');
    }
  }
  
  ###################################
  # options checking and initialize #
  ###################################
  
  #extracting the parameters
  max_k = floor(max_k);
  varsize = dim(dataset)[[2]];
  
  #option checking
  if((typeof(max_k)!="double") || max_k < 1)
  {
    stop('invalid max_k option');
  }
  if(max_k > varsize)
  {
    max_k = varsize;
  }
  if((typeof(threshold)!="double") || exp(threshold) <= 0 || exp(threshold) > 1)
  {
    stop('invalid threshold option');
  }
  if(typeof(equal_case)!="double")
  {
    stop('invalid equal_case option');
  }
  
  #######################################################################################

  if(!is.null(user_test))
  {
    ci_test = "user_test";
  }

  #call the main MMPC.temporal function after the checks and the initializations
  results = InternalMMPC.temporal(target, reps, group, dataset, max_k, threshold, test, ini, wei, user_test, dataInfo, hash, varsize, stat_hash, pvalue_hash, targetID, slopes = slopes, ncores = ncores);
  
  MMPC.temporal.output <-new("MMPC.temporal.output", selectedVars = results$selectedVars, selectedVarsOrder=results$selectedVarsOrder, hashObject=results$hashObject, pvalues=results$pvalues, stats=results$stats, univ = results$univ, max_k=results$max_k, threshold = results$threshold, runtime=results$runtime, test=ci_test, slope = slopes);
  
  return(MMPC.temporal.output);
  
}

#########################################################################################################

  InternalMMPC.temporal = function(target, reps, group, dataset, max_k, threshold, test = NULL, ini, wei, user_test = NULL, dataInfo = NULL, hash=FALSE, varsize, stat_hash, pvalue_hash, targetID, slopes, ncores)
{
  #get the current time
  runtime = proc.time();
  
  #######################################################################################
  
  rows = length(target)
  cols = ncol(dataset)
  
  #univariate feature selection test
  
  if ( is.null(ini) ) {
    univariateModels = univariateScore.temporal(target , reps, group, dataset , test, wei=wei, hash=hash, dataInfo, stat_hash=stat_hash, pvalue_hash=pvalue_hash, targetID, slopes = slopes, ncores = ncores);
  } else {
    univariateModels = ini
  }
  
  pvalues = univariateModels$pvalue;
  stats = univariateModels$stat;
  flags = univariateModels$flag;
  stat_hash = univariateModels$stat_hash;
  pvalue_hash = univariateModels$pvalue_hash;
  #if we dont have any associations , return
  if(min(pvalues , na.rm=TRUE) > threshold) #or min(pvalues, na.rm=TRUE)
  {
    cat('No associations!');
    
    results = NULL;
    results$selectedVars = c();
    class(results$selectedVars) = "numeric";
    results$selectedVarsOrder = c();
    class(results$selectedVarsOrder) = "numeric";
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
    results$slope = slopes
    
    return(results);
  }
  
  
  #Initialize the data structs
  selectedVars = numeric(varsize);
  selectedVarsOrder = numeric(varsize);
  
  #select the variable with the highest association
  selectedVar = which(flags == 1 & stats == stats[[which.max(stats)]]);
  selectedVars[selectedVar] = 1;
  selectedVarsOrder[selectedVar] = 1; #CHANGE
  
  #lets check the first selected var
  #cat('First selected var: %d, p-value: %.6f\n', selectedVar, pvalues[selectedVar]);
  
  #remaining variables to be considered
  remainingVars = numeric(varsize) + 1;
  remainingVars[selectedVar] = 0;
  remainingVars[pvalues > threshold] = 0;
  if (targetID > 0){
    remainingVars[targetID] = 0;
  }
  
  ################ main MMPC.temporal loop ################
  
  #main MMPC.temporal loop
  #loop until there are not remaining vars
  loop = any(as.logical(remainingVars));
  
  while(loop)
  {
    #lets find the variable with the max min association
    max_min_results = max_min_assoc.temporal(target, reps, group, dataset , test , wei, threshold , max_k , selectedVars , pvalues , stats , remainingVars , univariateModels, selectedVarsOrder, hash=hash, dataInfo, stat_hash=stat_hash, pvalue_hash=pvalue_hash, slopes = slopes);
    selectedVar = max_min_results$selected_var;
    selectedPvalue = max_min_results$selected_pvalue;
    remainingVars = max_min_results$remainingVars;
    pvalues = max_min_results$pvalues;
    stats = max_min_results$stats;
    stat_hash=max_min_results$stat_hash;
    pvalue_hash=max_min_results$pvalue_hash;
    
    #if the selected variable is associated with target , add it to the selected variables
    if(selectedPvalue <= threshold)
    {
      selectedVars[selectedVar] = 1;
      selectedVarsOrder[selectedVar] = max(selectedVarsOrder) + 1;
      remainingVars[selectedVar] = 0;
    }
    
    loop = any(as.logical(remainingVars));
  }
  
  selectedVarsOrder[which(!selectedVars)] = varsize;#
  numberofSelectedVars = sum(selectedVars);#
  selectedVarsOrder = sort(selectedVarsOrder);#
  #   selectedVars = selectedVarsOrder[1:numberofSelectedVars];
  
  #adjusting the results
  if(targetID > 0)
  {
    toAdjust <- which(selectedVars > targetID);
    selectedVars[toAdjust] = selectedVars[toAdjust] + 1;
  }
  
  results = NULL;
  results$selectedVars = which(selectedVars == 1);
  
  svorder = sort(pvalues[results$selectedVars] , index.return = TRUE);
  svorder = results$selectedVars[svorder$ix];
  results$selectedVarsOrder = svorder;
  
  hashObject = NULL;
  hashObject$stat_hash = stat_hash;
  hashObject$pvalue_hash = pvalue_hash;
  results$hashObject = hashObject;
  class(results$hashObject) = 'list';
  class(results$univ) = 'list';
  
  results$pvalues = exp(pvalues);
  results$stats = stats;
  results$univ = univariateModels
  
  results$max_k = max_k;
  results$threshold = exp(threshold);
  
  runtime = proc.time() - runtime;
  results$runtime = runtime;
  results$slope = slopes
  
  
  return(results);
}



# 
# # .onAttach <- function(libname, pkgname){
# #   # do whatever needs to be done when the package is loaded
# #   packageStartupMessage( "Loading MXM package version 0.2, thanks for downloading." )
# #   #load the dll files from the fortran code for the package
# #   #library.dynam("MXM", pkgname, libname)
# # }
