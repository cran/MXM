SES.temporal = function(target, reps = NULL, group, dataset, max_k = 3, threshold = 0.05, test = NULL, ini = NULL, wei = NULL, user_test = NULL, 
                        hash=FALSE, hashObject=NULL, slopes = FALSE, ncores = 1, logged = FALSE)
{
  #get the log threshold
  threshold = log(threshold)
  ##############################
  # initialization part of SES #
  ##############################
  equal_case = 3;
  stat_hash = NULL;
  pvalue_hash = NULL;
  
  if ( hash )  {
    if ( requireNamespace("hash") ) {
      if(is.null(hashObject))
      {
        stat_hash = hash();
        pvalue_hash = hash();
      } else if(class(hashObject) == "list"){
        stat_hash = hashObject$stat_hash;
        pvalue_hash = hashObject$pvalue_hash;
      } else   stop('hashObject must be a list of two hash objects (stat_hash, pvalue_hash)')
      
    } else {
      cat('The hash version of SES requires the hash package');
      return(NULL);
    }
  }
  
  dataInfo = NULL;
  ###################################
  # dataset checking and initialize #
  ###################################
  if (is.null(dataset) || is.null(target) ) {
      stop('invalid dataset or target (class feature) arguments.');
  } else  target = target;
  #check for NA values in the dataset and replace them with the variable mean
  if ( any(is.na(dataset)) ) {
    warning("The dataset contains missing values and they were replaced automatically by the variable (column) median.")
    dataset = apply(dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) });
  }
  ##################################
  # target checking and initialize #
  ##################################
  targetID = -1;
  #check if the target is a string
  if (is.character(target) & length(target) == 1) {
    findingTarget <- target == colnames(dataset);#findingTarget <- target %in% colnames(dataset);
    if ( !sum(findingTarget) == 1 ) {
      warning('Target name not in colnames or it appears multiple times');
      return(NULL);
    }
    targetID <- which(findingTarget);
    target <- dataset[ , targetID];
  }
  
  test = "testIndGLMM";
  #checking if target is a single number
  if (is.numeric(target) & length(target) == 1) {
    if (target > dim(dataset)[2]){
      warning('Target index larger than the number of variables');
      return(NULL);
    }
    targetID <- target;
    target <- dataset[ , targetID];
  }
  ################################
  # test checking and initialize #
  ################################
  la <- length( unique( as.numeric(target) ) )
  
  if (typeof(user_test) == "closure") {
    test = user_test;
  } else {
    #auto detect independence test in case of not defined by the user and the test is null or auto
    if (is.null(test) || test == "auto")  test = "testIndGLMM"
    #if target is a factor then use the Logistic test
    
    if ( la == 2) {
      dataInfo$target_type = "binary";
      cat('\nTarget variable type: Binary')
    } else if (sum( floor(target) - target ) == 0 ) {
      dataInfo$target_type = "discrete";
      cat('\nTarget variable type: Discrete')
    } else {
      dataInfo$target_type = "normal";
      cat('\nTarget variable type: Continuous')
    }
    cat("\nConditional independence test used: ");cat(test);cat("\n");
    #available conditional independence tests
    av_tests = c("testIndGLMM", "auto",  NULL);
    
    ci_test = test
    if (length(test) == 1) {   #avoid vectors, matrices etc
      test = match.arg(test , av_tests ,TRUE);
      #convert to closure type
      if (test == "testIndGLMM") {
        test = testIndGLMM;
      }
      
    } else   stop('invalid test option');
  }
  ###################################
  # options checking and initialize #
  ###################################
  #extracting the parameters
  max_k = floor(max_k);
  varsize = dim(dataset)[[2]];
  #option checking
  if ( (typeof(max_k) != "double") || max_k < 1 )   stop('invalid max_k option');
  if (max_k > varsize)   max_k = varsize;
  if ( (typeof(threshold) != "double") || exp(threshold) == 0 || exp(threshold) > 1 )   stop('invalid threshold option');
  if ( typeof(equal_case) != "double" )   stop('invalid equal_case option');
  ######################################################################################
  if( !is.null(user_test) )  ci_test = "user_test";
  #######################################################################################
  #call the main SES function after the checks and the initializations
  results = InternalSES.temporal(target, reps, group, dataset, max_k, threshold, test, ini, wei, equal_case, user_test, dataInfo, hash, varsize, stat_hash, pvalue_hash, targetID, slopes, ncores, logged = logged);
  SES.temporal.output <-new("SES.temporal.output", selectedVars = results$selectedVars, selectedVarsOrder=results$selectedVarsOrder, queues=results$queues, signatures=results$signatures, hashObject=results$hashObject, pvalues=results$pvalues, stats=results$stats, univ = results$uni, max_k=results$max_k, threshold = results$threshold, runtime=results$runtime, test=ci_test, slope = results$slope);
  
  return(SES.temporal.output);
}




