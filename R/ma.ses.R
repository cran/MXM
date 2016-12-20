# R Implementation of the Statistical equivalence signatures algorithm (SES)
# as described in the "Discovering multiple, equivalent biomarker signatures" paper
# by Ioannis Tsamardinos, Vincenzo Lagani and Dimitris Pappas
# R Implementation by Giorgos Athineou (2013-2014)
# VERSION: 17/3/2014

# INPUTS

# target : the class variable , provide either a vector, an 1D matrix, an 1D array (same length with the data rows), a factor 
# or a formula. The data can be either continuous data in R, values within (0,1), binary [0,1], nominal or ordinal data.
# dataset : tha dataset , provide a data frame (columns = variables , rows = samples) or a matrix or an ExpressionSet.
# max_k : the maximum conditioning set which is used in the current conditional indepedence test used.
# threshold : the significance threshold ( must be in (0,1) ) for testing the null hypothesis of the generated pvalues.
# test : the conditional independence test we are going to use. the available conditional independence tests so far in this 
# implementation are:
#   "maFisher" : Fisher conditional independence test for continous targets (or proportions) and continuous predictors only
#   "testIndSpearman" : Fisher conditional independence test for continous targets (or proportions) and continuous predictors only (Spearman correlation is calculated first)
#   "testIndReg" : Conditional independence test based on regression for continous targets (or proportions) and mixed predictors using the F test
#   "testIndRQ" : Conditional Independence Test based on quantile (median) regression for numerical class variables and mixed predictors (F test)
#   "testIndLogistic" : Conditional Independence Test based on logistic regression for binary,categorical or ordinal class variables and mixed predictors
#   "testIndPois" : Conditional Independence Test based on Poisson regression for discrete class variables and mixed predictors (log-likelihood ratio test)
#   "testIndSpeedglm" : Conditional Independence Test based on linear, binary logistic and poisson regression with mixed predictors (log-likelihood ratio test)
#   "testIndZIP" : Conditional Independence Test based on zero inflated poisson regression for discrete class variables and mixed predictors (log-likelihood ratio test)
#   "testIndNB" : Conditional Independence Test based on negative binomial regression for discrete class variables and mixed predictors (log-likelihood ratio test)
#   "testIndBeta" : Conditional Independence Test based on beta regression for proportions and mixed predictors (log likelihood ratio test)
#   "testIndMVreg" : Conditional Independence Test based on mu;ltivariate linear regression for Euclidean data and mixed predictors (log likelihood ratio test)
#   "gSquare" : Conditional Independence test based on the G test of independence (log likelihood ratio  test)
#   "censIndCR" : Conditional independence test for survival data based on the Log likelihood ratio test with mixed predictors (Cox regression)
# user_test : the user defined conditional independence test ( provide a closure type object )
# hash : a boolean variable whuch indicates whether (TRUE) or not (FALSE) to use the hash-based implementation of the statistics of SES.
# hashObject : a List with the hash objects (hash package) on which we have cached the generated statistics. 
#              SES requires this Object for the hash-based implementation of the statistics. This hashObject is produced or 
# updated by each run
#              of SES (if hash == TRUE) and it can be reused in next runs of SES.
# robust : Should the Fisher or the normal regression be robustly estimated? TRUE or FALSE 

# there are default values for all of the parameters of the algorithm.

# OUTPUT <LIST>
# The output of the algorithm is a LIST with the following quantities (14) :

# selectedVars : the selected variables i.e. the dependent of the target variables.
# selectedVarsOrder : the increasing order of the selected variables due to their pvalues
# queues : the selected variable queues with the multiple statistically equivalent variables 
# (if you want to make multiple statistically equivalent signatures you 
# have to take one variable from each queue).
# signatures : all the possible combinations of the variables in the queues. One variable per queue in each signature. 
# (signature ~ the minimum subset with the most relevant features).
# hashObject : the hashObject with the cached statistic results of the current run.

# pvalues : the pvalues of all of the variables.
# stats : the stats of all of the variables.
# all_queues : the queues of all of the variables.

# data : the dataset used in the current run.
# target : the class variable used in the current run.
# test : the conditional independence test used in the current run.
# max_k : the max_k option used in the current run.
# threshold : the threshold option used in the current run.

# runtime : the run time of the algorithm.

# Conditional independence test arguments have to be in this exact fixed order : 
# target(target variable), data(dataset), xIndex(x index), csIndex(cs index), dataInfo(list), 
# univariateModels(cached statistics for the univariate indepence test), hash(hash booleab), stat_hash(hash object), 
# pvalue_hash(hash object), robust=robust
# example: test(target, data, xIndex, csIndex, dataInfo=NULL, univariateModels=NULL, hash=FALSE, stat_hash=NULL, pvalue_hash=NULL, robust)
# output of each test: LIST of the generated pvalue, stat, flag and the updated hash objects.

# equal_case variable inside the code : it determines the method of the equivalent estimation
#   if equal_case = 1 then, if we have more than one equivalent vars in z , we select the one with the most closer pvalue to the pvalue of cvar
#   if equal_case = 2 then, if we have more than one equivalent vars in z , we select the one with the most minimum pvalue (>a)
#   else in any other case, if we have more than one equivalent vars in z , we select the first one
# In this version we support the equal_case = 3.
# 
# 
# #hashObject
# library(hash)

ma.ses = function(target, dataset, ina, statistic = FALSE, max_k = 3 , threshold = 0.05 , test = NULL , ini = NULL, user_test = NULL, hash=FALSE, hashObject=NULL, robust = FALSE, ncores = 1)
{
  #get the log threshold
  threshold = log(threshold)
  
  ##############################
  # initialization part of SES #
  ##############################
  
  equal_case = 3;
  stat_hash = NULL;
  pvalue_hash = NULL;
  
  if( hash )
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
      }else  stop('hashObject must be a list of two hash objects (stat_hash, pvalue_hash)')
    }else{
      cat('The hash version of ma.ses requires the hash package');
      return(NULL);
    }
  }
  
  dataInfo = NULL;
  
  ###################################
  # dataset checking and initialize #
  ###################################
  

  ## hey user, make sure everything is OK

    if(is.null(dataset) || is.null(target)) #|| (dim(as.matrix(target))[2] != 1 & class(target) != "Surv" ))
    {
      stop('invalid dataset or target (class feature) arguments.');
    }else{
      target = target;
    }
  
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
  
 
    #cat("\nConditional independence test used: ");cat(test);cat("\n");
    
    #available conditional independence tests
    av_tests = c("testIndFisher", "testIndSpearman", "testIndReg", "testIndRQ", "testIndBeta", "censIndCR","censIndWR", "testIndClogit", "testIndLogistic", "testIndPois", "testIndNB", "testIndBinom", "gSquare", "auto" , "testIndZIP" , "testIndSpeedglm", "testIndMVreg", NULL);
    
    ci_test = test
    #cat(test)
    
    if(length(test) == 1)  {#avoid vectors, matrices etc

      test = match.arg(test, av_tests, TRUE);
      #convert to closure type
      if(test == "testIndFisher") 
      {
        #an einai posostiaio target
        if ( min(target) > 0  &  max(target) < 1 )   target = log( target/(1 - target) ) ## logistic normal 
        
        if(class(dataset) == "data.frame")
        {
          if ( length( which_isFactor(dataset) ) > 0 ) {
            warning("Dataset contains categorical variables (factors). A regression model is advised to be used instead.")
          } else dataset <- as.matrix(dataset)
        }
            
        test = testIndFisher;
      }
      else if(test == "testIndSpearman")
      {
        #an einai posostiaio target
        if ( min(target) > 0  &  max(target) < 1 )   target = log( target / (1 - target) ) ## logistic normal 
        
        if(class(dataset) == "data.frame")
        {
          if ( length( which_isFactor(dataset) ) > 0 ) {
            warning("Dataset contains categorical variables (factors). A regression model is advised to be used instead.")
          } else dataset <- as.matrix(dataset)
        }
        target = rank(target)
        dataset = apply(dataset, 2, rank)  
        test = testIndSpearman;  ## Spearman is Pearson on the ranks of the data
      }
	    
    } else   stop('invalid test option');
  
  
  ###################################
  # options checking and initialize #
  ###################################
  
  #extracting the parameters
  max_k = floor(max_k);
  varsize = ncol(dataset);
  
  #option checking
  if((typeof(max_k)!="double") || max_k < 1)  stop('invalid max_k option');
  if(max_k > varsize)  max_k = varsize;
  if((typeof(threshold)!="double") || exp(threshold) <= 0 || exp(threshold) > 1)  stop('invalid threshold option');
  if(typeof(equal_case)!="double")  stop('invalid equal_case option');
  
  #######################################################################################
  
  if ( !is.null(user_test) )  ci_test = "user_test";

  ## end of checking
  D = max(ina)
  
  targ = list()
  data = list()
  for (l in 1:D) {  
   targ[[ l ]] = target[ ina == l ]
   data[[ l ]] = dataset[ ina == l, ]
  }   

  #call the main SES function after the checks and the initializations
  results = Internalmases(targ, data, statistic, max_k, threshold, test, ini, equal_case, user_test, dataInfo, hash, varsize, stat_hash, pvalue_hash, targetID, robust = robust, ncores = ncores);
  
  masesoutput <-new("mases.output", selectedVars = results$selectedVars, selectedVarsOrder=results$selectedVarsOrder, queues=results$queues, signatures=results$signatures, hashObject=results$hashObject, pvalues=results$pvalues, stats=results$stats, univ = results$univ, max_k=results$max_k, threshold = results$threshold, runtime=results$runtime, test=ci_test, rob = robust);
  
  return(masesoutput);
  
}


########################
########################


Internalmases = function(targ , data, statistic, max_k, threshold , test = NULL , ini, equal_case = 3 , user_test = NULL , dataInfo = NULL , hash=FALSE, varsize, stat_hash, pvalue_hash, targetID, robust, ncores)
{
  #get the current time
  runtime = proc.time();
  
  #######################################################################################
  
  ## if testIndFisher is selected and the user has put robust =TRUE, the robust is not taken into consideration. 
  #univariate feature selection test
  
  D = length(targ) 
    
  if ( is.null(ini) ) { 

    cols = ncol( data[[ 1 ]] )  
    
	sta = pval = matrix(0, D, cols)
	
    univariateModels = list()
    
    if ( identical(test, testIndFisher)  &  !robust )  ## Pearson's correlation 
    {
	
	  for ( l in 1:D ) {  
    
	  target = targ[[ l ]]
	  dataset = data[[ l ]]
      rows = length(target)
      a = as.vector( cor(target, dataset) )
      dof = rows - 3; #degrees of freedom
      wa = abs( 0.5 * log( (1 + a) / (1 - a) ) * sqrt(dof) )
      
      if ( targetID != - 1 )  wa[ targetID ] = 0  
       sta[l, ] = wa
       pval[l, ] = log(2) + pt(-wa, dof, log.p = TRUE)
	 }  
      univariateModels$stat = - 2 * Rfast::colsums(pval);
      univariateModels$pvalue = pchisq(univariateModels$stat, 2 * D, lower.tail = FALSE, log.p = TRUE) ;
      univariateModels$flag = numeric(cols) + 1;
      univariateModels$stat_hash = stat_hash;
      univariateModels$pvalue_hash = pvalue_hash;
	 
    } else if ( identical(test, testIndSpearman) ) {  ## Spearman's correlation
	
	 for ( l in 1:D ) {  
    
	  target = targ[[ l ]]
	  dataset = data[[ l ]]
      rows = length(target) 
      a = as.vector( cor(target, dataset) )
      dof = rows - 3; #degrees of freedom
      wa = abs( 0.5 * log( (1 + a) / (1 - a) ) * sqrt(dof) / 1.029563 )
      
      if ( targetID != - 1 )  wa[ targetID ] = 0
	  sta[l, ] = wa
      pval[l, ] = log(2) + pt(-wa, dof, lower.tail = FALSE, log.p = TRUE)
	 }  
      univariateModels$stat = - 2 * Rfast::colsums(pval);
      univariateModels$pvalue = pchisq(univariateModels$stat, 2 * D, lower.tail = FALSE, log.p = TRUE) ;
      univariateModels$flag = numeric(cols) + 1;
      univariateModels$stat_hash = stat_hash;
      univariateModels$pvalue_hash = pvalue_hash;    
	 
    } else if ( identical(test, testIndReg)  & robust ) {  ## M (Robust) linear regression
      
      univariateModels = list();
      fit1 = MASS::rlm(target ~ 1)
      lik2 = numeric(cols)
      dof = numeric(cols)
      
      if ( ncores <= 1 | is.null(ncores) ) {
        
      for (l in 1:D) {

      target = targ[[ l ]]
	  dataset = data[[ l ]]

        for ( i in 1:cols ) {
          
          if ( i != targetID ) {
            
            fit2 = MASS::rlm(target ~ dataset[, i] )
            lik2[i] = as.numeric( logLik(fit2) )
            dof[i] = length( coef(fit2) ) - 1
          } else {
            lik2[i] = lik1
          }
          
        } 
        
        lik1 = as.numeric( logLik(fit1) )
        sta[l, ] = as.vector( 2 * abs(lik1 - lik2) )
        pval[l, ] = pchisq(sta[l, ], dof, lower.tail = FALSE, log.p = TRUE)

       } 

        univariateModels$stat = - 2 * Rfast::colsums(pval)
        univariateModels$pvalue = pchisq(univariateModels$stat, 2 * D, lower.tail = FALSE, log.p = TRUE) ;
        univariateModels$flag = numeric(cols) + 1;
        univariateModels$stat_hash = stat_hash;
        univariateModels$pvalue_hash = pvalue_hash;
        
      } else {
       
       for ( l in 1:D )	{
	   
	   target = targ[[ l ]]
	   dataset = data[[ l ]]
	  
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mod <- foreach(i = 1:cols, .combine = rbind, .packages = "MASS") %dopar% {
          
          if ( i != targetID ) {
            fit2 = rlm(target ~ dataset[, i] )
            lik2 = as.numeric( logLik(fit2) )
            
            return( c(lik2, length( coef(fit2) ) ) )
          } else{
            return( c(0, 0) )
          }
          
        }
        stopCluster(cl)
        
        lik1 = as.numeric( logLik(fit1) )
        dof = as.vector( mod[, 2] ) - 1 
        sta[l, ] = as.vector( 2 * abs(lik1 - mod[, 1]) )
        pval[l, ] = pchisq(sta[l, ], dof, lower.tail = FALSE, log.p = TRUE)
		
       }
	   
        univariateModels$stat = - 2 * Rfast::colsums(pval)
        univariateModels$pvalue = pchisq(univariateModels$stat, 2 * D, lower.tail = FALSE, log.p = TRUE) ;
        univariateModels$flag = numeric(cols) + 1
        univariateModels$stat_hash = NULL
        univariateModels$pvalue_hash = NULL   
      }   
      
      
    } else if ( identical(test, testIndLogistic)  &  is.ordered(target) ) {  ## 
      
      lik2 = numeric(cols)
      dof = numeric(cols)
      univariateModels = list();
      
      fit1 = ordinal::clm(target ~ 1)
      df1 = length( coef(fit1) )
      
      if ( ncores <= 1 | is.null(ncores) ) {
       
       for ( l in 1:D )	{
	   
	   target = targ[[ l ]]
	   dataset = data[[ l ]]
	   
        for ( i in 1:cols ) {
          
          if (i != targetID){
            
            x <- model.matrix(target ~ dataset[, i] )
            fit2 <- ordinal::clm.fit(target, x)
            lik2[i] <- as.numeric( fit2$logLik )
            dof[i] <- length( coef(fit2) ) - df1
          } else {
            lik2[i] <- lik1
          }   
        }
        
        lik1 = as.numeric( logLik(fit1) )
        sta[l, ] = as.vector( 2 * abs(lik1 - lik2) )
        pval[l, ] = pchisq(sta[l, ], dof, lower.tail = FALSE, log.p = TRUE)
       
       }
 	   
        univariateModels$stat = - 2 * Rfast::colsums(pval)
        univariateModels$pvalue = pchisq(univariateModels$stat, 2 * D, lower.tail = FALSE, log.p = TRUE) ;
        univariateModels$flag = numeric(cols) + 1;
        univariateModels$stat_hash = stat_hash;
        univariateModels$pvalue_hash = pvalue_hash;
        
      } else {
        
       for ( l in 1:D )	{
	   
	   target = targ[[ l ]]
	   dataset = data[[ l ]]
	   
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mod <- foreach(i = 1:cols, .combine = rbind, .packages = "ordinal") %dopar% {
          ## arguments order for any CI test are fixed
          if ( i != targetID ) {
            x <- model.matrix(target ~ dataset[, i] )
            fit2 <- ordinal::clm.fit(target, x)
            lik2 <- as.numeric( fit2$logLik )
            
            return( c(lik2, length( coef(fit2) ) ) )
          } else{
            return( c(0, 0) )
          }
        }
        stopCluster(cl)
        
        lik1 = as.numeric( logLik(fit1) )
        dof = as.vector( mod[, 2] ) - df1 
        sta[l, ] = as.vector( 2 * abs(lik1 - mod[, 1]) )
        pval[l, ] = pchisq(sta[l, ], dof, lower.tail = FALSE, log.p = TRUE)
		
	   }
	   
        univariateModels$stat = - 2 * Rfast::colsums(pval)
        univariateModels$pvalue = pchisq(univariateModels$stat, 2 * D, lower.tail = FALSE, log.p = TRUE) ;
        univariateModels$flag = numeric(cols) + 1
        univariateModels$stat_hash = NULL
        univariateModels$pvalue_hash = NULL   
      }
      
    } else {  
      univariateModels = univariateScore.ma(targ, data, test, statistic = FALSE, dataInfo = dataInfo, hash=hash, stat_hash=stat_hash, pvalue_hash=pvalue_hash, targetID=targetID, robust=robust, ncores=ncores);
    }
    
  } else {
    univariateModels = ini
  } 
    
  pvalues = univariateModels$pvalue;      
  stats = univariateModels$stat;
  flags = univariateModels$flag;
#   stat_hash = univariateModels$stat_hash;
#   pvalue_hash = univariateModels$pvalue_hash;
  #if we dont have any associations , return
  if ( min(pvalues, na.rm=TRUE) > threshold ) 
  {
    cat('No associations!');
    
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
  if (targetID > 0){
    remainingVars[targetID] = 0;
  }
  
  ################ main ses loop ################
  
  #main SES loop
  #loop until there are not remaining vars
  loop = any(as.logical(remainingVars));
  
  #rep = 1;
  while(loop)
  {
    #lets find the equivalences
    IdEq_results <- IdentifyEquivalence.ma(equal_case, queues, targ, data, test, threshold, statistic, max_k, selectedVars, pvalues, stats, remainingVars, univariateModels, selectedVarsOrder, hash=hash, dataInfo, stat_hash=stat_hash, pvalue_hash=pvalue_hash, robust = robust, ncores = ncores);
    queues = IdEq_results$queues;
    selectedVars = IdEq_results$selectedVars;
    remainingVars = IdEq_results$remainingVars;
    pvalues = IdEq_results$pvalues;
    stats = IdEq_results$stats;
    stat_hash=IdEq_results$stat_hash;
    pvalue_hash=IdEq_results$pvalue_hash;
    
    #lets find the variable with the max min association
    max_min_results = max_min_assoc.ma(targ, data, test, threshold, statistic, max_k, selectedVars, pvalues, stats, remainingVars, univariateModels, selectedVarsOrder, hash=hash, dataInfo, stat_hash=stat_hash, pvalue_hash=pvalue_hash, robust = robust, ncores = ncores);
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
      #print(paste("rep: ",rep,", selected var: ",selectedVar,", pvalue = ",exp(selectedPvalue)))
      #rep = rep + 1;
      
      selectedVars[selectedVar] = 1;
      selectedVarsOrder[selectedVar] = max(selectedVarsOrder) + 1;
      remainingVars[selectedVar] = 0;
    }
    
    loop = any(as.logical(remainingVars));
  }
  
  #lets find the variables to be discarded
  IdEq_results <- IdentifyEquivalence.ma(equal_case, queues, targ, data, test, threshold, statistic, max_k, selectedVars, pvalues, stats, remainingVars, univariateModels, selectedVarsOrder, hash=hash, dataInfo, stat_hash=stat_hash, pvalue_hash=pvalue_hash, robust = robust, ncores = ncores);
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
  
  queues <- lapply( 1:length(queues) , function(i) { queues[[ i ]] = unique( queues[[i]] ) } )
  
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
  
  
  return(results);
}






















#univariate feature selection ( uncoditional independence )

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

    for(i in 1:nTests)
    {
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



#########################################################################################################

IdentifyEquivalence.ma = function(equal_case, queues, target, dataset, test, threshold, statistic, max_k, selectedVars, pvalues, stats, remainingVars, univariateModels, selectedVarsOrder, hash, dataInfo, stat_hash, pvalue_hash, robust = robust, ncores = ncores)
{ 
  varsToBeConsidered = which(selectedVars==1 | remainingVars==1); #CHANGE
  lastvar = which(selectedVarsOrder == max(selectedVarsOrder))[1]; #CHANGE
  
  #for every variable to be considered
  for(cvar in varsToBeConsidered)
  {
    #if var is the last one added, no sense to perform the check
    if (cvar == lastvar)  next;
    
    #the maximum conditioning set
    selectedVarsIDs = which(selectedVars == 1);
    cs = setdiff(selectedVarsIDs , cvar);
    k = min( c(max_k , length(cs)) );
    
    #for all possible subsets at most k size
    #problem in 1:k if K=0 - solve with while temporary
    klimit = 1;
    while (klimit <= k)
    {
      #set to investigate
      tempCS = setdiff(cs, lastvar)#CHANGE
      if(klimit == 1) #CHANGE
      {
        subsetcsk = as.matrix(lastvar); #CHANGE
      }else{
        subsetcsk = as.matrix( nchoosek(tempCS, klimit - 1) )   
        numSubsets = dim(subsetcsk)[2]; 
        subsetcsk = rbind(subsetcsk, lastvar*rep(1,numSubsets));#CHANGE
      }
      
      #flag to get out from the outermost loop
      breakFlag = FALSE;
      

      for (i in 1:ncol(subsetcsk) )
      {
        z = subsetcsk[,i];
        z = t(t(z));
        
        cur_results = test(target = target, dataset = dataset, xIndex = cvar, csIndex = z, statistic = statistic, dataInfo=dataInfo, univariateModels = univariateModels, hash = hash, stat_hash = stat_hash, pvalue_hash = pvalue_hash, robust=robust);
        stat_hash = cur_results$stat_hash;
        pvalue_hash = cur_results$pvalue_hash;
        
        #check if the pvalues and stats should be updated
        if (compare_p_values(pvalues[[cvar]] , cur_results$pvalue , stats[[cvar]] , cur_results$stat) )
        {
          pvalues[[cvar]] = cur_results$pvalue;
          stats[[cvar]] = cur_results$stat;
        }
        
        #if there is any subset that make var independent from y,
        #then let's throw away var; moreover, we also look for
        #equivalent variables. Note that we stop after the first 
        #z such that pvalue_{var, y | z} > threshold
        if (cur_results$flag & cur_results$pvalue > threshold)
        {
          remainingVars[[cvar]] = 0;
          selectedVars[[cvar]] = 0;
          queues = identifyTheEquivalent.ma(equal_case, queues, target, dataset, cvar, z, test, threshold, statistic, univariateModels, pvalues, hash, dataInfo, stat_hash, pvalue_hash, robust = robust, ncores = ncores);
          breakFlag = TRUE;
          break;
        }
      }
      
      if( breakFlag )
      { 
        break;
      }else
      {
        klimit = klimit + 1;
      }
    }
  }
  results <- list(pvalues = pvalues, stats = stats, queues = queues, selectedVars = selectedVars, remainingVars = remainingVars, stat_hash = stat_hash, pvalue_hash = pvalue_hash, rob = robust);
  return(results);
}

#########################################################################################################

identifyTheEquivalent.ma = function(equal_case, queues, target, dataset, cvar, z, test, threshold, statistic, univariateModels, pvalues, hash, dataInfo, stat_hash, pvalue_hash, robust = robust, ncores = ncores)
{
  z = t(z);
  
  #case 3
  #if we have more than one equivalent vars in z , we select the first one
  #for loop if we dont use the below lapply function
  #equalsY = 0;
  for(i in 1:ncol(z))
  {
    w = z[,i];
    w = t( t(w) );
    zPrime = c(setdiff(z , w) , cvar);
    cur_results = test(target = target, dataset = dataset, xIndex = w, csIndex = zPrime, statistic = statistic, dataInfo = dataInfo, univariateModels = univariateModels, hash = hash, stat_hash = stat_hash, pvalue_hash = pvalue_hash, robust=robust);
    
    if ( cur_results$flag & (cur_results$pvalue > threshold) )
    {  
      queues[[w]] = as.matrix( c(queues[[w]] , queues[[cvar]]) );
      break;
      #equalsY = equalsY+1;
    }
  }
  #cat("\nEquals Ys in the current z subset = %d",equalsY);
  return(queues);
}


#########################################################################################################

apply_ideq.ma = function(i, queues, target, dataset, cvar, z, test, threshold, statistic, univariateModels, hash, dataInfo, stat_hash, pvalue_hash, robust = robust, ncores = ncores)
{
  w = z[,i];
  w = t(t(w));
  zPrime = c(setdiff(z , w) , cvar);
  
  cur_results = test(target = target, dataset = dataset, xIndex = w, csIndex = zPrime, statistic = statistic, dataInfo = dataInfo, univariateModels = univariateModels, hash = hash, stat_hash = stat_hash, pvalue_hash = pvalue_hash, robust = robust);
  
  if(cur_results$flag & (cur_results$pvalue > threshold))
  {
    queues[[w]] = as.matrix(c(queues[[w]] , queues[[cvar]]));
    return(queues[[w]]);
  }else{
    return(NA);
  }
}


#########################################################################################################

max_min_assoc.ma = function(target, dataset, test, threshold, statistic, max_k, selectedVars, pvalues , stats , remainingVars , univariateModels, selectedVarsOrder, hash, dataInfo, stat_hash, pvalue_hash, robust = robust, ncores = ncores)
{
  #Initialize
  selected_var = -1;
  selected_pvalue = 2;
  selected_stat = 0;
  
  varsToIterate = which(remainingVars==1);
  
  for(cvar in varsToIterate)
  {
    mma_res = min_assoc.ma(target, dataset, test, max_k, cvar, statistic, selectedVars , pvalues , stats , univariateModels , selectedVarsOrder, hash, dataInfo, stat_hash, pvalue_hash, robust = robust, ncores = ncores);
    pvalues = mma_res$pvalues;
    stats = mma_res$stats;
    stat_hash = mma_res$stat_hash;
    pvalue_hash = mma_res$pvalue_hash;
    
    
    if(mma_res$pvalue > threshold)
    {
      remainingVars[[cvar]] = 0;
    }
    
    if(compare_p_values(mma_res$pvalue , selected_pvalue , mma_res$stat , selected_stat))
    {
      selected_var = cvar;
      selected_pvalue = mma_res$pvalue;
      selected_stat = mma_res$stat;
    }
  }
  
  results <- list(selected_var = selected_var , selected_pvalue = selected_pvalue , remainingVars = remainingVars , pvalues = pvalues , stats = stats, stat_hash=stat_hash, pvalue_hash = pvalue_hash, rob = robust);
  return(results); 
}

#########################################################################################################

min_assoc.ma = function(target, dataset, test, max_k, cvar, statistic, selectedVars, pvalues, stats, univariateModels , selectedVarsOrder, hash, dataInfo, stat_hash, pvalue_hash, robust = robust, ncores = ncores)
{
  #initialization
  #baseline values
  #   ma_pvalue = univariateModels$pvalue[[cvar]];
  #   ma_stat = univariateModels$stat[[cvar]];
  ma_pvalue = pvalues[[cvar]]; #CHANGE
  ma_stat = stats[[cvar]]; #CHANGE
  
  selectedVars = which(selectedVars==1);
  #max size of the condiotioning test
  k = min(c(max_k , length(selectedVars)));
  
  ck = 1;
  while ( ck <= k )
  {
    #lastvar = unique(which(selectedVarsOrder == max(selectedVarsOrder)));
    lastvar = which(selectedVarsOrder == max(selectedVarsOrder))[1]; #CHANGE
    
    tempCS = setdiff(selectedVars, lastvar) #CHANGE
    if( ck == 1 ) #CHANGE
    {
      subsetcsk = as.matrix(lastvar); #CHANGE
    }else{
      subsetcsk = as.matrix( nchoosek(tempCS, ck - 1) )
      numSubsets = dim(subsetcsk)[2]; #CHANGE
      subsetcsk = rbind(subsetcsk, lastvar*rep(1,numSubsets)); #CHANGE
    }


    for(i in 1:ncol(subsetcsk))
    {
      s = subsetcsk[,i];
      s = t(t(s));
      
      cur_results = test(target = target, dataset = dataset, xIndex = cvar, csIndex = s, statistic = statistic, dataInfo = dataInfo, univariateModels = univariateModels, hash = hash, stat_hash = stat_hash, pvalue_hash = pvalue_hash, robust = robust);
      stat_hash = cur_results$stat_hash;
      pvalue_hash = cur_results$pvalue_hash;
      
      #check if the pvalues and stats should be updated
      if(cur_results$flag == 1 & !compare_p_values(cur_results$pvalue, ma_pvalue, cur_results$stat , ma_stat))
      {
        ma_pvalue = cur_results$pvalue;
        pvalues[[cvar]] = cur_results$pvalue;
        
        ma_stat = cur_results$stat;
        stats[[cvar]] = cur_results$stat;
      }
    }
    ck = ck + 1;
  }
  results <- list(pvalue = ma_pvalue , stat = ma_stat , pvalues = pvalues , stats = stats, stat_hash=stat_hash, pvalue_hash  = pvalue_hash, rob = robust);
  return(results);
}



