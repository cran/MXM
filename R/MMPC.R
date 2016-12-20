# R Implementation of MMPC algorithm
# as described in the paper "The max-min hill-climbing Bayesian network structure learning algorithm". Machine learning, 65(1), 31-78, 2006
# by Tsamardinos, Brown and Aliferis
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
#   "testIndFisher" : Fisher conditional independence test for continous targets (or proportions) and continuous predictors only
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
# hash : a boolean variable whuch indicates whether (TRUE) or not (FALSE) to use the hash-based implementation of the statistics of MMPC.
# hashObject : a List with the hash objects (hash package) on which we have cached the generated statistics. 
#              MMPC requires this Object for the hash-based implementation of the statistics. This hashObject is produced or 
# updated by each run
#              of MMPC (if hash == TRUE) and it can be reused in next runs of MMPC.
# robust : Should the Fisher or the normal regression be robustly estimated? TRUE or FALSE 

# there are default values for all of the parameters of the algorithm.

# OUTPUT <LIST>
# The output of the algorithm is a LIST with the following quantities (14) :

# selectedVars : the selected variables i.e. the dependent of the target variables.
# selectedVarsOrder : the increasing order of the selected variables due to their pvalues
# hashObject : the hashObject with the cached statistic results of the current run.

# pvalues : the pvalues of all of the variables.
# stats : the stats of all of the variables.

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

# 
# 
# #hashObject
# library(hash)

MMPC = function(target, dataset, max_k = 3, threshold = 0.05, test = NULL, ini = NULL, wei = NULL, user_test = NULL, hash=FALSE, hashObject=NULL, robust = FALSE, ncores = 1, backward = FALSE)
{
  #get the log threshold
  threshold = log(threshold)
  
  ##############################
  # initialization part of MMPC #
  ##############################

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
      }else{
        stop('hashObject must be a list of two hash objects (stat_hash, pvalue_hash)')
      }
    }else{
      cat('The hash version of MMPC requires the hash package');
      return(NULL);
    }
  }
  
  dataInfo = NULL;
  
  ###################################
  # dataset checking and initialize #
  ###################################
  
  if(!is.null(dataset))
  {
    if( sum( class(target) == "matrix") == 1 )
    {
      if( sum( class(target) == "Surv") == 1 )  stop('Invalid dataset class. For survival analysis provide a dataframe-class dataset');      
    }   
    #check if dataset is an ExpressionSet object of Biobase package
    # if(class(dataset) == "ExpressionSet") {
      #get the elements (numeric matrix) of the current ExpressionSet object.
    #  dataset = Biobase::exprs(dataset);
    #  dataset = t(dataset);#take the features as columns and the samples as rows
      #     } else if(is.data.frame(dataset)) {
      #       if(class(target) != "Surv")
      #       {
      #         dataset = as.matrix(dataset);
      #       }
   # } else if((class(dataset) != "matrix") & (is.data.frame(dataset) == FALSE) ) {
   #   stop('Invalid dataset class. It must be either a matrix, a dataframe or an ExpressionSet');
   # }	
  }
  
  if( is.null(dataset) || is.null(target) ) #|| (dim(as.matrix(target))[2] != 1 & class(target) != "Surv" ))
  {
    stop('invalid dataset or target (class feature) arguments.');
  }else  target = target;
  
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
  
  
  if( sum( class(target) == "matrix") == 1 )
  {
    if(ncol(target) >= 2 & class(target) != "Surv")
    {
      if( (is.null(test) || test == "auto") & (is.null(user_test)) )
      {
        test = "testIndMVreg"
        warning("Multivariate target (ncol(target) >= 2) requires a multivariate test of conditional independence. The testIndMVreg was used. For a user-defined multivariate test, please provide one in the user_test argument.");
      }
    }
    
  }
  
  ################################
  # test checking and initialize #
  ################################
  
  la <- length( Rfast::sort_unique( as.numeric(target) ) )
  
  if(typeof(user_test) == "closure")
  {
    test = user_test;
  }else{
    #auto detect independence test in case of not defined by the user and the test is null or auto
    if(is.null(test) || test == "auto")
    {
      
      if ( la == 2 )   target <- as.factor(target)
      if ( sum( class(target) == "matrix") == 1 )  test = "testIndMVreg"
      
      #if target is a factor then use the Logistic test
      if("factor" %in% class(target))
      {
        test = "testIndLogistic";
        if( is.ordered(target) )
        {
          dataInfo$target_type = "ordinal";
          cat('\nTarget variable type: Ordinal')
        }else{
          if( la == 2 )
          {
            dataInfo$target_type = "binary"
            cat('\nTarget variable type: Binomial')
          }else{
            dataInfo$target_type = "nominal"
            cat('\nTarget variable type: Nominal')
          }
        }
        
      }else if ( class(target) == "numeric" || class(target) == "integer" ){
        
        if( sum( floor(target) - target ) == 0  &  la > 2 )
        {
          test = "testIndPois";
        }else{
          if(class(dataset) == "matrix")
          {
            test = "testIndFisher";
          }
          else if(class(dataset) == "data.frame")
          {
            if ( length( Rfast::which_isFactor(dataset)  ) > 0  )
            {
              test = "testIndReg";
            }else   test = "testIndFisher";
          }
        }
        
      }else if ( survival::is.Surv(target) ){
        test = "censIndCR";
      }else   stop('Target must be a factor, vector, matrix with at least 2 columns column or a Surv object');
    }
    
    if(test == "testIndLogistic")
    {
      if( is.ordered(target) )
      {
        dataInfo$target_type = "ordinal";
        cat('\nTarget variable type: Ordinal')
        
        if(requireNamespace("ordinal", quietly = TRUE, warn.conflicts = FALSE)==FALSE)
        {
          cat("The testIndLogistic test requires the ordinal package for the ordered logistic regression method. Please install it.");
          return(NULL);
        }
        
      }else{
        if( la == 2 )
        {
          dataInfo$target_type = "binary"
          cat('\nTarget variable type: Binomial')
        }else{
          dataInfo$target_type = "nominal"
          cat('\nTarget variable type: Nominal')
          
          if(requireNamespace("nnet", quietly = TRUE, warn.conflicts = FALSE)==FALSE)
          {
            cat("The testIndLogistic test requires the nnet package for the multinomial logistic regression method. Please install it.");
            return(NULL);
          }
        }
      }
    }
    
    #cat("\nConditional independence test used: ");cat(test);cat("\n");
    
    #available conditional independence tests
    av_tests = c("testIndFisher", "testIndSpearman", "testIndReg", "testIndRQ", "testIndBeta", "censIndCR", "censIndWR", "censIndER", "testIndClogit", 
                 "testIndLogistic", "testIndPois", "testIndNB", "testIndBinom", "gSquare", "auto" , "testIndZIP" , "testIndSpeedglm", "testIndMVreg", 
                 "testIndIGreg", NULL);
    
    ci_test = test
    #cat(test)
    
    if(length(test) == 1) #avoid vectors, matrices etc
    {
      test = match.arg(test , av_tests ,TRUE);
      #convert to closure type
      if(test == "testIndFisher")
      {
        #an einai posostiaio target
        if ( min(target) > 0  &  max(target) < 1 ) {
          target = log( target/(1 - target) ) ## logistic normal 
        }       
        test = testIndFisher;
      }
      else if(test == "testIndSpearman")
      {
        #an einai posostiaio target
        if ( min(target) > 0 & max(target) < 1 ) {
          target = log( target / (1 - target) ) ## logistic normal 
        }
        target <- rank(target)
        dataset <- apply(dataset, 2, rank)  
        test <- testIndSpearman;  ## Spearman is Pearson on the ranks of the data
      }
      else if (test == "testIndReg") ## It uMMPC the F test
      {
        #an einai posostiaio target
        if ( min(target) > 0 & max(target) < 1 ) {
          target = log(target/(1-target)) ## logistic normal 
        }     
        test = testIndReg;
      }
      else if(test == "testIndMVreg")
      {
        if ( min(target) > 0  &  sd( Rfast::rowsums(target) ) == 0 )  target = log( target[, -1]/target[, 1] ) 
        test = testIndMVreg;
      }     
      else if(test == "testIndBeta") 
      {
        test = testIndBeta;
      }
      else if(test == "testIndRQ") ## quantile regression
      {
        #an einai posostiaio target
        if ( all( target>0 & target<1 ) ){
          target = log( target/(1 - target) ) ## logistic normal 
        }      
        test = testIndRQ;
        if(requireNamespace("quantreg", quietly = TRUE, warn.conflicts = FALSE)==FALSE)
        {
          cat("The testIndRQ requires the quantreg package. Please install it.");
          return(NULL);
        }
      }
      else if (test == "testIndIGreg") ## Inverse Gaussian regression
      {
        test = testIndIGreg;
      }
      else if (test == "testIndPois") ## Poisson regression
      {
        test = testIndPois;
      }
      else if (test == "testIndSpeedglm") ## Poisson regression
      {
        test = testIndSpeedglm;
      }
      else if (test == "testIndNB") ## Negative binomial regression
      {       
        test = testIndNB;
      }
      else if (test == "testIndZIP") ## Zero inflated Poisson regression
      {
        test = testIndZIP;
      }
      else if(test == "censIndCR")
      {
        test = censIndCR;
        if(requireNamespace("survival", quietly = TRUE, warn.conflicts = FALSE)==FALSE)
        {
          cat("The censIndCR requires the survival package. Please install it.");
          return(NULL);
        }
      }
      else if(test == "censIndWR")
      {
        test = censIndWR;
        if(requireNamespace("survival", quietly = TRUE, warn.conflicts = FALSE)==FALSE)
        {
          cat("The censIndWR requires the survival package. Please install it.");
          return(NULL);
        }
      }
	   else if(test == "censIndER")
      {
        test = censIndER;
        if(requireNamespace("survival", quietly = TRUE, warn.conflicts = FALSE)==FALSE)
        {
          cat("The censIndER requires the survival package. Please install it.");
          return(NULL);
        }
      }
      else if(test == "testIndClogit")
      {
        test = testIndClogit;
        if(requireNamespace("survival", quietly = TRUE, warn.conflicts = FALSE)==FALSE)
        {
          cat("The testIndClogit requires the survival package. Please install it.");
          return(NULL);
        }
      }
      else if(test == "testIndBinom")
      {
        test = testIndBinom;
      }
      else if(test == "testIndLogistic")
      {
        test = testIndLogistic;
      }
      else if(test == "gSquare")
      {
        test = gSquare;
      }
      #more tests here
    }else{
      stop('invalid test option');
    }
  }
  
  ###################################
  # options checking and initialize #
  ###################################
  
  #extracting the parameters
  max_k = floor(max_k);
  varsize = ncol(dataset);
  
  #option checking
  if((typeof(max_k)!="double") || max_k < 1)
  {
    stop('invalid max_k option');
  }
  if(max_k > varsize)
  {
    max_k = varsize;
  }
  if ( (typeof(threshold) != "double") || exp(threshold) == 0 || exp(threshold) > 1 )
  {
    stop('invalid threshold option');
  }
  # if(typeof(equal_case)!="double")
  # {
  #   stop('invalid equal_case option');
  # }
  
  #######################################################################################
  
  if(!is.null(user_test))
  {
    ci_test = "user_test";
  }
  
  #call the main MMPC function after the checks and the initializations
  results = InternalMMPC(target, dataset, max_k, threshold, test, ini, wei, user_test, dataInfo, hash, varsize, stat_hash, pvalue_hash, targetID, robust = robust, ncores = ncores);
  
  #for testing backward phase
  #   results$selectedVars = c(results$selectedVars,15)
  #   results$selectedVarsOrder = c(results$selectedVarsOrder,15)
  #   print(results$selectedVars)
  #   print(results$selectedVarsOrder)
  
  #backword phase
  if ( backward ) {
    
    varsToIterate = results$selectedVars;
    varsOrder = results$selectedVarsOrder;
    met = 1:length(varsToIterate)
    
    if(length(varsToIterate) > 0)
    {
      for(i in 1:length(met))
      {
        tar <- dataset[, varsToIterate[i]];
        datas <- cbind(target, dataset[, -varsToIterate[i]])
        res = InternalMMPC(target = tar, dataset = datas, max_k = max_k, threshold = threshold, ini = NULL, wei = wei, test = test, user_test = user_test, dataInfo = dataInfo, hash = hash, varsize = varsize, stat_hash=stat_hash, pvalue_hash=pvalue_hash, targetID=targetID, robust = robust, ncores = ncores);
        if(1 %in% res$selectedVars == FALSE){
          met[i] = 0;
          #           results$selectedVars = results$selectedVars[-which(results$selectedVars == results$selectedVars[cvar])]
          #           results$selectedVarsOrder = results$selectedVarsOrder[-which(results$selectedVarsOrder == results$selectedVarsOrder[cvar])]
        }
      }
    }
    
    results$selectedVars = varsToIterate[met]
    results$selectedVarsOrder = varsOrder[met]
  }
  

  MMPCoutput <-new("MMPCoutput", selectedVars = results$selectedVars, selectedVarsOrder=results$selectedVarsOrder, hashObject=results$hashObject, pvalues=results$pvalues, stats=results$stats, univ=results$univ, max_k=results$max_k, threshold = results$threshold, runtime=results$runtime, test=ci_test, rob = robust);
  
  return(MMPCoutput);
  
}

#########################################################################################################

InternalMMPC = function(target, dataset, max_k, threshold, test=NULL, ini=NULL, wei=NULL, user_test=NULL, dataInfo=NULL, hash=FALSE, varsize, stat_hash, pvalue_hash, targetID, robust, ncores)
{
  #get the current time
  runtime <- proc.time();
  
  #######################################################################################
  la <- length( Rfast::sort_unique(target) )
  
  rows <- length(target)
  cols <- ncol(dataset)
  ## if testIndSpearman is selected and the user has put robust =TRUE, the robust is not taken into consideration. 
  #univariate feature selection test
 
    if ( is.null(ini) ) { 
    
    univariateModels = list()
    
    if ( identical(test, testIndFisher)  &  !robust )  ## Pearson's correlation 
    {
      a <- as.vector( cor(target, dataset) )
      univariateModels = list();
      dof <- rows - 3; #degrees of freedom
      wa <- 0.5 * log( (1 + a) / (1 - a) ) * sqrt(dof)
      
      if ( targetID != - 1 )  wa[ targetID ] = 0
      
      univariateModels$stat = abs(wa);
      univariateModels$pvalue = log(2) + pt(-abs(wa), dof, log.p = TRUE) ;
      univariateModels$flag = numeric(cols) + 1;
      univariateModels$stat_hash = stat_hash;
      univariateModels$pvalue_hash = pvalue_hash;
      
    } else if ( identical(test, testIndSpearman) ) {  ## Spearman's correlation
      a = as.vector( cor(target, dataset) )
      univariateModels = list();
      dof = rows - 3; #degrees of freedom
      wa = 0.5 * log( (1 + a) / (1 - a) ) * sqrt(dof) / 1.029563
      
      if ( targetID != - 1 )  wa[ targetID ] = 0
      
      univariateModels$stat = abs(wa) 
      univariateModels$pvalue = log(2) + pt(-abs(wa), dof, log.p = TRUE);
      univariateModels$flag = numeric(cols) + 1;
      univariateModels$stat_hash = stat_hash;
      univariateModels$pvalue_hash = pvalue_hash;
      
    } else if ( identical(test, gSquare) ) {  ## Spearman's correlation
      
      z <- cbind(target, dataset)
      dc <- Rfast::colrange(z, cont = FALSE)
      a <- Rfast::g2Test_univariate(z, dc)
      stat <- a$statistic[ a$x == 1 ]
      if ( targetID != - 1 )  stat[ targetID ] = 0
      
      univariateModels$stat = stat
      univariateModels$pvalue = pchisq(stat, a$df[ a$x == 1 ], lower.tail = FALSE, log.p = TRUE)
      univariateModels$flag = numeric(cols) + 1;
      univariateModels$stat_hash = stat_hash;
      univariateModels$pvalue_hash = pvalue_hash;
      
    } else if ( identical(test, testIndBeta) ) {  ## Beta regression
      
      mod <- beta.regs(target, dataset, wei, logged = TRUE, ncores = ncores)
      
      univariateModels$stat = mod[, 1]
      univariateModels$pvalue = mod[, 2]
      univariateModels$flag = numeric(cols) + 1;
      univariateModels$stat_hash = stat_hash;
      univariateModels$pvalue_hash = pvalue_hash;
      
    } else if ( identical(test, testIndReg)  &  robust  ) {  ## M (Robust) linear regression
      
      univariateModels = list();
      fit1 = MASS::rlm(target ~ 1, weights = wei)
      lik1 = as.numeric( logLik(fit1) )
      lik2 = numeric(cols)
      dof = numeric(cols)
      
      if ( ncores <= 1 | is.null(ncores) ) {
        
        for ( i in 1:cols ) {
          
          if ( i != targetID ) {
            
            fit2 = MASS::rlm(target ~ dataset[, i], weights = wei )
            lik2[i] = as.numeric( logLik(fit2) )
            dof[i] = length( coef(fit2) ) - 1
          } else  lik2[i] = lik1
          
        } 
        
        stat = 2 * abs(lik1 - lik2)
        univariateModels$stat = stat
        univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
        univariateModels$flag = numeric(cols) + 1;
        univariateModels$stat_hash = stat_hash;
        univariateModels$pvalue_hash = pvalue_hash;
        
      } else {
        
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mod <- foreach(i = 1:cols, .combine = rbind, .packages = "MASS") %dopar% {
          
          if ( i != targetID ) {
            fit2 = rlm(target ~ dataset[, i], weights = wei )
            lik2 = as.numeric( logLik(fit2) )
            
            return( c(lik2, length( coef(fit2) ) - 1 ) )
          } else{
            return( c(0, 0) )
          }
          
        }
        stopCluster(cl)
        
        lik1 = as.numeric( logLik(fit1) )
        stat = as.vector( 2 * abs(lik1 - mod[, 1]) )
        univariateModels$stat = stat
        univariateModels$pvalue = pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
        univariateModels$flag = numeric(cols) + 1
        univariateModels$stat_hash = NULL
        univariateModels$pvalue_hash = NULL   
      }   
      
    } else if ( identical(test, testIndReg)  &  !robust  &  is.matrix(dataset)  &  is.null(wei) ) {  ## linear regression
        
        mod = Rfast::univglms(target, dataset, oiko = "normal", logged = TRUE) 
        univariateModels$stat = mod[, 1]
        univariateModels$pvalue = mod[, 2]
        univariateModels$flag = numeric(cols) + 1;
        univariateModels$stat_hash = stat_hash;
        univariateModels$pvalue_hash = pvalue_hash;

    } else if ( identical(test, testIndReg)  &  !robust  &  is.data.frame(dataset)  &  is.null(wei) ) {  ## linear regression
        
      mod <- Rfast::regression(dataset, target)
      univariateModels$stat = mod[1, ]
      univariateModels$pvalue = pf(mod[1, ], mod[2, ], rows - mod[2, ], lower.tail = FALSE, log.p = TRUE)
      univariateModels$flag = numeric(cols) + 1;
      univariateModels$stat_hash = stat_hash;
      univariateModels$pvalue_hash = pvalue_hash;

    } else if ( identical(test, testIndSpeedglm)  ) {  ## big glm regresssion
        
        if ( is.factor(target)  ||  la == 2 ) {
          
          target <- as.numeric( as.factor(target) ) - 1 

          if ( targetID == -1  &  is.matrix(dataset)  &  is.null(wei)  ) {
            mod <- Rfast::univglms(target, dataset, oiko = "binomial", logged = TRUE) 
            stat <- mod[, 1]
            pval <- mod[, 2]
            
        } else {
            stat <- dof <- numeric(cols)
            for ( i in 1:cols ) {
              if ( i != targetID ) {
                fit2 <- speedglm::speedglm(target ~., data = data.frame(dataset[, i]), family = binomial(logit) )
                stat[i] <- fit2$deviance
                dof[i] <- length( coef(fit2) ) - 1    
              } else stat[i] <- fit2$nulldev  ;  dof[i] = 1
            }
            stat <- abs( stat - fit2$nulldev )
            pval <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE )
          }  
        
        } else if ( la > 2  &  sum(round(target) - target) == 0 ) {

          if ( targetID == -1  &  is.matrix(dataset)  &  is.null(wei) ) {
            mod <- Rfast::univglms(target, dataset, oiko = "poisson", logged = TRUE) 
            stat <- mod[, 1]
            pval <- mod[, 2]
            
          } else {
            stat <- dof <- numeric(cols)
            for ( i in 1:cols ) {
              if ( i != targetID ) {
                fit2 <- speedglm::speedglm(target ~., data = data.frame(dataset[, i]), family = poisson(log), weights = wei )
                stat[i] <- fit2$deviance
                dof[i] <- length( coef(fit2) ) - 1
              } else stat[i] <- fit2$nulldev   ;  dof[i] = 1
            }
            stat <- abs( stat - fit2$nulldev )
            pval <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE )
          } 
          
        } else {
            
            if ( targetID == -1  &  is.null(wei) ) {
              if ( is.matrix(dataset) )  {
                mod <- Rfast::univglms(target, dataset, oiko = "normal", logged = TRUE) 
                stat <- mod[, 1]
			    pval <- mod[, 2]
              } else if ( is.data.frame(dataset) ) {
                mod <- Rfast::regression(dataset, target) 
                stat <- mod[1, ]
                pval <- pf(stat, mod[2, ], cols - mod[2, ] - 1, lower.tail = FALSE, log.p = TRUE)
              }  
              
            } else {
              stat <- numeric(cols)
              for ( i in 1:cols ) {
                if ( i != targetID ) {
                  fit2 <- speedglm::speedlm(target ~., data = data.frame(dataset[, i]), weights = wei )
                  suma <- summary(fit2)[[ 13 ]]
                  stat[i] <- suma[1]
                  dof[i] <- suma[3]
                } else stat[i] <- 0   ;  dof[i] = 1
              }
              pvalue <- pf(stat, dof, cols - dof, lower.tail = FALSE, log.p = TRUE)
              
            } 
            
        }
  
        univariateModels$stat = stat
        univariateModels$pvalue = pval
        univariateModels$flag = numeric(cols) + 1;
        univariateModels$stat_hash = stat_hash;
        univariateModels$pvalue_hash = pvalue_hash;

    } else if ( identical(test, testIndLogistic)  &  is.ordered(target) ) {  ## ordinal regression
      
      lik2 = numeric(cols)
      dof = numeric(cols)
      univariateModels = list();
      
      fit1 = ordinal::clm(target ~ 1, weights = wei)
      lik1 = as.numeric( logLik(fit1) )
      df1 = length( coef(fit1) )
      
      if ( ncores <= 1 | is.null(ncores) ) {
        
        for ( i in 1:cols ) {
          
          if (i != targetID){
            
            mat <- model.matrix(target ~ dataset[, i] )
            fit2 <- ordinal::clm.fit(target, mat, weights = wei)
            lik2[i] <- as.numeric( fit2$logLik )
            dof[i] <- length( coef(fit2) ) - df1
          } else  lik2[i] <- lik1

        }
        
        stat = as.vector( 2 * abs(lik1 - lik2) )
        
        univariateModels$stat = stat
        univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
        univariateModels$flag = numeric(cols) + 1;
        univariateModels$stat_hash = stat_hash;
        univariateModels$pvalue_hash = pvalue_hash;
        
      } else {
        
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mod <- foreach(i = 1:cols, .combine = rbind, .packages = "ordinal") %dopar% {
          if ( i != targetID ) {
            mat <- model.matrix(target ~ dataset[, i] )
            fit2 <- ordinal::clm.fit(target, mat, weights = wei)
            lik2 <- as.numeric( fit2$logLik )
            
            return( c(lik2, length( coef(fit2) ) ) )
          } else{
            return( c(0, 0) )
          }
        }
        stopCluster(cl)
        
        lik1 = as.numeric( logLik(fit1) )
        stat = as.vector( 2 * abs(lik1 - mod[, 1]) )
        univariateModels$stat = stat
        univariateModels$pvalue = pchisq(stat, mod[, 2] - df1, lower.tail = FALSE, log.p = TRUE)
        univariateModels$flag = numeric(cols) + 1
        univariateModels$stat_hash = NULL
        univariateModels$pvalue_hash = NULL   
      }
      
    } else if ( identical(test, testIndLogistic)  &  la > 2  ) {  ## multinomial regression
      
      target = as.factor( as.numeric( target ) );
      
      lik2 = numeric(cols)
      dof = numeric(cols)
      fit1 = nnet::multinom(target ~ 1, trace = FALSE, weights = wei)
      lik1 = as.numeric( logLik(fit1) )
      df1 = length( coef(fit1) )
      
      if ( ncores <= 1 | is.null(ncores) ) {
        
        for ( i in 1:cols ) {
          
          if (i != targetID){
            
            fit2 = nnet::multinom(target ~ dataset[, i], trace = FALSE, weights = wei )
            lik2[i] = as.numeric( logLik(fit2) )
            dof[i] = length( coef(fit2) ) - df1
          } else  lik2[i] = lik1
        }
        
        stat = as.vector( 2 * abs(lik1 - lik2) )
        
        univariateModels$stat = stat
        univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
        univariateModels$flag = numeric(cols) + 1;
        univariateModels$stat_hash = stat_hash;
        univariateModels$pvalue_hash = pvalue_hash;
        
      } else {
        
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mod <- foreach(i = 1:cols, .combine = rbind, .packages = "nnet") %dopar% {
          if ( i != targetID ) {
            
            
            fit2 = nnet::multinom(target ~ dataset[, i], weights = wei)
            lik2 = as.numeric( logLik(fit2 ) )
            
            return( c(lik2, length( coef(fit2) ) ) )
          } else{
            return( c(0, 0) )
          }
        }
        stopCluster(cl)
        
        lik1 = as.numeric( logLik(fit1) )
        stat = as.vector( 2 * abs(lik1 - mod[, 1]) )
        univariateModels$stat = stat
        univariateModels$pvalue = pchisq(stat, mod[, 2] - df1, lower.tail = FALSE, log.p = TRUE)
        univariateModels$flag = numeric(cols) + 1
        univariateModels$stat_hash = NULL
        univariateModels$pvalue_hash = NULL   
      }
        
    } else if ( identical(test, testIndLogistic)  &  la == 2  &  is.matrix(dataset)  &  is.null(wei) ) {  ## logistic regression
    
        if ( is.factor(target) )      target <- as.numeric(target) - 1

        mod <- Rfast::univglms( target, dataset, oiko = "binomial", logged = TRUE )
        univariateModels$stat = mod[, 1]
        univariateModels$pvalue = mod[, 2]
        univariateModels$flag = numeric(cols) + 1;
        univariateModels$stat_hash = stat_hash;
        univariateModels$pvalue_hash = pvalue_hash;
        
    } else if ( identical(test, testIndLogistic)  &  la == 2  &  ( !is.null(wei)  ||  is.data.frame(dataset)  ) ) {  ## Poisson regression
      
      univariateModels = list();
      fit1 = glm(target ~ 1, binomial, weights = wei)
      lik1 = as.numeric( logLik(fit1) )
      lik2 = numeric(cols)
      dof = numeric(cols)
      
      if ( ncores <= 1 | is.null(ncores) ) {
        
        for ( i in 1:cols ) {
          
          if ( i != targetID ) {
            fit2 = glm( target ~ dataset[, i], binomial, weights = wei )
            lik2[i] = as.numeric( logLik(fit2) )
            dof[i] = length( coef(fit2) ) - 1
          } else {
            lik2[i] = lik1
          }
        }
        
        stat = 2 * abs(lik1 - lik2)
        
        univariateModels$stat = stat
        univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
        univariateModels$flag = numeric(cols) + 1;
        univariateModels$stat_hash = stat_hash;
        univariateModels$pvalue_hash = pvalue_hash;
        
      } else {
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mod <- foreach(i = 1:cols, .combine = rbind) %dopar% {
          
          if ( i != targetID ) {
            fit2 = glm( target ~ dataset[, i], binomial, weights = wei )
            lik2 = as.numeric( logLik(fit2) )
            
            return( c(lik2, length( coef(fit2) ) - 1 ) )
          } else{
            return( c(0, 0) )
          }
          
        }
        stopCluster(cl)
        
        lik1 = as.numeric( logLik(fit1) )
        stat = as.vector( 2 * abs(lik1 - mod[, 1]) )
        univariateModels$stat = stat
        univariateModels$pvalue = pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
        univariateModels$flag = numeric(cols) + 1
        univariateModels$stat_hash = NULL
        univariateModels$pvalue_hash = NULL   
      }
	  
    } else if ( identical(test, testIndPois)  &  is.matrix(dataset)  &  is.null(wei) ) {  ## Poisson regression
      
      mod <- Rfast::univglms( target, dataset, oiko = "poisson", logged = TRUE ) 
      univariateModels$stat = mod[, 1]
      univariateModels$pvalue = mod[, 2]
      univariateModels$flag = numeric(cols) + 1;
      univariateModels$stat_hash = stat_hash;
      univariateModels$pvalue_hash = pvalue_hash;

    } else if ( identical(test, testIndPois)  &  ( !is.null(wei)  ||  is.data.frame(dataset)  ) ) {  ## Poisson regression
      
      univariateModels = list();
      fit1 = glm(target ~ 1, poisson, weights = wei)
      lik1 = as.numeric( logLik(fit1) )
      lik2 = numeric(cols)
      dof = numeric(cols)
      
      if ( ncores <= 1 | is.null(ncores) ) {
        
        for ( i in 1:cols ) {
          
          if ( i != targetID ) {
            fit2 = glm( target ~ dataset[, i], poisson, weights = wei )
            lik2[i] = as.numeric( logLik(fit2) )
            dof[i] = length( coef(fit2) ) - 1
          } else  lik2[i] = lik1
        }
        
        stat = 2 * abs(lik1 - lik2)
        
        univariateModels$stat = stat
        univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
        univariateModels$flag = numeric(cols) + 1;
        univariateModels$stat_hash = stat_hash;
        univariateModels$pvalue_hash = pvalue_hash;
        
      } else {
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mod <- foreach(i = 1:cols, .combine = rbind) %dopar% {
          
          if ( i != targetID ) {
            fit2 = glm( target ~ dataset[, i], poisson, weights = wei )
            lik2 = as.numeric( logLik(fit2) )
            
            return( c(lik2, length( coef(fit2) ) - 1 ) )
          } else{
            return( c(0, 0) )
          }
          
        }
        stopCluster(cl)
        
        lik1 = as.numeric( logLik(fit1) )
        stat = as.vector( 2 * abs(lik1 - mod[, 1]) )
        univariateModels$stat = stat
        univariateModels$pvalue = pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
        univariateModels$flag = numeric(cols) + 1
        univariateModels$stat_hash = NULL
        univariateModels$pvalue_hash = NULL   
      }
      
    } else if ( identical(test, testIndNB) ) {  ## Zero-inflated Poisson regression
      
      univariateModels = list();
      lik1 <- MASS::glm.nb( target ~ 1, weights = wei )$deviance
      if ( ncores <= 1 | is.null(ncores) ) {
        lik2 <- dof <- numeric(cols)
        
        for ( i in 1:cols ) {
          
          if ( i != targetID ) {
            fit2 = MASS::glm.nb( target ~ dataset[, i], weights = wei )
            lik2[i] = fit2$deviance
            dof[i] = length( coef(fit2) ) - 1
          } else  lik2[i] = fit2$null.deviance
        }
        
        stat = abs(lik1 - lik2)
        univariateModels$stat = stat
        univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
        univariateModels$flag = numeric(cols) + 1;
        univariateModels$stat_hash = stat_hash;
        univariateModels$pvalue_hash = pvalue_hash;
        
      } else {
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mod <- foreach(i = 1:cols, .combine = rbind, .packages = "MASS") %dopar% {
          
          if ( i != targetID ) {
            fit2 = MASS::glm.nb( target ~ dataset[, i], weights = wei )
            lik2 = fit2$deviance
            
            return( c(lik2, length( coef(fit2) ) - 1 ) )
          } else{
            return( c(0, 0) )
          }
          
        }
        stopCluster(cl)
        
        stat <- lik1 - as.vector(mod[, 1])
        univariateModels$stat = stat
        univariateModels$pvalue = pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
        univariateModels$flag = numeric(cols) + 1
        univariateModels$stat_hash = NULL
        univariateModels$pvalue_hash = NULL
      }
	  
    } else if ( identical(test, testIndZIP) ) {  ## Zero-inflated Poisson regression
      
      univariateModels = list();
      mod <- zip.regs(target, dataset, wei, logged = TRUE, ncores = ncores) 
      univariateModels$stat = mod[, 1]
      univariateModels$pvalue = mod[, 2]
      univariateModels$flag = numeric(cols) + 1;
      univariateModels$stat_hash = stat_hash;
      univariateModels$pvalue_hash = pvalue_hash;
      
    } else if ( identical(test, testIndRQ) ) {  ## Median (quantile) regression
      
      univariateModels = list();
      fit1 = quantreg::rq(target ~ 1, weights = wei)
      stat = pval = numeric(cols)
      
      if ( ncores <= 1 | is.null(ncores) ) {
        
        for ( i in 1:cols ) {
          
          if (i != targetID) {
            
            fit2 = quantreg::rq(target ~ dataset[, i], weights = wei )
            ww = anova(fit1, fit2, test = "rank")
            df1 = as.numeric( ww[[1]][1] )
            df2 = as.numeric( ww[[1]][2] )
            stat[i] = as.numeric( ww[[1]][3] )
            pval[i] = pf(stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE)
          } else {
            pval[i] = log(1);
            stat[i] = 0;
          }
        }
        
        univariateModels$stat = stat
        univariateModels$pvalue = pval
        univariateModels$flag = numeric(cols) + 1;
        univariateModels$stat_hash = stat_hash;
        univariateModels$pvalue_hash = pvalue_hash;
        
      } else {
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mod <- foreach(i = 1:cols, .combine = rbind, .packages = "quantreg") %dopar% {
          
          if (i != targetID) {
            fit2 = quantreg::rq(target ~ dataset[, i], weights = wei )
            ww = anova(fit1, fit2, test = "rank")
            df1 = as.numeric( ww[[1]][1] )
            df2 = as.numeric( ww[[1]][2] )
            stat = as.numeric( ww[[1]][3] )
            pval = pf(stat, df1, df2, lower.tail = FALSE, log.p = TRUE)
            
            return( c(stat, pval ) )
          } else{
            return( c(0, 0) )
          }
          
        }
        stopCluster(cl)
        
        
        univariateModels$stat = as.vector( mod[, 1] )
        univariateModels$pvalue = as.vector( mod[, 2] )
        univariateModels$flag = numeric(cols) + 1
        univariateModels$stat_hash = NULL
        univariateModels$pvalue_hash = NULL   
      }
    
    } else if ( identical(test, testIndIGreg) ) {  ## Poisson regression
      
      univariateModels = list();
      fit1 = glm(target ~ 1, family = inverse.gaussian(link = log), weights = wei)
      lik1 = as.numeric( logLik(fit1) )
      lik2 = numeric(cols)
      dof = numeric(cols)
      
      if ( ncores <= 1 | is.null(ncores) ) {
        
        for ( i in 1:cols ) {
          
          if ( i != targetID ) {
            fit2 = glm( target ~ dataset[, i], family = inverse.gaussian(link = log), weights = wei )
            lik2[i] = as.numeric( logLik(fit2) )
            dof[i] = length( coef(fit2) ) - 1
          } else {
            lik2[i] = lik1
          }
        }
        
        stat = 2 * abs(lik1 - lik2)
        univariateModels$stat = stat
        univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
        univariateModels$flag = numeric(cols) + 1;
        univariateModels$stat_hash = stat_hash;
        univariateModels$pvalue_hash = pvalue_hash;
        
      } else {
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mod <- foreach(i = 1:cols, .combine = rbind) %dopar% {
          
          if ( i != targetID ) {
            fit2 = glm( target ~ dataset[, i], family = inverse.gaussian(link = log), weights = wei )
            lik2 = as.numeric( logLik(fit2) )
            
            return( c(lik2, length( coef(fit2) ) - 1 ) )
          } else{
            return( c(0, 0) )
          }
          
        }
        stopCluster(cl)
        
        lik1 = as.numeric( logLik(fit1) )
        stat = as.vector( 2 * abs(lik1 - mod[, 1]) )
        univariateModels$stat = stat
        univariateModels$pvalue = pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
        univariateModels$flag = numeric(cols) + 1
        univariateModels$stat_hash = NULL
        univariateModels$pvalue_hash = NULL   
      }  
         
    } else if ( identical(test, censIndCR) ) {  ## Cox regression
      
      univariateModels = list();
      fit1 = survival::coxph(target ~ 1, weights = wei)
      lik1 <- fit1$loglik
      lik2 = numeric(cols)
      dof = numeric(cols)
      
      if ( ncores <= 1 | is.null(ncores) ) {
        
        for ( i in 1:cols ) {
          
          if ( i != targetID ) {
            fit2 = survival::coxph( target ~ dataset[, i], weights = wei )
            lik2[i] = as.numeric( logLik(fit2) )
            dof[i] = length( coef(fit2) )
          } else  lik2[i] = lik1
        }
        
        stat = 2 * abs(lik1 - lik2)
        
        univariateModels$stat = stat
        univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
        univariateModels$flag = numeric(cols) + 1;
        univariateModels$stat_hash = stat_hash;
        univariateModels$pvalue_hash = pvalue_hash;
        
      } else {
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mod <- foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
          
          if ( i != targetID ) {
            fit2 = survival::coxph( target ~ dataset[, i], weights = wei )
            lik2 = as.numeric( logLik(fit2) )
            
            return( c(lik2, length( coef(fit2) ) ) )
          } else{
            return( c(0, 0) )
          }
          
        }
        stopCluster(cl)
        
        lik1 = fit1$loglik
        stat = as.vector( 2 * abs(lik1 - mod[, 1]) )
        univariateModels$stat = stat
        univariateModels$pvalue = pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
        univariateModels$flag = numeric(cols) + 1
        univariateModels$stat_hash = NULL
        univariateModels$pvalue_hash = NULL   
      }
	        
    } else if ( identical(test, censIndWR) ) {  ## Weibull regression
      
      univariateModels = list();
      fit1 = survival::survreg(target ~ 1, weights = wei)
      lik1 = as.numeric( logLik(fit1) )
      lik2 = numeric(cols)
      dof = numeric(cols)
      
      if ( ncores <= 1 | is.null(ncores) ) {
        
        for ( i in 1:cols ) {
          
          if ( i != targetID ) {
            fit2 = survival::survreg( target ~ dataset[, i], weights = wei )
            lik2[i] = as.numeric( logLik(fit2) )
            dof[i] = length( coef(fit2) ) - 1
          } else  lik2[i] = lik1
        }
        
        stat = 2 * abs(lik1 - lik2)
        
        univariateModels$stat = stat
        univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
        univariateModels$flag = numeric(cols) + 1;
        univariateModels$stat_hash = stat_hash;
        univariateModels$pvalue_hash = pvalue_hash;
        
      } else {
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mod <- foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
          
          if ( i != targetID ) {
            fit2 = survival::survreg( target ~ dataset[, i], weights = wei )
            lik2 = as.numeric( logLik(fit2) )
            
            return( c(lik2, length( coef(fit2) ) - 1 ) )
          } else{
            return( c(0, 0) )
          }
          
        }
        stopCluster(cl)
        
        stat = as.vector( 2 * abs(lik1 - mod[, 1]) )
        univariateModels$stat = stat
        univariateModels$pvalue = pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
        univariateModels$flag = numeric(cols) + 1
        univariateModels$stat_hash = NULL
        univariateModels$pvalue_hash = NULL   
      }
	  
    } else{  
      univariateModels = univariateScore(target, dataset, test, wei = wei, dataInfo = dataInfo, hash=hash, stat_hash=stat_hash, pvalue_hash=pvalue_hash, targetID=targetID, robust=robust, ncores=ncores);
    }
    
  } else   univariateModels = ini

   
  pvalues = univariateModels$pvalue;      
  stats = univariateModels$stat;
  flags = univariateModels$flag;
  #   stat_hash = univariateModels$stat_hash;
  #   pvalue_hash = univariateModels$pvalue_hash;
  
  #if we dont have any associations , return
  if ( min(pvalues , na.rm=TRUE) > threshold ) 
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
    results$rob = robust
    
    return(results);
  }
  
  
  #Initialize the data structs
  selectedVars = numeric(varsize);
  selectedVarsOrder = numeric(varsize);

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
  
  ################ main MMPC loop ################
  
  #main MMPC loop
  #loop until there are not remaining vars
  loop = any(as.logical(remainingVars));
  
  #rep = 1;
  while(loop)
  {

    max_min_results = max_min_assoc(target, dataset, test, wei, threshold, max_k, selectedVars, pvalues, stats, remainingVars , univariateModels, selectedVarsOrder, hash=hash, dataInfo, stat_hash=stat_hash, pvalue_hash=pvalue_hash, robust = robust, ncores = ncores);
    selectedVar = max_min_results$selected_var;
    selectedPvalue = max_min_results$selected_pvalue;
    remainingVars = max_min_results$remainingVars;
    pvalues = max_min_results$pvalues;
    stats = max_min_results$stats;
    stat_hash=max_min_results$stat_hash;
    pvalue_hash=max_min_results$pvalue_hash;
    
    #if the selected variable is associated with target , add it to the selected variables
    if ( selectedPvalue <= threshold )
    {
      #print(paste("rep: ",rep,", selected var: ",selectedVar,", pvalue = ",exp(selectedPvalue)))
      #rep = rep + 1;
      
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
  
  # #queues correctness
  # all_queues = queues
  # queues = queues[which(selectedVars==1)];
  # 
  # queues <- lapply(1:length(queues) , function(i){queues[[i]] = unique(queues[[i]]);});
  # 
  #adjusting the results
  if(targetID > 0)
  {
    toAdjust <- which(selectedVars > targetID)
    selectedVars[toAdjust] = selectedVars[toAdjust] + 1
  }
  
  
  results = NULL;
  results$selectedVars = which( selectedVars == 1 )
  
  svorder = sort(pvalues[results$selectedVars], index.return = TRUE)
  svorder = results$selectedVars[svorder$ix]
  results$selectedVarsOrder = svorder
  
  hashObject = NULL
  hashObject$stat_hash = stat_hash
  hashObject$pvalue_hash = pvalue_hash
  results$hashObject = hashObject
  class(results$hashObject) = 'list'
  
  results$pvalues = exp(pvalues)
  results$stats = stats;
  results$univ = univariateModels
  
  #   results$all_queues = all_queues;
  #   already known
  #   results$data = dataset;
  #   results$target = target;
  #   results$test = test;
  results$max_k = max_k
  results$threshold = exp(threshold)
  
  runtime = proc.time() - runtime
  results$runtime = runtime
  results$rob = robust
  
  
  return(results)
  
}



