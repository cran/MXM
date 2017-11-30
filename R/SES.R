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
SES = function(target, dataset, max_k = 3, threshold = 0.05 , test = NULL, ini = NULL, wei = NULL, user_test = NULL, 
               hash = FALSE, hashObject = NULL, robust = FALSE, ncores = 1, logged = FALSE) {
  #get the log threshold
  threshold = log(threshold)
  ##############################
  # initialization part of SES #
  ##############################
   equal_case = 3;
  stat_hash = NULL;
  pvalue_hash = NULL;
  
  if ( hash )  {
    if ( requireNamespace("hash") )  {
      if (is.null(hashObject) )  {
        stat_hash = hash();
        pvalue_hash = hash();
      } else if ( class(hashObject) == "list" ) {
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
  if ( !is.null(dataset) ) {
    if ( sum( class(target) == "matrix") == 1 )  {
      if ( sum( class(target) == "Surv") == 1 )  stop('Invalid dataset class. For survival analysis provide a dataframe-class dataset');      
    }
    #check if dataset is an ExpressionSet object of Biobase package
    #if(class(dataset) == "ExpressionSet") {
      #get the elements (numeric matrix) of the current ExpressionSet object.
    #  dataset = Biobase::exprs(dataset);
    #  dataset = t(dataset);#take the features as columns and the samples as rows
#     }else if(is.data.frame(dataset)){
#       if(class(target) != "Surv")
#       {
#         dataset = as.matrix(dataset);
#       }
    #}else if((class(dataset) != "matrix") & (is.data.frame(dataset) == FALSE) ){
    #  stop('Invalid dataset class. It must be either a matrix, a dataframe or an ExpressionSet');
    #}
  }
  if ( is.null(dataset) || is.null(target) ) {  #|| (dim(as.matrix(target))[2] != 1 & class(target) != "Surv" ))
    stop('invalid dataset or target (class feature) arguments.');
  } else  target = target;
  #check for NA values in the dataset and replace them with the variable median or the mode
  if ( any(is.na(dataset) ) ) {
    #dataset = as.matrix(dataset);
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    if ( is.matrix(dataset) )  {
      dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
    } else {
      poia <- which( is.na(dataset), arr.ind = TRUE )[2]
      for ( i in poia )  {
        xi <- dataset[, i]
        if ( is.numeric(xi) ) {                    
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
  if ( is.character(target) & length(target) == 1 ) {
    findingTarget <- target == colnames(dataset);#findingTarget <- target %in% colnames(dataset);
    if ( !sum(findingTarget) == 1 ) {
      warning('Target name not in colnames or it appears multiple times');
      return(NULL);
    }
    targetID <- which(findingTarget);
    target <- dataset[ , targetID];
  }
  
  #checking if target is a single number
  if (is.numeric(target) & length(target) == 1){
    if (target > dim(dataset)[2]) {
      warning('Target index larger than the number of variables');
      return(NULL);
    }
    targetID <- target;
    target <- dataset[ , targetID];
  }

  if ( is.matrix(target) ) {
    if ( ncol(target) >= 2  &  class(target) != "Surv" ) {
      if ( (is.null(test) || test == "auto")  &  is.null(user_test) ) {
        test = "testIndMVreg"
        warning("Multivariate target (ncol(target) >= 2) requires a multivariate test of conditional independence. The testIndMVreg was used. For a user-defined multivariate test, please provide one in the user_test argument.");
      }
    }
  }
  ################################
  # test checking and initialize #
  ################################
  la <- length( unique( as.numeric(target) ) )
  
  if (typeof(user_test) == "closure") {
    test = user_test;
  } else {
    #auto detect independence test in case of not defined by the user and the test is null or auto
    if ( is.null(test) || test == "auto" ) {
      
      if ( la == 2 )   target <- as.factor(target)
      if ( sum( class(target) == "matrix") == 1 )  test = "testIndMVreg"
      #if target is a factor then use the Logistic test
      if ( "factor" %in% class(target) )  {
        test = "testIndLogistic";
        if( is.ordered(target) ) {
          dataInfo$target_type = "ordinal";
        } else {
          if ( la == 2 ) {
            dataInfo$target_type = "binary"
          } else {
            dataInfo$target_type = "nominal"
          }
        }
        
      } else if ( ( is.numeric(target) || is.integer(target) ) & survival::is.Surv(target) == FALSE ) {
        
        if ( sum( floor(target) - target ) == 0  &  la > 2 ) {
          test = "testIndPois";
        } else {
          if ( is.matrix(dataset) ) { 
            test = "testIndFisher";
          } else if ( is.data.frame(dataset) ) {
            if ( length( Rfast::which_isFactor(dataset)  ) > 0  ) {
              test = "testIndReg";
            } else   test = "testIndFisher";
          }
        }
        
      } else if ( survival::is.Surv(target) ) {
        test = "censIndCR";
      } else   stop('Target must be a factor, vector, matrix with at least 2 columns column or a Surv object');
    }
    
    if ( test == "testIndLogistic") {
      if ( is.ordered(target) )  {
        dataInfo$target_type = "ordinal";

      } else {
        if ( la == 2 )  {
          dataInfo$target_type = "binary"
        } else {
          dataInfo$target_type = "nominal"
        }
      }
    }
    
    #available conditional independence tests
    av_tests = c("testIndFisher", "testIndSpearman", "testIndReg", "testIndRQ", "testIndBeta", "censIndCR", "censIndWR", "censIndER", 
                 "testIndClogit", "testIndLogistic", "testIndPois", "testIndNB", "testIndBinom", "gSquare", "auto", "testIndZIP", 
                 "testIndSpeedglm", "testIndMVreg", "testIndIGreg", "testIndGamma", "testIndNormLog", "testIndTobit", NULL);
    ci_test = test
    #cat(test)
    if (length(test) == 1) {      #avoid vectors, matrices etc
      test = match.arg(test, av_tests, TRUE);
      #convert to closure type
      if (test == "testIndFisher") {
        #an einai posostiaio target
        test = testIndFisher;
      }
      else if (test == "testIndSpearman")  {
        #an einai posostiaio target
        target <- rank(target)
        dataset <- apply(dataset, 2, rank)  
        test <- testIndSpearman;  ## Spearman is Pearson on the ranks of the data
      }
      else if (test == "testIndReg") ## It uMMPC the F test
      {
        #an einai posostiaio target
        test = testIndReg;
      }
      else if (test == "testIndMVreg") {
        if ( min(target) > 0 & sd( Rfast::rowsums(target) ) == 0 )  target = log( target[, -1]/target[, 1] ) 
        test = testIndMVreg;
        
      } else if (test == "testIndBeta") {
        test = testIndBeta;
        
      } else if (test == "testIndRQ") {  ## quantile regression
        #an einai posostiaio target
        test = testIndRQ;
        
      } else if (test == "testIndIGreg") { ## Inverse Gaussian regression
        test = testIndIGreg;
        
      } else if (test == "testIndPois") { ## Poisson regression
        test = testIndPois;
        
      } else if (test == "testIndSpeedglm") {  
        test = testIndSpeedglm;
        
      } else if (test == "testIndNB") { ## Negative binomial regression
        test = testIndNB;
        
      } else if (test == "testIndGamma") {  ## Gamma regression
        test = testIndGamma;
        
      } else if (test == "testIndNormLog") { ## Normal regression with a log link
        test = testIndNormLog;
        
      } else if (test == "testIndZIP") { ## Zero inflated Poisson regression
        test = testIndZIP;
        
      } else if (test == "testIndTobit") { ## Tobit regression
        test = testIndTobit;
        
      } else if(test == "censIndCR") {
        test = censIndCR;
        
      } else if(test == "censIndWR") {
        test = censIndWR;

      } else if(test == "censIndER") {
        test = censIndER;
        
      } else if(test == "testIndClogit") {
        test = testIndClogit;
        
      } else if(test == "testIndBinom") {
        test = testIndBinom;
        
      } else if(test == "testIndLogistic") {
        test = testIndLogistic;
        
      } else if(test == "gSquare") {
        test = gSquare;
      }
      #more tests here
    } else {
      stop('invalid test option');
    }
  }
  ###################################
  # options checking and initialize #
  ###################################
  #extracting the parameters
  max_k = floor(max_k);
  varsize = dim(dataset)[2];
  #option checking
  if ( (typeof(max_k)!="double") || max_k < 1 )  stop('invalid max_k option');
  if ( max_k > varsize )   max_k = varsize;
  if ( (typeof(threshold) != "double" ) || exp(threshold) == 0  || exp(threshold) > 1 )   stop('invalid threshold option');
  if ( typeof(equal_case) != "double" )   stop('invalid equal_case option');

  if ( !is.null(user_test) )   ci_test = "user_test";
  #call the main SES function after the checks and the initializations
  results = InternalSES(target, dataset, max_k, threshold, test, ini, wei, equal_case, user_test, dataInfo, hash, varsize, stat_hash, pvalue_hash, targetID, robust = robust, ncores = ncores, logged = logged);
  SESoutput <- new("SESoutput", selectedVars = results$selectedVars, selectedVarsOrder=results$selectedVarsOrder, queues=results$queues, signatures=results$signatures, hashObject=results$hashObject, pvalues=results$pvalues, stats=results$stats, univ = results$univ, max_k=results$max_k, threshold = results$threshold, n.tests = results$n.tests, runtime=results$runtime, test=ci_test, rob = robust);
  return(SESoutput);
}




# .onAttach <- function(libname, pkgname){
#   # do whatever needs to be done when the package is loaded
#   packageStartupMessage( "Loading MXM package version 0.2, thanks for downloading." )
#   #load the dll files from the fortran code for the package
#   #library.dynam("MXM", pkgname, libname)
# }
