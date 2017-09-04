testIndReg = function(target, dataset, xIndex, csIndex, wei = NULL, dataInfo = NULL, univariateModels=NULL , hash = FALSE, stat_hash = NULL, pvalue_hash = NULL, robust=FALSE) 
  {
  # TESTINDREG Conditional Independence Test for continous class variables 
  # PVALUE = TESTINDREG(Y, DATA, XINDEX, CSINDEX, DATAINFO)
  # This test provides a p-value PVALUE for the NULL hypothesis H0 which is
  # X is independent by TARGET given CS. The pvalue is calculated following
  # nested models
  # This method requires the following inputs
  #   TARGET: a numeric vector containing the values of the target (continuous) variable. 
  #   Its support can be R or any number betweeen 0 and 1, i.e. it contains proportions.
  #   DATASET: a numeric data matrix containing the variables for performing the test. They can be mixed variables. 
  #   XINDEX: the index of the variable whose association with the target we want to test. 
  #   CSINDEX: the indices if the variable to condition on. 
  #   DATAINFO: information on the structure of the data
  # this method returns: the pvalue PVALUE, the statistic STAT and a control variable FLAG.
  # if FLAG == 1 then the test was performed succesfully
  # References
  # [1] Norman R. Draper and Harry Smith. Applied Regression 
  # Analysis, Wiley, New York, USA, third edition, May 1998.
  #initialization
  #if the test cannot performed succesfully these are the returned values
  pvalue = log(1);
  stat = 0;
  flag = 0;
  csIndex[ which( is.na(csIndex) ) ] = 0;
  
  if (hash) {
    csIndex2 = csIndex[which(csIndex!=0)]
    csIndex2 = sort(csIndex2)
    xcs = c(xIndex,csIndex2)
    key = paste(as.character(xcs) , collapse=" ");
    if (is.null(stat_hash[[key]]) == FALSE)  {
      stat = stat_hash[[key]];
      pvalue = pvalue_hash[[key]];
      flag = 1;
      results <- list(pvalue = pvalue, stat = stat, flag = flag , stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
  }
  #if the xIndex is contained in csIndex, x does not bring any new
  #information with respect to cs
  if ( !is.na( match(xIndex, csIndex) ) ) {
    if(hash) {  #update hash objects
      stat_hash[[key]] <- 0;#.set(stat_hash , key , 0)
      pvalue_hash[[key]] <- log(1);#.set(pvalue_hash , key , 1)
    }
    results <- list(pvalue = log(1), stat = 0, flag = 1 , stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  #check input validity
  if(xIndex < 0 || csIndex < 0) {
    message(paste("error in testIndReg : wrong input of xIndex or csIndex"))
    results <- list(pvalue = pvalue, stat = stat, flag = flag , stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  
  xIndex = unique(xIndex);
  csIndex = unique(csIndex);
  #extract the data
  x = dataset[ , xIndex];
  cs = dataset[ , csIndex];
  n = length(target)
  #if x = any of the cs then pvalue = 1 and flag = 1.
  #That means that the x variable does not add more information to our model due to an exact copy of this in the cs, so it is independent from the target
  if ( length(cs) != 0 ) {
    if ( is.null(dim(cs)[2]) ) {   #cs is a vector
      if (any(x != cs) == FALSE) {   #if(!any(x == cs) == FALSE)
        if (hash) {    #update hash objects
          stat_hash[[key]] <- 0;   #.set(stat_hash , key , 0)
          pvalue_hash[[key]] <- log(1);   #.set(pvalue_hash , key , 1)
        }
        results <- list(pvalue = log(1), stat = 0, flag = 1 , stat_hash=stat_hash, pvalue_hash=pvalue_hash);
        return(results);
      }
    } else {     #more than one var
      for ( col in 1:dim(cs)[2] ) {
        if (any(x != cs[,col]) == FALSE) {    #if(!any(x == cs) == FALSE)
          if (hash) {     #update hash objects
            stat_hash[[key]] <- 0;    #.set(stat_hash , key , 0)
            pvalue_hash[[key]] <- log(1);    #.set(pvalue_hash , key , 1)
          }
          results <- list(pvalue = log(1), stat = 0, flag = 1 , stat_hash=stat_hash, pvalue_hash=pvalue_hash);
          return(results);
        }
      }
    }
  }
  res <- tryCatch(
{
  #if the conditioning set (cs) is empty, we use a simplified formula
  if (length(cs) == 0) {
    if (robust) { ## robust estimation

       fit1 = MASS::rlm(target ~ 1, weights = wei, maxit= 2000, method = "MM")
       fit2 = MASS::rlm(target ~ x, weights = wei, maxit = 2000, method = "MM")
       stat = 2 * as.numeric( logLik(fit2) ) - 2 * as.numeric( logLik(fit1) )
       dof = length( coef(fit2) ) - 1
       pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
      
    } else {
      s0 = Rfast::Var(target) * (n - 1)
      mat = model.matrix(target~ ., as.data.frame(x) )
	    if ( !is.null(wei) ) {
	      fit2 = lm.wfit(mat, target, w = wei)
	    } else   fit2 = lm.fit(mat, target)
      s1 = sum( fit2$residuals^2 )
      dof = fit2$df.residual
      stat = ( ( s0 - s1 ) / (n - 1 - dof ) ) / (s1 / dof)
      pvalue = pf(stat, 1, dof, lower.tail = FALSE, log.p = TRUE)
    } 
    
  } else {
     if (robust) { ##robust estimation
       if ( length(csIndex) == 1 ) {
         fit1 = MASS::rlm(target ~ dataset[,  csIndex], maxit = 2000, weights = wei, method = "MM" )
       } else  fit1 = MASS::rlm(target ~ ., data = as.data.frame(dataset[,  csIndex]), maxit = 2000, method = "MM" )
       fit2 = MASS::rlm(target ~ ., data = as.data.frame(dataset[, c(csIndex, xIndex)]), maxit = 2000, method = "MM" )
       stat = 2 * as.numeric( logLik(fit2) ) - 2 * as.numeric( logLik(fit1) )
       dof = length( coef(fit2) ) - length( coef(fit1) ) 
       pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
     } else {
       fit2 = lm( target ~., data = as.data.frame( dataset[, c(csIndex, xIndex)] ), weights = wei )
       if ( any( is.na(fit2$coefficients) ) ) {
         stat = 0
         pvalue = log(1)
       } else {   
         fit1 = lm( target ~., data = as.data.frame( dataset[, csIndex] ), weights = wei )
         mod = anova(fit1, fit2)
         stat = as.numeric(mod[2, 5] )  
         df1 = abs(mod[2, 3])
         df2 = mod[2, 1]
         pvalue = pf(stat, df1, df2, lower.tail = FALSE, log.p = TRUE)  
       } 
     }
   }
  flag = 1;
  #last error check
  if (is.na(pvalue) || is.na(stat) )  {
    pvalue = log(1);
    stat = 0;
    flag = 0;
  } else {
    #update hash objects
    if ( hash )  {
      stat_hash[[key]] <- stat;#.set(stat_hash , key , stat)
      pvalue_hash[[key]] <- pvalue;#.set(pvalue_hash , key , pvalue)
    }
  }
  #testerrorcaseintrycatch(4);
  results <- list(pvalue = pvalue, stat = stat, flag = flag , stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  return(results);
},
error=function(cond) {
  pvalue = log(1);
  stat = 0;
  flag = 1;
  results <- list(pvalue = pvalue, stat = stat, flag = flag , stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  return(results);
},
finally={}
  )    
  return(res);
}