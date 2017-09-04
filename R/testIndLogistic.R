testIndLogistic = function(target, dataset, xIndex, csIndex, wei = NULL, dataInfo = NULL, univariateModels = NULL, hash = FALSE, stat_hash = NULL, pvalue_hash = NULL, target_type = 0, robust = FALSE)
{

  csIndex[ which( is.na(csIndex) ) ] = 0
  if (hash) {
    csIndex2 = csIndex[which(csIndex!=0)]
    csIndex2 = sort(csIndex2)
    xcs = c(xIndex,csIndex2)
    key = paste(as.character(xcs) , collapse=" ");
    if ( !is.null(stat_hash[[key]]) ) {
      stat = stat_hash[[key]];
      pvalue = pvalue_hash[[key]];
      flag = 1;
      results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
  }
  #if the test cannot performed succesfully these are the returned values
  pvalue = log(1);
  stat = 0;
  flag = 0;
  #if the xIndex is contained in csIndex, x does not bring any new
  #information with respect to cs
  if ( !is.na( match(xIndex, csIndex) ) ) {
    if ( hash) {  #update hash objects
      stat_hash[[key]] <- 0;#.set(stat_hash , key , 0)
      pvalue_hash[[key]] <- log(1);#.set(pvalue_hash , key , 1)
    }
    results <- list(pvalue = log(1), stat = 0, flag = 1, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  #check input validity
  if (xIndex < 0 || csIndex < 0) {
    message(paste("error in testIndLogistic : wrong input of xIndex or csIndex"))
    results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  x = dataset[ , xIndex];
  cs = dataset[ , csIndex];
  if(length(cs) == 0 || is.na(cs) )  cs = NULL;
  #if x = any of the cs then pvalue = 1 and flag = 1.
  #That means that the x variable does not add more information to our model due to an exact copy of this in the cs, so it is independent from the target
  if ( length(cs) != 0 )  {
    if ( is.null(dim(cs)[2]) ) {  #cs is a vector
      if ( any(x != cs) == FALSE) {  #if(!any(x == cs) == FALSE)
        if (hash) {  #update hash objects
          stat_hash[[key]] <- 0;    #.set(stat_hash , key , 0)
          pvalue_hash[[key]] <- log(1);    #.set(pvalue_hash , key , 1)
        }
        results <- list(pvalue = log(1), stat = 0, flag = 1 , stat_hash=stat_hash, pvalue_hash=pvalue_hash);
        return(results);
      }
    } else { #more than one var
      for ( col in 1:dim(cs)[2] ) {
        if ( any(x != cs[,col]) == FALSE ) {  #if(!any(x == cs) == FALSE)
          if (hash) {    #update hash objects
            stat_hash[[key]] <- 0;   #.set(stat_hash , key , 0)
            pvalue_hash[[key]] <- log(1);   #.set(pvalue_hash , key , 1)
          }
          results <- list(pvalue = log(1), stat = 0, flag = 1 , stat_hash=stat_hash, pvalue_hash=pvalue_hash);
          return(results);
        }
      }
    }
  }
  
  if ( target_type == 0 )  {
    target_type = dataInfo$target_type;
    if ( dataInfo$target_type == "nominal" )  {     #cat("Multinomial Logistic Regression")
      target_type = 2;
    } else if ( dataInfo$target_type == "ordinal") {     #cat("Ordinal Logistic Regression")
      target_type = 3;
    } else if (dataInfo$target_type == "binary") {     #cat("Binary Logistic Regression")
      target_type = 1;
    } else {       #cat("Multinomial Logistic Regression")
      target_type = 2;     #default value in case of bad definition
    }
  } else {
    target_type = floor(target_type);
    if (target_type < 1 || target_type > 3) {
      message(paste("error in testIndLogistic : wrong input of target_type"))
      results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
  }
  
  if (length(cs) == 0) {
    #if the univariate models have been already compute
    if ( !is.null(univariateModels) ) {
      pvalue = univariateModels$pvalue[[xIndex]];
      stat = univariateModels$stat[[xIndex]];
      flag = univariateModels$flag[[xIndex]];
      results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
    
    if ( target_type == 1 ) { # binomial
        fit2 = glm(target ~ x, binomial, weights = wei, model = FALSE)
        stat = fit2$null.deviance - fit2$deviance
        dof = length( coef(fit2) ) - 1
    } else if ( target_type == 2 ) {  # nominal
      #Fitting multinomial Logistic regression
      fit1 <- nnet::multinom(target ~ 1, trace = FALSE, weights = wei)
      fit2 <- nnet::multinom(target ~ x, trace = FALSE, weights = wei)
      stat <- fit1$deviance - fit2$deviance
      dof <- length( coef(fit2) ) - length( coef(fit1) )
    }
    else if ( target_type == 3 ) { # ordinal
      #Fitting ordinal Logistic regression
      fit1 <- ordinal::clm(target ~ 1, weights = wei )
      fit2 <- ordinal::clm(target ~ x, weights = wei )
      stat <- 2 * fit2$logLik - 2 * fit1$logLik
      dof <- length( coef(fit2) ) - length( coef(fit1) )
    }
  } else {
    if ( target_type == 1 ) { # binomial
      fit1 <- glm(target ~., data = as.data.frame( dataset[, csIndex] ), binomial, weights = wei, model = FALSE)
      fit2 <- glm(target ~., data = as.data.frame( dataset[, c(csIndex, xIndex)] ), binomial, weights = wei, model = FALSE)
      stat <- fit1$deviance - fit2$deviance  
      dof <- length( coef(fit2) ) - length( coef(fit1) )
    }
    else if ( target_type == 2 ) { # nominal
      #Fitting multinomial Logistic regression
      fit1 <- nnet::multinom( target ~., data = as.data.frame( dataset[, csIndex] ), trace = FALSE, weights = wei)
      fit2 <- nnet::multinom(target ~.,  data = as.data.frame( dataset[, c(xIndex, csIndex)] ), trace = FALSE, weights = wei)
      stat <- deviance(fit1) - deviance(fit2)
      dof <- length( coef(fit2) ) - length( coef(fit1) )
    }
    else if ( target_type == 3 ) { # ordinal
      #Fitting ordinal Logistic regression
      fit1 = ordinal::clm( target ~., data = as.data.frame( dataset[, csIndex] ), weights = wei )
      fit2 <- ordinal::clm(target ~., data = as.data.frame( dataset[, c(csIndex, xIndex)] ), weights = wei )
      stat <- 2 * fit2$logLik - 2 * fit1$logLik
      dof <- length( coef(fit2) ) - length( coef(fit1) )
    }
  }
  
  pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE); 
  flag = 1;
  #update hash objects
  if (hash) {
    stat_hash[[key]] <- stat;      #.set(stat_hash , key , stat)
    pvalue_hash[[key]] <- pvalue;     #.set(pvalue_hash , key , pvalue)
  }
  #last error check
  if ( is.na(pvalue) || is.na(stat) ) {
    pvalue = log(1);
    stat = 0;
    flag = 0;
  } else {
    #update hash objects
    if (hash) {
      stat_hash[[key]] <- stat;      #.set(stat_hash , key , stat)
      pvalue_hash[[key]] <- pvalue;      #.set(pvalue_hash , key , pvalue)
    }
  }
  #testerrorcaseintrycatch(4);
  results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  return(results);
}