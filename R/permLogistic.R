permLogistic = function(target, dataset, xIndex, csIndex, wei = NULL, dataInfo = NULL, univariateModels = NULL, hash = FALSE, stat_hash = NULL, 
                           pvalue_hash = NULL, robust = FALSE, threshold = 0.05, R = 999) {

  csIndex[ which( is.na(csIndex) ) ] = 0
  thres <- threshold * R + 1
  
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
  pvalue = 1;
  stat = 0;
  flag = 0;
  #if the xIndex is contained in csIndex, x does not bring any new
  #information with respect to cs
  if ( !is.na( match(xIndex, csIndex) ) ) {
    if ( hash) {  #update hash objects
      stat_hash[[key]] <- 0;#.set(stat_hash , key , 0)
      pvalue_hash[[key]] <- 1;#.set(pvalue_hash , key , 1)
    }
    results <- list(pvalue = 1, stat = 0, flag = 1, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  #check input validity
  if (xIndex < 0 || csIndex < 0) {
    message(paste("error in testIndLogistic : wrong input of xIndex or csIndex"))
    results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  #xIndex = unique(xIndex);
  #csIndex = unique(csIndex);
  #extract the data
  x = dataset[ , xIndex];
  cs = dataset[ , csIndex];
  if (length(cs) == 0 || is.na(cs) )  cs = NULL;
  #if x = any of the cs then pvalue = 1 and flag = 1.
  #That means that the x variable does not add more information to our model due to an exact copy of this in the cs, so it is independent from the target
  if ( length(cs) != 0 )  {
    if ( is.null(dim(cs)[2]) ) {  #cs is a vector
      if ( any(x != cs) == FALSE) {  #if(!any(x == cs) == FALSE)
        if (hash) {  #update hash objects
          stat_hash[[key]] <- 0;    #.set(stat_hash , key , 0)
          pvalue_hash[[key]] <- 1;    #.set(pvalue_hash , key , 1)
        }
        results <- list(pvalue = 1, stat = 0, flag = 1 , stat_hash=stat_hash, pvalue_hash=pvalue_hash);
        return(results);
      }
    } else { #more than one var
      for ( col in 1:dim(cs)[2] ) {
        if ( any(x != cs[,col]) == FALSE ) {  #if(!any(x == cs) == FALSE)
          if (hash) {    #update hash objects
            stat_hash[[key]] <- 0;   #.set(stat_hash , key , 0)
            pvalue_hash[[key]] <- 1;   #.set(pvalue_hash , key , 1)
          }
          results <- list(pvalue = 1, stat = 0, flag = 1 , stat_hash=stat_hash, pvalue_hash=pvalue_hash);
          return(results);
        }
      }
    }
  }
  
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

    target_type = floor(target_type);
    if (target_type < 1 || target_type > 3) {
      message(paste("error in testIndLogistic : wrong input of target_type"))
      results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }

  #trycatch for dealing with errors
  if (length(cs) == 0) {
    #binomial
    if ( target_type == 1 ) { 
      fit2 = glm(target ~ x, binomial, weights = wei, model = FALSE)
      dev2 <- fit2$deviance
      stat = fit2$null.deviance - dev2
	  if ( stat > 0 ) {
        step <- 0
        j <- 1		
        n <- length(target)
        while (j <= R & step < thres ) {
          xb <- sample(x, n)  
          bit2 <- glm(target ~ xb, binomial, weights = wei)  
          step <- step + ( bit2$deviance < dev2 )
          j <- j + 1
        }
        pvalue <- (step + 1) / (R + 1)  
	  }	
    } else if ( target_type == 2 ) {  
      # Fitting multinomial Logistic regression
      fit1 <- nnet::multinom(target ~ 1, trace = FALSE, weights = wei)
      fit2 <- nnet::multinom(target ~ x, trace = FALSE, weights = wei)
      dev2 <- fit2$deviance
      stat <- fit1$deviance - dev2
	  if ( stat > 0 ) {  
        step <- 0
        j <- 1		
        n <- length(target)
        while (j <= R & step < thres ) {
          xb <- sample(x, n)  
          bit2 <- nnet::multinom(target, xb, trace = FALSE, weights = wei)  
          step <- step + ( bit2$deviance < dev2 )
          j <- j + 1
        }
        pvalue <- (step + 1) / (R + 1)
	  }	
    }
    else if ( target_type == 3 ) { # ordinal
      #Fitting ordinal Logistic regression
      fit1 <- ordinal::clm(target ~ 1, weights = wei )
      fit2 <- ordinal::clm(target ~ x, weights = wei )
      lik2 <- as.numeric( fit2$logLik )
      stat <- 2 * lik2 - 2 * fit1$logLik
 	  if ( stat > 0 ) {
        step <- 0
        j <- 1		
        n <- length(target)
        while (j <= R & step < thres ) {
          xb <- sample(x, n)  
          bit2 <- ordinal::clm.fit(target ~ xb, weights = wei)            
          step <- step + ( bit2$logLik > lik2 )
          j <- j + 1
        }
        pvalue <- (step + 1) / (R + 1)
	  }	
    }
  } else {
    if ( target_type == 1 ) { # binomial
      fit1 <- glm(target ~ cs, binomial, weights = wei, model = FALSE)
      fit2 <- glm(target ~ cs + x, binomial, weights = wei, model = FALSE)
      dev2 <- fit2$deviance  
      stat <- fit1$deviance - dev2
	  if (stat > 0) {
        step <- 0
        j <- 1		
        n <- length(target)
        while (j <= R & step < thres ) {
          xb <- sample(x, n)  
          bit2 <- glm(target ~ cs + xb, binomial, weights = wei, model = FALSE)  
          step <- step + ( bit2$deviance < dev2 )
          j <- j + 1
        }
        pvalue <- (step + 1) / (R + 1) 
	  }  
    }
    else if ( target_type == 2 ) { # nominal
      #Fitting multinomial Logistic regression
      fit1 <- nnet::multinom( target ~ cs, trace = FALSE, weights = wei)
      fit2 <- nnet::multinom(target ~ cs + x, trace = FALSE, weights = wei)
      dev2 <- deviance(fit2)
      stat <- deviance(fit1) - dev2
 	  if (stat > 0) {
        step <- 0
        j <- 1		
        n <- length(x)
        while (j <= R & step < thres ) {
          xb <- sample(x, n)  
          bit2 <- nnet::multinom(target ~ cs + xb, trace = FALSE, weights = wei)  
          step <- step + ( bit2$deviance < dev2 )
          j <- j + 1
        }
        pvalue <- (step + 1) / (R + 1) 
      }   
    }
    else if ( target_type == 3 ) { # ordinal
      #Fitting ordinal Logistic regression
      fit1 = ordinal::clm( target ~ cs, weights = wei )
      fit2 <- ordinal::clm(target ~ cs + x, weights = wei )
      lik2 <- 2 * fit2$logLik 
      stat <- lik2 - 2 * fit1$logLik
  	  if (stat > 0) {
        step <- 0
        j <- 1		
        n <- length(x)
        while (j <= R & step < thres ) {
          xb <- sample(x, n)  
          bit2 <- ordinal::clm.fit(target ~ cs + xb, weights = wei)            
          step <- step + ( bit2$logLik > lik2 )
          j <- j + 1
        }
        pvalue <- (step + 1) / (R + 1) 
	  } 
    }
  }
  flag = 1;
  #update hash objects
  if (hash) {
    stat_hash[[key]] <- stat;      #.set(stat_hash , key , stat)
    pvalue_hash[[key]] <- pvalue;     #.set(pvalue_hash , key , pvalue)
  }
  #last error check
  if ( is.na(pvalue) || is.na(stat) ) {
    pvalue = 1;
    stat = 0;
    flag = 0;
  } else {
    #update hash objects
    if (hash) {
      stat_hash[[key]] <- stat;      #.set(stat_hash , key , stat)
      pvalue_hash[[key]] <- pvalue;      #.set(pvalue_hash , key , pvalue)
    }
  }
  results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  return(results);
  
}