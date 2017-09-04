waldSpeedglm = function(target, dataset, xIndex, csIndex, wei = NULL, dataInfo = NULL, univariateModels = NULL, hash = FALSE, 
                           stat_hash = NULL, pvalue_hash = NULL, target_type = 1, robust = FALSE) {
  #   TESTINDSPEEDGLM Conditional independence test for large sample sized data (tens and hundreds of thousands) for normal, binary discrete or ordinal class variables
  #   provides a p-value PVALUE for the null hypothesis: X independent by target
  #   given CS. The pvalue is calculated by comparing a logistic model based 
  #   on the conditioning set CS against a model containing both X and CS. 
  #   The comparison is performed through a chi-square test with one degree 
  #   of freedom on the difference between the deviances of the two models. 
  #   TESTINDSPEEDGLM requires the following inputs:
  #   target: a vector containing the values of the target variable. 
  #   target must be a vector with percentages, binay data, numerical values or integers  
  #   dataset: a numeric data matrix containing the variables for performing
  #   the conditional independence test. They can be mixed variables, either continous or categorical
  #   xIndex: the index of the variable whose association with the target
  #   must be tested. Can be any type of variable, either continous or categorical.
  #   csIndex: the indices of the variables to condition on. They can be mixed variables, either continous or categorical
  #   target_Type: the type of the target
  #   target_type == 1 (normal target)
  #   target_type == 2 (binary target)
  #   target_type == 3 (discrete target)
  #   default target_type=0
  #   this method returns: the pvalue PVALUE, the statistic STAT and a control variable FLAG.
  #   if FLAG == 1 then the test was performed succesfully 
  target = as.numeric( as.vector(target) );
  csIndex[which(is.na(csIndex))] = 0
  
  if ( hash )  {
    csIndex2 = csIndex[which(csIndex!=0)]
    csIndex2 = sort(csIndex2)
    xcs = c(xIndex,csIndex2)
    key = paste(as.character(xcs) , collapse=" ");
    if (is.null(stat_hash[[key]]) == FALSE)   {
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
    if ( hash )  {         #update hash objects
      stat_hash$key <- 0;       #.set(stat_hash , key , 0)
      pvalue_hash$key <- log(1);        #.set(pvalue_hash , key , 1)
    }
    results <- list(pvalue = log(1), stat = 0, flag = 1, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  #check input validity
  if (xIndex < 0 || csIndex < 0)  {
    message(paste("error in testIndSpeedglm : wrong input of xIndex or csIndex"))
    results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  x = dataset[ , xIndex];
  cs = dataset[ , csIndex];
  if(length(cs) == 0 || is.na(cs) )  cs = NULL;
  #if x = any of the cs then pvalue = log(1) and flag = 1.
  #That means that the x variable does not add more information to our model due to an exact copy of this in the cs, so it is independent from the target
  if ( length(cs) != 0 )  {
    if ( is.null(dim(cs)[2]) )  {  #cs is a vector
      if ( any(x != cs) == FALSE )  {    #if(!any(x == cs) == FALSE)
        if( hash )  {      #update hash objects
          stat_hash$key <- 0;     #.set(stat_hash , key , 0)
          pvalue_hash$key <- log(1);       #.set(pvalue_hash , key , 1)
        }
        results <- list(pvalue = log(1), stat = 0, flag = 1 , stat_hash=stat_hash, pvalue_hash=pvalue_hash);
        return(results);
      }
    } else { #more than one var
      for (col in 1:dim(cs)[2]) {
        if (any(x != cs[,col]) == FALSE)  {    #if(!any(x == cs) == FALSE)
          if ( hash )  {     #update hash objects
            stat_hash$key <- 0;      #.set(stat_hash , key , 0)
            pvalue_hash$key <- log(1);    #.set(pvalue_hash , key , 1)
          }
          results <- list(pvalue = log(1), stat = 0, flag = 1 , stat_hash=stat_hash, pvalue_hash=pvalue_hash);
          return(results);
        }
      }
    }
  }
      #binomial or multinomial target?
      yCounts = length( unique(target) );
      if (yCounts == 2) {
        target_type = 1;
      } else  target_type = 2

      if ( length(cs) == 0 ) {
       if ( target_type == 1 ) {
          fit = speedglm::speedglm(target ~ x, weights = wei, data = as.data.frame(x), family = binomial(logit) )
		  if ( any( is.na( fit$coefficients ) ) ) {
		    stat <- 0
			pvalue <- log(1)
            flag <- 1  			
		  } else {
            stat = summary(fit)[[17]][2,3]^2
            pvalue = pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
            flag = 1;
	      }		
          
        } else {
          fit = speedglm::speedglm(target ~ x, weights = wei, data = as.data.frame(x), family = poisson(log) )
		  if ( any( is.na( fit$coefficients ) ) ) {
		    stat <- 0
			pvalue <- log(1)
            flag <- 1  			
		  } else {
            stat = summary(fit)[[17]][2,3]^2
            pvalue = pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
            flag = 1;
		  }	
        }
        
      } else {
        
        if ( target_type == 1 ) {
          fit = speedglm::speedglm( target ~ ., data = as.data.frame( dataset[, c(xIndex, csIndex)] ), family = binomial(logit), weights = wei )
		  if ( any( is.na( fit$coefficients ) ) ) {
		    stat <- 0
			pvalue <- log(1)
            flag <- 1  			
		  } else {
            stat = summary(fit)[[17]][2, 3]^2
            pvalue = pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
            flag = 1;
          }
        } else {
          fit = speedglm::speedglm( target ~ ., data = as.data.frame( dataset[, c(xIndex, csIndex)] ), family = poisson(log), weights = wei )
		  if ( any( is.na( fit$coefficients ) ) ) {
		    stat <- 0
			pvalue <- log(1)
            flag <- 1  			
		  } else {
            stat = summary(fit)[[17]][2, 3]^2
            pvalue = pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
            flag = 1;
		  }	 
        }   
      }
      #update hash objects
      if( hash )  {
        stat_hash$key <- stat;       #.set(stat_hash , key , stat)
        pvalue_hash$key <- pvalue;        #.set(pvalue_hash , key , pvalue)
      }
      
      if ( is.na(pvalue) || is.na(stat) )   {
        pvalue = log(1);
        stat = 0;
        flag = 0;
      } else {
        #update hash objects
        if( hash )  {
          stat_hash[[key]] <- stat;       #.set(stat_hash , key , stat)
          pvalue_hash[[key]] <- pvalue;            #.set(pvalue_hash , key , pvalue)
        }
      }
      results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
}