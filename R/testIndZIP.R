testIndZIP = function(target, dataset, xIndex, csIndex, wei = NULL, dataInfo=NULL, univariateModels=NULL , hash = FALSE, stat_hash=NULL, pvalue_hash=NULL,robust=FALSE) 
{
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
    if (is.null(stat_hash[[key]]) == FALSE) {
      stat = stat_hash[[key]];
      pvalue = pvalue_hash[[key]];
      flag = 1;
      results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
  }
  #if the xIndex is contained in csIndex, x does not bring any new
  #information with respect to cs
  if ( !is.na( match(xIndex, csIndex) ) ) {
    if (hash) { #update hash objects
      stat_hash[[key]] <- 0;     #.set(stat_hash , key , 0)
      pvalue_hash[[key]] <- log(1);     #.set(pvalue_hash , key , 1)
    }
    results <- list(pvalue = log(1), stat = 0, flag = 1 , stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  #check input validity
  if (xIndex < 0 || csIndex < 0) {
    message(paste("error in testIndZIP : wrong input of xIndex or csIndex"))
    results <- list(pvalue = pvalue, stat = stat, flag = flag , stat_hash=stat_hash, pvalue_hash=pvalue_hash);
    return(results);
  }
  
  xIndex = unique(xIndex);
  csIndex = unique(csIndex);
  #extract the data
  x = dataset[ , xIndex];
  cs = dataset[ , csIndex];
  #if x = any of the cs then pvalue = log(1) and flag = 1.
  #That means that the x variable does not add more information to our model due to an exact copy of this in the cs, so it is independent from the target
  if ( length(cs) != 0 ) {
    if ( is.null(dim(cs)[2]) ) { #cs is a vector
      if (any(x != cs) == FALSE) { #if(!any(x == cs) == FALSE)
        if (hash) {#update hash objects
          stat_hash[[key]] <- 0;    #.set(stat_hash , key , 0)
          pvalue_hash[[key]] <- log(1);     #.set(pvalue_hash , key , 1)
        }
        results <- list(pvalue = log(1), stat = 0, flag = 1 , stat_hash=stat_hash, pvalue_hash=pvalue_hash);
        return(results);
      }
    } else { #more than one var
      for (col in 1:dim(cs)[2]) {
        if (any(x != cs[,col]) == FALSE) { #if(!any(x == cs) == FALSE)
          if (hash) {  #update hash objects
            stat_hash[[key]] <- 0;    #.set(stat_hash , key , 0)
            pvalue_hash[[key]] <- log(1);    #.set(pvalue_hash , key , 1)
          }
          results <- list(pvalue = log(1), stat = 0, flag = 1 , stat_hash=stat_hash, pvalue_hash=pvalue_hash);
          return(results);
        }
      }
    }
  }
  #if x or target is constant then there is no point to perform the test
  #if ( Rfast::Var( as.numeric(x) ) == 0 ) {
  #  if (hash) {     #update hash objects
  #    stat_hash[[key]] <- 0;   #.set(stat_hash , key , 0)
  #    pvalue_hash[[key]] <- log(1);    #.set(pvalue_hash , key , 1)
  #  }
  #  results <- list(pvalue = log(1), stat = 0, flag = 1 , stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  #  return(results);
  #}
  #trycatch for dealing with errors
  res <- tryCatch(
    {
      #if the conditioning set (cs) is empty, we use a simplified formula
      if (length(cs) == 0) {
        #compute the relationship between x,target directly
        if ( !is.null(wei) )  lik1 <- zipmle.wei(target, wei)  else  lik1 <- Rfast::zip.mle(target)
        fit2 <- zip.reg(target, x, wei = wei, lgy = NULL) 
        stat <- 2 * abs( lik1 - fit2$loglik )
        dof <- length( fit2$be ) - 1
        pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)   
        flag <- 1;

      } else {
        lgy <- sum( lgamma(target + 1) )
        fit1 <- zip.reg(target, as.data.frame( dataset[, csIndex] ), wei = wei, lgy = lgy) 
        fit2 <- zip.reg(target, as.data.frame( dataset[, c(csIndex, xIndex)] ), wei = wei, lgy = lgy) 
        stat <- 2 * abs( fit1$loglik - fit2$loglik )
        dof <- length( fit2$be ) - length( fit1$be )
        pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)   
        flag <- 1;
      } 
      #last error check
      if ( is.na(pvalue) || is.na(stat) ) {
        pvalue = log(1);
        stat = 0;
        flag = 0;
      } else {
        #update hash objects
        if( hash )  {
          stat_hash[[key]] <- stat;      #.set(stat_hash , key , stat)
          pvalue_hash[[key]] <- pvalue;     #.set(pvalue_hash , key , pvalue)
        }
      }
      #testerrorcaseintrycatch(4);
      results <- list(pvalue = pvalue, stat = stat, flag = flag , stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
      
    },
    error=function(cond) {
      #   message(paste("error in try catch of the testIndZIP test"))
      #   message("Here's the original error message:")
      #   message(cond)
      #   #        #for debug
      #   #        print("\nxIndex = \n");
      #   #        print(xIndex);
      #   #        print("\ncsindex = \n");
      #   #        print(csIndex);
      #   stop();
      #error case
      pvalue = log(1);
      stat = 0;
      flag = 0;
      
      results <- list(pvalue = pvalue, stat = stat, flag = flag , stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    },
    finally={}
  )    
  return(res);
}