testIndTobit = function(target, dataset, xIndex, csIndex, wei = NULL, univariateModels=NULL, hash = FALSE, stat_hash=NULL, pvalue_hash=NULL){
  # Conditional independence test based on the Log Likelihood ratio test
  if (!survival::is.Surv(target) )   stop('The survival test can not be performed without a Surv object target');
  csIndex[which(is.na(csIndex))] = 0;
  if ( hash )  {
    csIndex2 = csIndex[which(csIndex!=0)]
    csindex2 = sort(csIndex2)
    xcs = c(xIndex,csIndex2)
    key = paste(as.character(xcs) , collapse=" ");
    if (is.null(stat_hash[key]) == FALSE) {
      stat = stat_hash[key];
      pvalue = pvalue_hash[key];
      results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
  }
  #initialization: these values will be returned whether the test cannot be carried out
  pvalue = log(1);
  stat = 0;
  results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  tob <- NULL;
  tob_full <- NULL;
  #if the censored indicator is empty, a dummy variable is created
  if  ( length(csIndex) == 0 || sum(csIndex == 0, na.rm = TRUE) > 0 ) {
    tob <- survival::survreg(target ~ ., data = as.data.frame( dataset[, xIndex] ), weights = wei, dist = "gaussian")
    dof <- length( coef(tob) ) - 1
    stat <- 2 * abs( diff(tob$loglik) )
    pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE);
  } else {
    tob <- survival::survreg(target ~ ., data = as.data.frame( dataset[, csIndex] ), weights = wei, dist = "gaussian") 
    tob_full <- survival::survreg(target ~ ., data = as.data.frame(  dataset[, c(csIndex, xIndex)] ), weights = wei, dist = "gaussian")
    res <- anova(tob, tob_full)
    stat <- res[2, 6] ;
    dF <- res[2, 5];
    pvalue = pchisq(stat, dF, lower.tail = FALSE, log.p = TRUE)
  }  
  if ( is.na(pvalue) || is.na(stat) ) {
    pvalue <- log(1);
    stat <- 0;
  } else {
    #update hash objects
    if( hash )  {
      stat_hash[key] <- stat;      #.set(stat_hash , key , stat)
      pvalue_hash[key] <- pvalue;     #.set(pvalue_hash , key , pvalue)
    }
  }
  results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  return(results);
}
