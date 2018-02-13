censIndCR = function(target, dataset, xIndex, csIndex, wei = NULL, univariateModels=NULL, hash = FALSE, stat_hash=NULL, pvalue_hash=NULL){
  # Conditional independence test based on the Log Likelihood ratio test
  if ( !survival::is.Surv(target) )   stop('The survival test can not be performed without a Surv object target');
  csIndex[ which( is.na(csIndex) ) ] = 0;
  
  if ( hash ) {
    csIndex2 = csIndex[which(csIndex!=0)]
    csindex2 = sort(csIndex2)
    xcs = c(xIndex,csIndex2)
    key = paste(as.character(xcs) , collapse=" ");
    if ( is.null(stat_hash[[key]]) == FALSE ) {
      stat = stat_hash[[key]];
      pvalue = pvalue_hash[[key]];
      results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
  }
  #initialization: these values will be returned whether the test cannot be carried out
  pvalue = log(1);
  stat = 0;
  results <- list(pvalue = pvalue, stat = stat, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  cox_results = NULL;
  cox_results_full = NULL;
  #timeIndex = dim(dataset)[2];
  event = target[,2]
  numCases = dim(dataset)[1];
  if ( length(event) == 0 )  event = vector('numeric',numCases) + 1;
      if (is.na(csIndex) || length(csIndex) == 0 || csIndex == 0) {
	    options(warn = -1)
        cox_results <- survival::coxph(target ~ dataset[, xIndex], weights = wei )
        res <- anova(cox_results)
        dof <- res[2, 3]
        stat <- res[2, 2]
        pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE);
      } else {
        options(warn = -1)
        cox_results <- survival::coxph(target ~ ., data = as.data.frame(  dataset[ , csIndex] ), weights = wei) 
        cox_results_full <- survival::coxph(target ~ ., data = as.data.frame(  dataset[ , c(csIndex, xIndex)] ), weights = wei) 
        res = anova(cox_results_full, cox_results)
        stat = res[2, 2]
        dF = res[2, 3]
        pvalue = pchisq(stat, dF, lower.tail = FALSE, log.p = TRUE)
      }  
      if ( is.na(pvalue) || is.na(stat) ) {
        pvalue = log(1);
        stat = 0;
      } else {
        #update hash objects
        if( hash )  {
          stat_hash[[key]] <- stat;      #.set(stat_hash , key , stat)
          pvalue_hash[[key]] <- pvalue;     #.set(pvalue_hash , key , pvalue)
        }
      }
      #testerrorcaseintrycatch(4);
      results <- list(pvalue = pvalue, stat = stat,  stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
}