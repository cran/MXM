censIndWR = function(target, dataset, xIndex, csIndex, wei = NULL, dataInfo=NULL, univariateModels=NULL, hash = FALSE, stat_hash=NULL, pvalue_hash=NULL,robust=FALSE){
  # Conditional independence test based on the Log Likelihood ratio test
  if (!survival::is.Surv(target) )   stop('The survival test can not be performed without a Surv object target');
  csIndex[which(is.na(csIndex))] = 0;
  if ( hash )  {
    csIndex2 = csIndex[which(csIndex!=0)]
    csindex2 = sort(csIndex2)
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
  #initialization: these values will be returned whether the test cannot be carried out
  pvalue = log(1);
  stat = 0;
  flag = 0;
  results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  weibull_results = NULL;
  weibull_results_full = NULL;
  #timeIndex = dim(dataset)[2];
  event = target[,2]#dataInfo$event;
  x = dataset[ , xIndex];
  #if the censored indicator is empty, a dummy variable is created
  numCases = dim(dataset)[1];
  if (length(event) == 0)  event = vector('numeric',numCases) + 1;
      if (is.na(csIndex) || length(csIndex) == 0 || csIndex == 0) {
        weibull_results <- survival::survreg(target ~ x, weights = wei)
        dof <- length( coef(weibull_results) ) - 1
        stat = 2 * abs( diff(weibull_results$loglik) )
        pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE);
      } else {
        weibull_results <- survival::survreg(target ~ ., data = as.data.frame( dataset[ , c(csIndex)] ), weights = wei) 
        weibull_results_full <- survival::survreg(target ~ ., data = as.data.frame(  dataset[ , c(csIndex, xIndex)] ), weights = wei )
        res = anova(weibull_results, weibull_results_full)
        stat = abs( res[2, 6] );
        dF = abs( res[2, 5] );
        pvalue = pchisq(stat, dF, lower.tail = FALSE, log.p = TRUE)
      }  
      flag = 1;
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

}
