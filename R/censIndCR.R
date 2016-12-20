censIndCR = function(target, dataset, xIndex, csIndex, wei = NULL, dataInfo=NULL, univariateModels=NULL, hash = FALSE, stat_hash=NULL, pvalue_hash=NULL,robust=FALSE){
  # Conditional independence test based on the Log Likelihood ratio test
  
  if(!survival::is.Surv(target))
  {
    stop('The survival test can not be performed without a Surv object target');
  }
  
  csIndex[which(is.na(csIndex))] = 0;
  
  if( hash )
  {
    csIndex2 = csIndex[which(csIndex!=0)]
    csindex2 = sort(csIndex2)
    xcs = c(xIndex,csIndex2)
    key = paste(as.character(xcs) , collapse=" ");
    if(is.null(stat_hash[[key]]) == FALSE)
    {
      stat = stat_hash[[key]];
      pvalue = pvalue_hash[[key]];
      flag = 1;
      
      results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
      return(results);
    }
  }
  
  #initialization: these values will be returned whether the test cannot be carried out
  
  #in this test dataset must be a dataframe
  #dataset = as.data.frame(dataset);
  #dataset = cbind(dataset,target[,1]);#dataset$timeToEvent = target[,1];#dataInfo$timeToEvent;
  
  pvalue = log(1);
  stat = 0;
  flag = 0;
  results <- list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  cox_results = NULL;
  cox_results_full = NULL;
  
  #timeIndex = dim(dataset)[2];
  event = target[,2]#dataInfo$event;
  
  #retrieving the data
  x = dataset[ , xIndex];
  
  #if the censored indicator is empty, a dummy variable is created
  numCases = dim(dataset)[1];
  if (length(event)==0){
    event = vector('numeric',numCases) + 1;
  }
  #if the conditioning set (cs) is empty, lets take it easy.
  if (is.na(csIndex) || length(csIndex) == 0 || csIndex == 0){
    
    #perform the test. If the coxph function launches a warning, the
    #function returns "flag=0", that means "the test cannot be performed"
    
    #fitting the model
    tryCatch(
      cox_results <- survival::coxph(target ~ x, weights = wei),
      warning=function(w) {
        #Do nothing...
      }
    )
    if (is.null(cox_results)){
      return(results);
    }
    
    #retrieve the p value and stat. 
    res <- anova(cox_results)
    dof <- res[2, 3]
    stat <- res[2, 2]
    pvalue <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE);
    
    if( hash )#update hash objects
    {
      stat_hash[[key]] <- stat;#.set(stat_hash , key , stat)
      pvalue_hash[[key]] <- pvalue;#.set(pvalue_hash , key , pvalue)
    }
    
    flag = 1;
    
  }else{
 
    tryCatch(
      
      #fitting the full model
      cox_results_full <- survival::coxph(target ~ ., data = as.data.frame(  dataset[ , c(csIndex, xIndex)] ), weights = wei ),
      
      warning=function(w) {
        #Do nothing
      }
    )
    if ( is.null(cox_results_full) ){
      return(results);
    }

    res = anova(cox_results_full)
    pr = nrow(res)
    stat = res[pr, 2]
    dF = res[pr, 3]
    pvalue = pchisq(stat, dF, lower.tail = FALSE, log.p = TRUE)
    
    if( hash )#update hash objects
    {
      stat_hash[[key]] <- stat;#.set(stat_hash , key , stat)
      pvalue_hash[[key]] <- pvalue;#.set(pvalue_hash , key , pvalue)
    }
    
    flag = 1;
    
  }
  
  results = list(pvalue = pvalue, stat = stat, flag = flag, stat_hash=stat_hash, pvalue_hash=pvalue_hash);
  return(results);
  
}