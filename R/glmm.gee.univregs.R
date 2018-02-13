glmm.gee.univregs <- function(target, reps = NULL, id, dataset, targetID = -1, test, wei = NULL, 
                              slopes = FALSE, correl = "echangeable", se = "jack", ncores = 1) {
  
  if ( identical(test, testIndGLMMReg)  |  identical(test, testIndGLMMPois)  |  identical(test, testIndGLMMLogistic) ) {
    results <- univariateScore.temporal(target = target, reps = reps, group = id, dataset = dataset, test = test, wei = wei, 
                                        targetID = targetID, slopes = slopes, ncores = ncores) 
    
  } else  if ( identical(test, testIndGEEReg)  |  identical(test, testIndGEEPois)  |  identical(test, testIndGEELogistic)  |
               identical(test, testIndGEENormLog)  |  identical(test, testIndGEEGamma)  |  identical(test, testIndGEEOrdinal) ) {
    results <- univariateScore.gee(target = target, reps = reps, group = id, dataset = dataset, test = test, wei = wei, 
                                   targetID = targetID, correl = correl, se = se, ncores = ncores)
  }
  
  results
}
    