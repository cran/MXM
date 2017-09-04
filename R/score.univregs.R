score.univregs <- function(target, dataset, test) {
  
  univariateModels <- list();
  dm <- dim(dataset)
  rows <- dm[1]
  cols <- dm[2]
  id <- NULL
  ina <- NULL
  if ( identical(test, testIndLogistic) )  ina <- target
  id <- Rfast::check_data(dataset, ina)
  if ( sum(id > 0) )  dataset[, id] <- rnorm(rows * length(id) )

  ## Beta regression 
  if ( identical(test, testIndBeta) ) {
    mod <- Rfast::score.betaregs(target, dataset, logged = TRUE )
    univariateModels$stat = mod[, 1]
    univariateModels$pvalue = mod[, 2]
  ## Negative Binomial 
  } else if ( identical(test, testIndNB) ) {
    mod <- Rfast::score.negbinregs(target, dataset, logged = TRUE )
    univariateModels$stat = mod[, 1]
    univariateModels$pvalue = mod[, 2]
  ## Poisson
  } else if ( identical(test, testIndPois) ) {
    mod <- Rfast::score.glms(target, dataset, oiko = "poisson", logged = TRUE )
    univariateModels$stat = mod[, 1]
    univariateModels$pvalue = mod[, 2]
  ## logistic or multinomial regression
  } else if ( identical(test, testIndLogistic) ) { 
    mod <- Rfast::score.multinomregs(target, dataset, logged = TRUE)  
    univariateModels$stat = mod[, 1]
    univariateModels$pvalue = mod[, 2]
  } else if ( identical(test, testIndGamma) ) { 
    mod <- Rfast::score.gammaregs(target, dataset, logged = TRUE)  
    univariateModels$stat = mod[, 1]
    univariateModels$pvalue = mod[, 2]
  ## Weibull
  #} else if ( identical(test, censIndWR) ) {
  #  mod <- Rfast::score.weibregs(target, dataset, logged = TRUE )
  #  univariateModels$stat = mod[, 1]
  #  univariateModels$pvalue = mod[, 2]
  #  univariateModels$flag = numeric(cols) + 1; 
  } else univariateModels <- NULL
  
  if ( !is.null(univariateModels) )  {
    univariateModels$flag = numeric(cols) + 1  
    if ( sum(id>0) > 0 ) {
      univariateModels$stat[id] = 0
      univariateModels$pvalue[id] = 1
    }
  }
  univariateModels
  
}