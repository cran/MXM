univregs <- function(target, dataset, test, wei = NULL, robust = FALSE, ncores = 1) {
  
  la <- length( unique(target) )
  univariateModels <- list();
  
  dm <- dim(dataset)
  rows <- dm[1]
  cols <- dm[2]

if ( identical(test, testIndFisher)  &  robust == FALSE )  { ## Pearson's correlation 
  a <- as.vector( cor(target, dataset) )
  dof <- rows - 3; #degrees of freedom
  wa <- 0.5 * log( (1 + a) / (1 - a) ) * sqrt(dof)
  univariateModels$stat = abs(wa);
  univariateModels$pvalue = log(2) + pt(-abs(wa), dof, log.p = TRUE) ;
  univariateModels$flag = numeric(cols) + 1;
  ## univariateModels$stat_hash = stat_hash;
  ## univariateModels$pvalue_hash = pvalue_hash;
  
} else if ( identical(test, testIndSpearman) ) {  ## Spearman's correlation
  a <- as.vector( cor(target, dataset) )
  dof <- rows - 3; #degrees of freedom
  wa <- 0.5 * log( (1 + a) / (1 - a) ) * sqrt(dof) / 1.029563
  univariateModels$stat = abs(wa) 
  univariateModels$pvalue = log(2) + pt(-abs(wa), dof, log.p = TRUE);
  univariateModels$flag = numeric(cols) + 1;

} else if ( identical(test, gSquare) ) {  ## G^2 test
  z <- cbind(target, dataset)
  dc <- Rfast::colrange(z, cont = FALSE)
  a <- Rfast::g2Test_univariate(z, dc)
  stat <- a$statistic[ a$x == 1 ]
  univariateModels$stat = stat
  univariateModels$pvalue = pchisq(stat, a$df[ a$x == 1 ], lower.tail = FALSE, log.p = TRUE)
  univariateModels$flag = numeric(cols) + 1;

} else if ( identical(test, testIndBeta) ) {  ## Beta regression
  mod <- beta.regs(target, dataset, wei, logged = TRUE, ncores = ncores)
  univariateModels$stat = mod[, 1]
  univariateModels$pvalue = mod[, 2]
  univariateModels$flag = numeric(cols) + 1;

} else if ( identical(test, testIndReg)  &  robust  ) {  ## M (Robust) linear regression
  fit1 = MASS::rlm(target ~ 1, weights = wei)
  lik1 = as.numeric( logLik(fit1) )
  lik2 = numeric(cols)
  dof = numeric(cols)
  ina <- 1:cols

  if ( ncores <= 1 | is.null(ncores) ) {
    
    for ( i in ina ) {
      fit2 = MASS::rlm(target ~ dataset[, i], weights = wei )
      lik2[i] = as.numeric( logLik(fit2) )
      dof[i] = length( coef(fit2) ) - 1
    } 
    
    stat = 2 * abs(lik1 - lik2)
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    univariateModels$flag = numeric(cols) + 1;

  } else {
    
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = ina, .combine = rbind, .packages = "MASS") %dopar% {
      fit2 = rlm(target ~ dataset[, i], weights = wei )
      lik2 = as.numeric( logLik(fit2) )
      return( c(lik2, length( coef(fit2) ) ) )
    } 
    
    stopCluster(cl)
    
    lik1 = as.numeric( logLik(fit1) )
    stat = as.vector( 2 * abs(lik1 - mod[, 1]) )
    dof = as.vector( mod[, 2] ) - 1 
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    univariateModels$flag = numeric(cols) + 1
 
  }   
  
} else if ( identical(test, testIndReg)  &  robust == FALSE  &  is.matrix(dataset)  &  is.null(wei) ) {  ## linear regression
  mod = Rfast::univglms(target, dataset, logged = TRUE) 
  univariateModels$stat = mod[, 1]
  univariateModels$pvalue = mod[, 2]
  univariateModels$flag = numeric(cols) + 1;

} else if ( identical(test, testIndReg)  &  robust == FALSE  &  is.data.frame(dataset)  &  is.null(wei) ) {  ## linear regression
  mod <- Rfast::regression(dataset, target)
  univariateModels$stat = mod[1, ]
  univariateModels$pvalue = pf(mod[1, ], mod[2, ], rows - mod[2, ], lower.tail = F, log.p = FALSE)
  univariateModels$flag = numeric(cols) + 1;

} else if ( identical(test, testIndSpeedglm)  ) {  ## big glm regresssion
  if ( is.factor(target)  ||  la == 2 ) {
    
    target <- as.numeric( as.factor(target) ) - 1 
    ina <- 1:cols

    if ( is.matrix(dataset)  &  is.null(wei)  ) {
      mod <- Rfast::univglms(target, dataset, oiko = "binomial", logged = TRUE)
      stat <- mod[, 1]
      pval <- mod[, 2]
      
    } else {
      stat <- dof <- numeric(cols)
      for ( i in ina ) {
        fit2 <- speedglm::speedglm(target ~., data = data.frame(dataset[, i]), family = binomial(logit), weights = wei )
        stat[i] <- fit2$deviance
        dof[i] <- length( coef(fit2) )     
      }
      stat <- abs( stat - fit2$nulldev )
      pval <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE )
    }  
    
  } else if ( la > 2  &  sum(round(target) - target) == 0 ) {
    ina <- 1:cols

    if (  is.matrix(dataset)  &  is.null(wei) ) {
      mod <- Rfast::univglms(target, dataset, oiko = "poisson", logged = TRUE)
      stat <- mod[, 1]
      pval <- mod[, 2]
      
    } else {
      stat <- dof <- numeric(cols)
      for ( i in ina ) {
        fit2 <- speedglm::speedglm(target ~., data = data.frame(dataset[, i]), family = poisson(log), weights = wei )
        stat[i] <- fit2$deviance
        dof[i] <- length( coef(fit2) )
      }
      stat <- abs( stat - fit2$nulldev )
      pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE )
    } 
    
  } else {
     
    if ( is.null(wei) ) {
      if ( is.matrix(dataset) )  mod <- Rfast::univglms(target, dataset, oiko = "normal", logged = TRUE) 
      if ( is.data.frame(dataset) )  mod <- Rfast::regression(dataset, target)
      stat <- mod[1, ]
      pval <- pf(stat, mod[2, ], cols - mod[2, ] - 1, lower.tail = FALSE, log.p = TRUE)
      
    } else {
      stat <- dof <- numeric(cols)
      ina <- 1:cols

      for ( i in ina ) {
        fit2 <- speedglm::speedlm(target ~., data = data.frame(dataset[, i]), weights = wei )
        suma <- summary(fit2)[[ 13 ]]
        stat[i] <- suma[1]
        dof[i] <- suma[3]
      }
      pvalue <- pf(stat, dof, cols - dof, lower.tail = FALSE, log.p = TRUE)
      
    } 
    
  }
  
  univariateModels$stat = stat
  univariateModels$pvalue = pval
  univariateModels$flag = numeric(cols) + 1;

} else if ( identical(test, testIndLogistic)  &  is.ordered(target) ) {  ## ordinal regression
  lik2 <- numeric(cols)
  dof <- numeric(cols)
  ina <- 1:cols

  fit1 <- ordinal::clm(target ~ 1, weights = wei)
  lik1 <- as.numeric( logLik(fit1) )
  df1 <- length( coef(fit1) )
  
  if ( ncores <= 1 | is.null(ncores) ) {
    
    for ( i in ina ) {
      mat <- model.matrix(target ~ dataset[, i] )
      fit2 <- ordinal::clm.fit(target, mat, weights = wei)
      lik2[i] <- as.numeric( fit2$logLik )
      dof[i] <- length( coef(fit2) ) - df1
    }
    
    stat = as.vector( 2 * abs(lik1 - lik2) )
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    univariateModels$flag = numeric(cols) + 1;

  } else {
    
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = ina, .combine = rbind, .packages = "ordinal") %dopar% {
      mat <- model.matrix(target ~ dataset[, i] )
      fit2 <- ordinal::clm.fit(target, mat, weights = wei)
      lik2 <- as.numeric( fit2$logLik )
      return( c(lik2, length( coef(fit2) ) ) )
    } 
    stopCluster(cl)
    
    lik1 = as.numeric( logLik(fit1) )
    stat = as.vector( 2 * abs(lik1 - mod[, 1]) )
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, mod[, 2] - df1, lower.tail = FALSE, log.p = TRUE)
    univariateModels$flag = numeric(cols) + 1

  }
  
} else if ( identical(test, testIndLogistic) == TRUE  &  la > 2  ) {  ## multinomial regression
  
  target = as.factor( as.numeric( target ) );
  lik2 = numeric(cols)
  dof = numeric(cols)
  fit1 = nnet::multinom(target ~ 1, trace = FALSE, weights = wei)
  lik1 = as.numeric( logLik(fit1) )
  df1 = length( coef(fit1) )
  ina <- 1:cols

  if ( ncores <= 1 | is.null(ncores) ) {
    
    for ( i in ina ) {
      fit2 = nnet::multinom(target ~ dataset[, i], trace = FALSE, weights = wei )
      lik2[i] = as.numeric( logLik(fit2) )
      dof[i] = length( coef(fit2) ) - df1
    }
    
    stat = as.vector( 2 * abs(lik1 - lik2) )
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    univariateModels$flag = numeric(cols) + 1;

  } else {
    
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = ina, .combine = rbind, .packages = "nnet") %dopar% {
      fit2 = nnet::multinom(target ~ dataset[, i], weights = wei)
      lik2 = as.numeric( logLik(fit2 ) )
      return( c(lik2, length( coef(fit2) ) ) )
    }
    stopCluster(cl)
    
    lik1 = as.numeric( logLik(fit1) )
    stat = as.vector( 2 * abs(lik1 - mod[, 1]) )
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, mod[, 2] - df1, lower.tail = FALSE, log.p = TRUE)
    univariateModels$flag = numeric(cols) + 1
 
  }
  
} else if ( identical(test, testIndLogistic)  &  la == 2  &  is.matrix(dataset)  &  is.null(wei) ) {  ## logistic regression
  if ( is.factor(target) )   target <- as.numeric(target) - 1
  
  mod <- Rfast::univglms( target, dataset, logged = TRUE ) 
  univariateModels$stat = mod[, 1]
  univariateModels$pvalue = mod[, 2]
  univariateModels$flag = numeric(cols) + 1;

} else if ( identical(test, testIndLogistic)  &  la == 2  &  ( !is.null(wei)  ||  is.data.frame(dataset)  ) ) {  ## Logistic regression
  fit1 = glm(target ~ 1, binomial, weights = wei)
  lik1 = as.numeric( logLik(fit1) )
  lik2 = numeric(cols)
  dof = numeric(cols)
  ina <- 1:cols

  if ( ncores <= 1 | is.null(ncores) ) {

    for ( i in ina ) {
      fit2 = glm( target ~ dataset[, i], binomial, weights = wei )
      lik2[i] = as.numeric( logLik(fit2) )
      dof[i] = length( coef(fit2) ) - 1
    }
    
    stat = 2 * abs(lik1 - lik2)
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    univariateModels$flag = numeric(cols) + 1;

  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = ina, .combine = rbind) %dopar% {
      fit2 = glm( target ~ dataset[, i], binomial, weights = wei )
      lik2 = as.numeric( logLik(fit2) )
      return( c(lik2, length( coef(fit2) ) - 1 ) )
    }
    stopCluster(cl)
    
    lik1 = as.numeric( logLik(fit1) )
    stat = as.vector( 2 * abs(lik1 - mod[, 1]) )
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
    univariateModels$flag = numeric(cols) + 1

  }
  
} else if ( identical(test, testIndBinom) ) {  ## Logistic regression
  wei <- target[, 2] 
  y <- target[, 1] / wei
  fit1 = glm(y ~ 1, binomial, weights = wei)
  lik1 = as.numeric( logLik(fit1) )
  lik2 = numeric(cols)
  dof = numeric(cols)
  ina <- 1:cols

  if ( ncores <= 1 | is.null(ncores) ) {
    
    for ( i in ina ) {
      fit2 = glm( y ~ dataset[, i], binomial, weights = wei )
      lik2[i] = as.numeric( logLik(fit2) )
      dof[i] = length( coef(fit2) ) - 1
    }
    
    stat = 2 * abs(lik1 - lik2)
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    univariateModels$flag = numeric(cols) + 1;

  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    wei <- target[, 2] 
    y <- target[, 1] / wei
    mod <- foreach(i = ina, .combine = rbind) %dopar% {
      fit2 = glm( y ~ dataset[, i], binomial, weights = wei )
      lik2 = as.numeric( logLik(fit2) )
      return( c(lik2, length( coef(fit2) ) - 1 ) )
    }
    stopCluster(cl)
    
    lik1 = as.numeric( logLik(fit1) )
    stat = as.vector( 2 * abs(lik1 - mod[, 1]) )
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
    univariateModels$flag = numeric(cols) + 1

  }
  
} else if ( identical(test, testIndPois)  &  is.matrix(dataset)  &  is.null(wei) ) {  ## Poisson regression
  mod <- Rfast::univglms( target, dataset, logged = TRUE ) 
  univariateModels$stat = mod[, 1]
  univariateModels$pvalue = mod[, 2]
  univariateModels$flag = numeric(cols) + 1;

} else if ( identical(test, testIndPois)  &  ( !is.null(wei)  ||  is.data.frame(dataset)  ) ) {  ## Poisson regression
  fit1 = glm(target ~ 1, poisson, weights = wei)
  lik1 = as.numeric( logLik(fit1) )
  lik2 = numeric(cols)
  dof = numeric(cols)
  ina <- 1:cols

  if ( ncores <= 1 | is.null(ncores) ) {
    
    for ( i in ina ) {
      fit2 = glm( target ~ dataset[, i], poisson, weights = wei )
      lik2[i] = as.numeric( logLik(fit2) )
      dof[i] = length( coef(fit2) ) - 1
    }
    
    stat = 2 * abs(lik1 - lik2)
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    univariateModels$flag = numeric(cols) + 1;

  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = ina, .combine = rbind) %dopar% {
      fit2 = glm( target ~ dataset[, i], poisson, weights = wei )
      lik2 = as.numeric( logLik(fit2) )
      return( c(lik2, length( coef(fit2) ) - 1 ) )
    }
    stopCluster(cl)
    
    lik1 = as.numeric( logLik(fit1) )
    stat = as.vector( 2 * abs(lik1 - mod[, 1]) )
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
    univariateModels$flag = numeric(cols) + 1
  }
  
} else if ( identical(test, testIndNB) ) {  ## Zero-inflated Poisson regression
  lik1 <- MASS::glm.nb( target ~ 1, weights = wei )$deviance
  ina <- 1:cols

  if ( ncores <= 1 | is.null(ncores) ) {
    lik2 <- dof <- numeric(cols)
    
    for ( i in ina ) {
      fit2 = MASS::glm.nb( target ~ dataset[, i], weights = wei )
      ik2[i] = fit2$deviance
      dof[i] = length( coef(fit2) ) - 1
    }
    
    stat = abs(lik1 - lik2)
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    univariateModels$flag = numeric(cols) + 1;

  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = ina, .combine = rbind, .packages = "MASS") %dopar% {
      fit2 = MASS::glm.nb( target ~ dataset[, i], weights = wei )
      ik2 = fit2$deviance
      return( c(lik2, length( coef(fit2) ) - 1 ) )
    }
    stopCluster(cl)
    
    stat <- lik1 - as.vector(mod[, 1])
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
    univariateModels$flag = numeric(cols) + 1
  }
  
} else if ( identical(test, testIndZIP) ) {  ## Zero-inflated Poisson regression
  moda <- zip.regs(target, dataset, wei, logged = TRUE, ncores = ncores) 
  univariateModels$stat = moda[, 1]
  univariateModels$pvalue = moda[, 2]
  univariateModels$flag = numeric(cols) + 1;

} else if ( identical(test, testIndRQ) ) {  ## Median (quantile) regression
  
  fit1 = quantreg::rq(target ~ 1, weights = wei)
  stat = pval = numeric(cols)
  ina <- 1:cols

  if ( ncores <= 1 | is.null(ncores) ) {
    
    for ( i in ina ) {
      fit2 = quantreg::rq(target ~ dataset[, i], weights = wei )
      ww = anova(fit1, fit2, test = "rank")
      df1 = as.numeric( ww[[1]][1] )
      df2 = as.numeric( ww[[1]][2] )
      stat[i] = as.numeric( ww[[1]][3] )
      pval[i] = pf(stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE)
    }
    
    univariateModels$stat = stat
    univariateModels$pvalue = pval
    univariateModels$flag = numeric(cols) + 1;

  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = ina, .combine = rbind, .packages = "quantreg") %dopar% {
      fit2 = quantreg::rq(target ~ dataset[, i], weights = wei )
      ww = anova(fit1, fit2, test = "rank")
      df1 = as.numeric( ww[[1]][1] )
      df2 = as.numeric( ww[[1]][2] )
      stat = as.numeric( ww[[1]][3] )
      pval = pf(stat, df1, df2, lower.tail = FALSE, log.p = TRUE)
      return( c(stat, pval ) )
    }
    stopCluster(cl)
    
    univariateModels$stat = as.vector( mod[, 1] )
    univariateModels$pvalue = as.vector( mod[, 2] )
    univariateModels$flag = numeric(cols) + 1
  }
  
} else if ( identical(test, testIndIGreg) ) {  ## Poisson regression
  fit1 = glm(target ~ 1, family = inverse.gaussian(link = log), weights = wei)
  lik1 = as.numeric( logLik(fit1) )
  lik2 = numeric(cols)
  dof = numeric(cols)
  ina <- 1:cols

  if ( ncores <= 1 | is.null(ncores) ) {
    
    for ( i in ina ) {
      fit2 = glm( target ~ dataset[, i], family = inverse.gaussian(link = log), weights = wei )
      lik2[i] = as.numeric( logLik(fit2) )
      dof[i] = length( coef(fit2) ) - 1
    }
    
    stat = 2 * abs(lik1 - lik2)
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    univariateModels$flag = numeric(cols) + 1;

  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = ina, .combine = rbind) %dopar% {
      fit2 = glm( target ~ dataset[, i], family = inverse.gaussian(link = log), weights = wei )
      lik2 = as.numeric( logLik(fit2) )
      return( c(lik2, length( coef(fit2) ) - 1 ) )
    }
    stopCluster(cl)
    
    lik1 = as.numeric( logLik(fit1) )
    stat = as.vector( 2 * abs(lik1 - mod[, 1]) )
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
    univariateModels$flag = numeric(cols) + 1
  }
  
} else if ( identical(test, censIndCR) ) {  ## Cox regression
  stat = numeric(cols)
  dof = numeric(cols)
  ina <- 1:cols

  if ( ncores <= 1 | is.null(ncores) ) {
    
    for ( i in ina ) {
      fit2 = survival::coxph( target ~ dataset[, i], weights = wei )
      res <- anova(fit2)
      dof[i] <- res[2, 3]
      stat[i] <- res[2, 2]
    }
    
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    univariateModels$flag = numeric(cols) + 1;

  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = ina, .combine = rbind, .packages = "survival") %dopar% {
      fit2 = survival::coxph( target ~ dataset[, i], weights = wei )
      res <- anova(fit2)
      return( c(res[2, 2], res[2, 3] ) )
    }
    stopCluster(cl)

    univariateModels$stat = mod[, 1]
    univariateModels$pvalue = pchisq(mod[, 1], mod[, 2], lower.tail = FALSE, log.p = TRUE)
    univariateModels$flag = numeric(cols) + 1
  }
  
} else if ( identical(test, censIndWR) ) {  ## Weibull regression
  fit1 = survival::survreg(target ~ 1, weights = wei)
  lik1 = as.numeric( logLik(fit1) )
  lik2 = numeric(cols)
  dof = numeric(cols)
  ina <- 1:cols

  if ( ncores <= 1 | is.null(ncores) ) {
    
    for ( i in ina ) {
      fit2 = survival::survreg( target ~ dataset[, i], weights = wei )
      lik2[i] = as.numeric( logLik(fit2) )
      dof[i] = length( coef(fit2) ) - 1
    }
    
    stat = 2 * abs(lik1 - lik2)
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    univariateModels$flag = numeric(cols) + 1;

  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = ina, .combine = rbind, .packages = "survival") %dopar% {
      fit2 = survival::survreg( target ~ dataset[, i], weights = wei )
      lik2 = as.numeric( logLik(fit2) )
      return( c(lik2, length( coef(fit2) ) - 1 ) )
    }
    stopCluster(cl)
    
    stat = as.vector( 2 * abs(lik1 - mod[, 1]) )
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
    univariateModels$flag = numeric(cols) + 1
  }
  
} else if ( identical(test, testIndClogit) ) {  ## Weibull regression
  id = target[, 2] #the patient id
  case = as.logical(target[, 1]);  ## case 
  
  stat = numeric(cols)
  dof = numeric(cols)
  ina <- 1:cols

  if ( ncores <= 1 | is.null(ncores) ) {
    
    for ( i in ina ) {
      fit2 = survival::clogit( case ~ dataset[, i] + strata(id) ) 
      dof[i] = length( coef(fit2) ) 
      stat[i] = 2 * abs( diff(fit2$loglik) )
    }
    
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    univariateModels$flag = numeric(cols) + 1;

  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = ina, .combine = rbind, .packages = "survival") %dopar% {
      fit2 = survival::clogit(case ~ dataset[, i] + strata(id) ) 
      return( c(2 * abs( diff(fit2$loglik) ), length( coef(fit2) ) ) )
    }
    stopCluster(cl)
    
    stat = as.vector( 2 * abs(lik1 - mod[, 1]) )
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
    univariateModels$flag = numeric(cols) + 1
  }
  
} else if ( identical(test, censIndER) ) {  ## Weibull regression
  fit1 = survival::survreg(target ~ 1, dist = "exponential", weights = wei)
  lik1 = as.numeric( logLik(fit1) )
  lik2 = numeric(cols)
  dof = numeric(cols)
  ina <- 1:cols

  if ( ncores <= 1 | is.null(ncores) ) {
    
    for ( i in ina ) {
      fit2 = survival::survreg( target ~ dataset[, i], dist = "exponential", weights = wei )
      lik2[i] = as.numeric( logLik(fit2) )
      dof[i] = length( coef(fit2) ) - 1
    }
    
    stat = 2 * abs(lik1 - lik2)
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    univariateModels$flag = numeric(cols) + 1;

  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = ina, .combine = rbind, .packages = "survival") %dopar% {
      fit2 = survival::survreg( target ~ dataset[, i], dist = "exponential", weights = wei )
      return( c(as.numeric( logLik(fit2) ), length( coef(fit2) ) - 1 ) )
    }
    stopCluster(cl)
    
    stat = as.vector( 2 * abs(lik1 - mod[, 1]) )
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, mod[ ,2], lower.tail = FALSE, log.p = TRUE)
    univariateModels$flag = numeric(cols) + 1
  }
  
}
  
univariateModels

}