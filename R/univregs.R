univregs <- function(target, dataset, targetID = -1, test = NULL, user_test = NULL, wei = NULL, dataInfo = NULL, robust = FALSE, ncores = 1) {
  
  univariateModels <- list();
  dm <- dim(dataset)
  rows <- dm[1]
  cols <- dm[2]
  if (targetID != -1 ) {
     target <- dataset[, targetID]
     dataset[, targetID] <- rbinom(rows, 1, 0.5)
  }   
  id <- NULL
  if ( !identical(test, testIndFisher) & !identical(test, testIndSpearman) ) {
    ina <- NULL
    id <- Rfast::check_data(dataset)
    if ( sum(id > 0) )  dataset[, id] <- rnorm(rows * length(id) )
  }  
  la <- length( unique(target) )
  
if ( !is.null(user_test) ) {
  univariateModels <- univariateScore(target, dataset, test = user_test, wei, dataInfo = dataInfo, targetID, robust)

} else if ( identical(test, testIndFisher)  &  robust == FALSE )  { ## Pearson's correlation 
  a <- as.vector( cor(target, dataset) )
  dof <- rows - 3; #degrees of freedom
  wa <- 0.5 * log( (1 + a) / (1 - a) ) * sqrt(dof)
  id <- which( is.na(a) )
  if ( length(id) > 0)  wa[id] <- 0
  univariateModels$stat = wa;
  univariateModels$pvalue = log(2) + pt( abs(wa), dof, lower.tail = FALSE, log.p = TRUE) ;

} else if ( identical(test, testIndSpearman) ) {  ## Spearman's correlation
  a <- as.vector( cor(target, dataset) )
  dof <- rows - 3; #degrees of freedom
  wa <- 0.5 * log( (1 + a) / (1 - a) ) * sqrt(dof) / 1.029563
  id <- which( is.na(a) )
  if ( length(id) > 0)  wa[id] <- 0
  univariateModels$stat = wa 
  univariateModels$pvalue = log(2) + pt( abs(wa), dof, lower.tail = FALSE, log.p = TRUE);

} else if ( identical(test, gSquare) ) {  ## G^2 test
  z <- cbind(dataset, target)
  if ( !is.matrix(z) )   z <- as.matrix(z)
  dc <- Rfast::colrange(z, cont = FALSE)
  a <- Rfast::g2tests(data = z, x = 1:cols, y = cols + 1, dc = dc)
  stat <- a$statistic
  univariateModels$stat = stat
  univariateModels$pvalue = pchisq(stat, a$df, lower.tail = FALSE, log.p = TRUE)

} else if ( identical(test, testIndBeta) ) {  ## Beta regression
  mod <- beta.regs(target, dataset, wei, logged = TRUE, ncores = ncores)
  univariateModels$stat = mod[, 1]
  univariateModels$pvalue = mod[, 2]

} else if ( identical(test, testIndReg)  &  robust  ) {  ## M (Robust) linear regression
  fit1 = MASS::rlm(target ~ 1, maxit = 2000, method = "MM")
  lik1 = as.numeric( logLik(fit1) )
  lik2 = numeric(cols)
  dof = numeric(cols)
  ina <- 1:cols

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 = MASS::rlm(target ~ dataset[, i], maxit = 2000, method = "MM" )
      lik2[i] = as.numeric( logLik(fit2) )
      dof[i] = length( coef(fit2) ) - 1
    } 
    stat = 2 * (lik2 - lik1)
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)

  } else {
    
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "MASS") %dopar% {
       fit2 = MASS::rlm(target ~ dataset[, i], maxit = 2000, method = "MM" )
       lik2 = as.numeric( logLik(fit2) )
       return( c(lik2, length( coef(fit2) ) ) )
    }  
    stopCluster(cl)
    stat = as.vector( 2 * abs(lik1 - mod[, 1]) )
    dof = as.vector( mod[, 2] ) - 1 
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
  }   
  
} else if ( identical(test, testIndReg)  &  !robust  &  !is.null(wei) ) {  ## Weighted linear regression
  
  univariateModels = list();
  stat = pval = numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 = lm( target ~ dataset[, i], weights = wei, y = FALSE, model = FALSE )
      tab = anova(fit2)
      stat[i] = tab[1, 4] 
      df1 = tab[1, 1]    ;   df2 = tab[2, 1]
      pval[i] = pf( stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE )
    }
    univariateModels$stat = stat
    univariateModels$pvalue = pval

  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind) %dopar% {
        ww <- lm( target ~ dataset[, i], weights = wei, y = FALSE, model = FALSE )
        tab <- anova( ww )
        stat <- tab[1, 4] 
        df1 <- tab[1, 1]   ;  df2 = tab[2, 1]
        pval <- pf( stat, df1, df2, lower.tail = FALSE, log.p = TRUE )
        return( c(stat, pval) )
    }
    stopCluster(cl)
    univariateModels$stat = mod[, 1]
    univariateModels$pvalue = mod[, 2]
  }   
  
} else if ( identical(test, testIndReg)  &  robust == FALSE  &  is.matrix(dataset)  &  is.null(wei) ) {  ## linear regression
  mod = Rfast::univglms(target, dataset, logged = TRUE) 
  univariateModels$stat = mod[, 1]
  univariateModels$pvalue = mod[, 2]

} else if ( identical(test, testIndReg)  &  robust == FALSE  &  is.data.frame(dataset)  &  is.null(wei) ) {  ## linear regression
  #mod <- Rfast::regression(dataset, target)
  #univariateModels$stat = mod[1, ]
  #univariateModels$pvalue = pf(mod[1, ], mod[2, ], rows - mod[2, ], lower.tail = F, log.p = TRUE)
  stat <- numeric(cols)
  pval <- numeric(cols)  
  for (i in 1:cols) {
    mod <- anova( lm(y ~ dataset[, i]) )
	  stat[i] <- mod[1, 4]
    pval[i] <- pf(stat[i], mod[1, 1], mod[2, 1], lower.tail = FALSE, log.p = TRUE)  		   
  }
  univariateModels$stat = stat
  univariateModels$pvalue = pval
  
} else if ( identical(test, testIndMVreg)  &  !robust  &  !is.null(wei) ) {  ## Weighted linear regression
  
  univariateModels = list();
  stat = pval = numeric(cols)
  
  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 = lm( target ~ dataset[, i], weights = wei, y = FALSE, model = FALSE )
      tab = anova(fit2)
      stat[i] = tab[2, 3] 
      df1 = tab[2, 4]    ;   df2 = tab[2, 5]
      pval[i] = pf( stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE )
    }
    univariateModels$stat = stat
    univariateModels$pvalue = pval
    
  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind) %dopar% {
      ww <- lm( target ~ dataset[, i], weights = wei, y = FALSE, model = FALSE )
      tab <- anova( ww )
      stat <- tab[2, 3] 
      df1 <- tab[2, 4]   ;  df2 = tab[2, 5]
      pval <- pf( stat, df1, df2, lower.tail = FALSE, log.p = TRUE )
      return( c(stat, pval) )
    }
    stopCluster(cl)
    univariateModels$stat = mod[, 1]
    univariateModels$pvalue = mod[, 2]
  }   
  
} else if ( identical(test, testIndSpeedglm) ) {  ## big glm regresssion
  if ( is.factor(target)  ||  la == 2 ) {
    
    target <- as.numeric( as.factor(target) ) - 1  
    if ( is.matrix(dataset)  &  is.null(wei)  ) {
      mod <- Rfast::univglms(target, dataset, oiko = "binomial", logged = TRUE)
      stat <- mod[, 1]
      pval <- mod[, 2]
      
    } else {
      stat <- dof <- numeric(cols)
      for ( i in 1:cols ) {
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
      pval <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE )
    } 
    
  } else {
     
    if ( is.null(wei) ) {
      if ( is.matrix(dataset) )  {
        mod <- Rfast::univglms(target, dataset, oiko = "normal", logged = TRUE) 
      }  
      #if ( is.data.frame(dataset) )   mod <- Rfast::regression(dataset, target)
      stat <- mod[1, ]
      pval <- pf(stat, mod[2, ], rows - mod[2, ] - 1, lower.tail = FALSE, log.p = TRUE)
	  if ( is.data.frame(dataset) ) {
      stat <- numeric(cols)
		  pval <- numeric(cols)  
	    for (i in 1:cols) {
        mod <- anova( lm(y ~ dataset[, i]) )
		    stat[i] <- mod[1, 5]
        pval[i] <- pf(stat[i], mod[1, 1], mod[2, 1], lower.tail = FALSE, log.p = TRUE)  		   
      } 
      mat <- cbind(1:cols, pval, stat)	 
	  } 	 
    } else {
      stat <- dof <- numeric(cols)
      for ( i in 1:cols ) {
        fit2 <- speedglm::speedlm(target ~., data = data.frame(dataset[, i]), weights = wei )
        suma <- summary(fit2)[[ 13 ]]
        stat[i] <- suma[1]
        dof[i] <- suma[2]
      }
      pval <- pf(stat, dof, rows - dof - 1, lower.tail = FALSE, log.p = TRUE)
    } 
    
  }
  univariateModels$stat = stat
  univariateModels$pvalue = pval

} else if ( identical(test, testIndLogistic)  &  is.ordered(target) ) {  ## ordinal regression
  lik2 <- numeric(cols)
  dof <- numeric(cols)
  fit1 <- ordinal::clm(target ~ 1, weights = wei)
  lik1 <- as.numeric( logLik(fit1) )
  df1 <- length( coef(fit1) )
  
  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      mat <- model.matrix(target ~ dataset[, i] )
      fit2 <- ordinal::clm.fit(target, mat, weights = wei)
      lik2[i] <- as.numeric( fit2$logLik )
      dof[i] <- length( coef(fit2) ) - df1
    }
    stat = 2 * (lik2 - lik1)
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)

  } else {
    
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "ordinal") %dopar% {
        mat <- model.matrix(target ~ dataset[, i] )
        fit2 <- ordinal::clm.fit(target, mat, weights = wei)
        lik2 <- as.numeric( fit2$logLik )
        return( c(lik2, length( coef(fit2) ) ) )
    } 
    stopCluster(cl)
    stat = 2 * (mod[, 1] - lik1)
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, mod[, 2] - df1, lower.tail = FALSE, log.p = TRUE)

  }
  
} else if ( identical(test, testIndLogistic) == TRUE  &  la > 2  ) {  ## multinomial regression
  
  target = as.factor( as.numeric( target ) );
  lik2 = numeric(cols)
  dof = numeric(cols)
  fit1 = nnet::multinom(target ~ 1, trace = FALSE, weights = wei)
  lik1 = as.numeric( logLik(fit1) )
  df1 = length( coef(fit1) )

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 = nnet::multinom(target ~ dataset[, i], trace = FALSE, weights = wei )
      lik2[i] = as.numeric( logLik(fit2) )
      dof[i] = length( coef(fit2) ) - df1
    }
    stat = 2 * (lik2 - lik1)
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)

  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "nnet") %dopar% {
        fit2 = nnet::multinom(target ~ dataset[, i], weights = wei)
        lik2 = as.numeric( logLik(fit2 ) )
        return( c(lik2, length( coef(fit2) ) ) )

    }
    stopCluster(cl)
    stat = 2 * (mod[, 1] - lik1)
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, mod[, 2] - df1, lower.tail = FALSE, log.p = TRUE)
  }
  
} else if ( identical(test, testIndLogistic)  &  la == 2  &  is.matrix(dataset)  &  is.null(wei) ) {  ## logistic regression
  if ( is.factor(target) )   target <- as.numeric(target) - 1
  mod <- Rfast::univglms( target, dataset, oiko = "binomial", logged = TRUE )
  univariateModels$stat = mod[, 1]
  univariateModels$pvalue = mod[, 2]

} else if ( identical(test, testIndLogistic)  &  la == 2  &  ( !is.null(wei)  ||  is.data.frame(dataset)  ) ) {  ## Logistic regression
  fit1 = glm(target ~ 1, binomial, weights = wei)
  lik1 = fit1$deviance
  lik2 = numeric(cols)
  dof = numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 = glm( target ~ dataset[, i], binomial, weights = wei )
      lik2[i] = fit2$deviance
      dof[i] = length( coef(fit2) ) - 1
    }
    stat = lik1 - lik2
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)

  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind) %dopar% {
        fit2 = glm( target ~ dataset[, i], binomial, weights = wei )
        lik2 = fit2$deviance
        return( c(lik2, length( coef(fit2) ) - 1 ) )
    }
    stopCluster(cl)
    stat = lik1 - mod[, 1]
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
  }
  
} else if ( identical(test, testIndBinom) ) {  ## Logistic regression
  wei <- target[, 2] 
  y <- target[, 1] / wei
  fit1 = glm(y ~ 1, binomial, weights = wei)
  lik1 = fit1$deviance
  lik2 = numeric(cols)
  dof = numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 = glm( y ~ dataset[, i], binomial, weights = wei )
      lik2[i] = fit2$deviance
      dof[i] = length( coef(fit2) ) - 1
    }
    stat = lik1 - lik2
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)

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
    stat = as.vector( lik1 - mod[, 1] )
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
  }
  
} else if ( identical(test, testIndPois)  &  is.matrix(dataset)  &  is.null(wei) ) {  ## Poisson regression
  mod <- Rfast::univglms( target, dataset, logged = TRUE ) 
  univariateModels$stat = mod[, 1]
  univariateModels$pvalue = mod[, 2]

} else if ( identical(test, testIndPois)  &  ( !is.null(wei)  ||  is.data.frame(dataset)  ) ) {  ## Poisson regression
  fit1 = glm(target ~ 1, poisson, weights = wei)
  lik1 = fit1$deviance
  lik2 = numeric(cols)
  dof = numeric(cols)
  ina <- 1:cols

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in ina ) {
      fit2 = glm( target ~ dataset[, i], poisson, weights = wei )
      lik2[i] = fit2$deviance
      dof[i] = length( coef(fit2) ) - 1
    }
    stat = lik1 - lik2
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)

  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = ina, .combine = rbind) %dopar% {
        fit2 = glm( target ~ dataset[, i], poisson, weights = wei )
        return( c(fit2$deviance, length( coef(fit2) ) - 1 ) )
    }
    stopCluster(cl)
    stat = as.vector( lik1 - mod[, 1] )
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
  }
  
} else if ( identical(test, testIndNB) ) {  ## Negative binomial regression
  lik1 <- MASS::glm.nb( target ~ 1, weights = wei )$twologlik

  if ( ncores <= 1 | is.null(ncores) ) {
    lik2 <- dof <- numeric(cols)
    for ( i in 1:cols ) {
      fit2 = MASS::glm.nb( target ~ dataset[, i], weights = wei )
      lik2[i] = fit2$twologlik
      dof[i] = length( coef(fit2) ) - 1
    }
    stat = lik2 - lik1
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)

  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = ina, .combine = rbind, .packages = "MASS") %dopar% {
        fit2 = MASS::glm.nb( target ~ dataset[, i], weights = wei )
        return( c(fit2$twologlik, length( coef(fit2) ) - 1 ) )
    }
    stopCluster(cl)
    stat <- as.vector(mod[, 1]) - lik1
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
  }

} else if ( identical(test, testIndNormLog)   ) {  ## Normal log link regression
  fit1 = glm(target ~ 1, family = gaussian(link = log), weights = wei)
  lik1 = fit1$deviance
  lik2 = numeric(cols)
  dof = numeric(cols)
  ina <- 1:cols
  
  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in ina ) {
      fit2 = glm( target ~ dataset[, i], family = gaussian(link = log), weights = wei )
      lik2[i] = fit2$deviance
      dof[i] = length( coef(fit2) ) - 1
    }
    stat = lik1 - lik2
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    
  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = ina, .combine = rbind) %dopar% {
      fit2 = glm( target ~ dataset[, i], family = gaussian(link = log), weights = wei )
      return( c(fit2$deviance, length( coef(fit2) ) - 1 ) )
    }
    stopCluster(cl)
    stat = as.vector( lik1 - mod[, 1] )
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
  }
   
} else if ( identical(test, testIndGamma)   ) {  ## Gamma regression
  fit1 = glm(target ~ 1, family = Gamma(link = log), weights = wei)
  lik1 = fit1$deviance
  lik2 = numeric(cols)
  dof = numeric(cols)
  ina <- 1:cols
  
  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in ina ) {
      fit2 = glm( target ~ dataset[, i], family = Gamma(link = log), weights = wei )
      lik2[i] = fit2$deviance
      dof[i] = length( coef(fit2) ) - 1
    }
    stat = lik1 - lik2
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    
  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = ina, .combine = rbind) %dopar% {
      fit2 = glm( target ~ dataset[, i], family = Gamma(link = log), weights = wei )
      return( c(fit2$deviance, length( coef(fit2) ) - 1 ) )
    }
    stopCluster(cl)
    stat = as.vector( lik1 - mod[, 1] )
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
  } 
  
} else if ( identical(test, testIndZIP) ) {  ## Zero-inflated Poisson regression
  moda <- zip.regs(target, dataset, wei, logged = TRUE, ncores = ncores) 
  univariateModels$stat = moda[, 1]
  univariateModels$pvalue = moda[, 2]

} else if ( identical(test, testIndRQ) ) {  ## Median (quantile) regression
  
  fit1 = quantreg::rq(target ~ 1, weights = wei)
  stat = pval = numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 = quantreg::rq(target ~ dataset[, i], weights = wei )
      ww = anova(fit1, fit2, test = "rank")
      df1 = as.numeric( ww[[1]][1] )
      df2 = as.numeric( ww[[1]][2] )
      stat[i] = as.numeric( ww[[1]][3] )
      pval[i] = pf(stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE)
    }
    univariateModels$stat = stat
    univariateModels$pvalue = pval

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
  }
  
} else if ( identical(test, testIndIGreg) ) {  ## Inverse Gaussian regression
  fit1 = glm(target ~ 1, family = inverse.gaussian(link = log), weights = wei)
  lik1 = as.numeric( logLik(fit1) )
  lik2 = numeric(cols)
  dof = numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
      for ( i in 1:cols ) {
        fit2 = glm( target ~ dataset[, i], family = inverse.gaussian(link = log), weights = wei )
        lik2[i] = as.numeric( logLik(fit2) )
        dof[i] = length( coef(fit2) ) - 1
      }
    stat = 2 * (lik2 - lik1)
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)

  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind) %dopar% {
        fit2 = glm( target ~ dataset[, i], family = inverse.gaussian(link = log), weights = wei )
        lik2 = as.numeric( logLik(fit2) )
        return( c(lik2, length( coef(fit2) ) - 1 ) )
    }
    stopCluster(cl)
    stat = as.vector( 2 * (mod[, 1] - lik1) )
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
  }
  
} else if ( identical(test, censIndCR) ) {  ## Cox regression
  stat = numeric(cols)
  dof = numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 = survival::coxph( target ~ dataset[, i], weights = wei)
      res <- anova(fit2)
      dof[i] <- res[2, 3]
      stat[i] <- res[2, 2]
    }
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)

  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
        fit2 = survival::coxph( target ~ dataset[, i], weights = wei )
        res <- anova(fit2)
        return( c(res[2, 2], res[2, 3] ) )
    }
    stopCluster(cl)
    univariateModels$stat = mod[, 1]
    univariateModels$pvalue = pchisq(mod[, 1], mod[, 2], lower.tail = FALSE, log.p = TRUE)
  }
  
} else if ( identical(test, censIndWR) ) {  ## Weibull regression
  fit1 = survival::survreg(target ~ 1, weights = wei)
  lik1 = as.numeric( logLik(fit1) )
  lik2 = numeric(cols)
  dof = numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 = survival::survreg( target ~ dataset[, i], weights = wei )
      lik2[i] = as.numeric( logLik(fit2) )
      dof[i] = length( coef(fit2) ) - 1
    }
    stat = 2 * (lik2 - lik1)
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)

  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
        fit2 = survival::survreg( target ~ dataset[, i], weights = wei )
        lik2 = as.numeric( logLik(fit2) )
        return( c(lik2, length( coef(fit2) ) - 1 ) )
    }
    stopCluster(cl)
    stat = as.vector( 2 * (mod[, 1] - lik1) )
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
  }
  
} else if ( identical(test, testIndTobit) ) {  ## Tobit regression
  fit1 = survival::survreg(target ~ 1, weights = wei, dist = "gaussian")
  lik1 = as.numeric( logLik(fit1) )
  lik2 = numeric(cols)
  dof = numeric(cols)
  
  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 = survival::survreg( target ~ dataset[, i], weights = wei, dist = "gaussian" )
      lik2[i] = as.numeric( logLik(fit2) )
      dof[i] = length( coef(fit2) ) - 1
    }
    stat = 2 * (lik2 - lik1)
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    
  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
      fit2 = survival::survreg( target ~ dataset[, i], weights = wei, dist = "gaussian" )
      lik2 = as.numeric( logLik(fit2) )
      return( c(lik2, length( coef(fit2) ) - 1 ) )
    }
    stopCluster(cl)
    stat = as.vector( 2 * (mod[, 1] - lik1) )
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
  }
  
} else if ( identical(test, testIndClogit) ) {  ## Conditional logistic regression
  subject = target[, 2] #the patient id
  case = as.logical(target[, 1]);  ## case 
  stat = numeric(cols)
  dof = numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 = survival::clogit( case ~ dataset[, i] + strata(subject) ) 
      dof[i] = length( coef(fit2) ) 
      stat[i] = 2 * abs( diff(fit2$loglik) )
    }
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)

  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
        fit2 = survival::clogit(case ~ dataset[, i] + strata(subject) ) 
        return( c(2 * abs( diff(fit2$loglik) ), length( coef(fit2) ) ) )
    }
    stopCluster(cl)
    stat = as.vector( 2 * (mod[, 1] - lik1) )
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, mod[, 2], lower.tail = FALSE, log.p = TRUE)
  }
  
} else if ( identical(test, censIndER) ) {  ## Exponential regression
  fit1 = survival::survreg(target ~ 1, dist = "exponential", weights = wei)
  lik1 = as.numeric( logLik(fit1) )
  lik2 = numeric(cols)
  dof = numeric(cols)

  if ( ncores <= 1 | is.null(ncores) ) {
    for ( i in 1:cols ) {
      fit2 = survival::survreg( target ~ dataset[, i], dist = "exponential", weights = wei )
      lik2[i] = as.numeric( logLik(fit2) )
      dof[i] = length( coef(fit2) ) - 1
    }
    stat = 2 * (lik2 - lik1)
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)

  } else {
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    mod <- foreach(i = 1:cols, .combine = rbind, .packages = "survival") %dopar% {
        fit2 = survival::survreg( target ~ dataset[, i], dist = "exponential", weights = wei )
        return( c(as.numeric( logLik(fit2) ), length( coef(fit2) ) - 1 ) )
    }
    stopCluster(cl)
    stat = as.vector( 2 * (mod[, 1] - lik1) )
    univariateModels$stat = stat
    univariateModels$pvalue = pchisq(stat, mod[ ,2], lower.tail = FALSE, log.p = TRUE)
  }
  
}  else   univariateModels <- NULL
  
  if ( !is.null(univariateModels) )  {
    univariateModels$flag = numeric(cols) + 1  
    if (targetID != - 1) {
      univariateModels$stat[targetID] = 0
      univariateModels$pvalue[targetID] = log(1)
    }
    if ( sum(id>0) > 0 ) {
      univariateModels$stat[id] = 0
      univariateModels$pvalue[id] = log(1)
    }
  }
  
  univariateModels
}