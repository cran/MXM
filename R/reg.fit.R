##########
### A generic regression model accepting many regression models
###
##########

reg.fit <- function(y, dataset, event = NULL, reps = NULL, group = NULL, slopes = FALSE, 
                    reml = FALSE, model = NULL, robust = FALSE, wei = NULL, xnew = NULL) {
  
  ## possible models are "gaussian" (default), "binary", "binomial", "multinomial", "poisson",
  ## "ordinal", "Cox", "Weibull", "exponential", "zip", "beta", "median", "negbin",
  ## "longitudinal" or "grouped".
  ## robust is either TRUE or FALSE
  ## y is the target variable, can be a numerical variable, a matrix, a factor, ordinal factor, percentages, or time to event
  ## dataset is the indendent variable(s). It can be a vector, a matrix or a dataframe with continuous only variables, 
  ## a data frame with mixed or only categorical variables
  ## event is NULL unless you have time to event data (survival regression)
  ## reps is NULL unless you have time measurements (longitudinal data)
  ## group is NULL unless you have grouped (or clustered) data or longitudinal data (needs reps)
  ## slopes is for the longitudinal data only, TRUE or FALSE 
  ## reml is for mixed models. If TRUE, REML will be used, otherwise ML will be used

  x <- as.data.frame(dataset)  ## just in case
  if ( is.null( colnames(x) ) )  colnames(x) <- paste("X", 1:ncol(x), sep = "")
  la <- length( unique(target) )

  if ( is.null(model) ) {

    ## linear regression 
    if ( sum( class(y) == "numeric" ) == 1 & is.null(event) & is.null(reps)  &  is.null(group) )  model <- "gaussian"  
    
    ## multivariate data
    if ( sum( class(y) == "matrix" ) == 1 ) { 
      if ( min(y) > 0 &  sd(Rfast::rowsums(y) == 0 ) )  y <- log(y[, -1] / y[, 1])  ## compositional data
      model <- "gaussian"
    }
    
    ## percentages
    if ( is.vector(y) )  {
      if ( all(y >0 & y < 1) )  y <- log( y / (1 - y) ) 
      model <- "gaussian"
    }
    
    ## surival data
    if ( !is.null(event) ) {
      target <- survival::Surv(time = y, event = event)
      model <- "Cox"
    }

    ## longitudinal data
    if ( !is.null(reps) & !is.null(group) )  model <- "longitudinal"
    
    ## grouped data
    if ( is.null(reps) & !is.null(group) )  model <- "grouped"
    
    ## binary data
    if ( la == 2 )  model <- "binary"   
    
    ## ordinal, multinomial or perhaps binary data
    if ( is.factor(y) ) {
      if ( !is.ordered(y) ) {
        if ( la == 2 ) {
          y <- as.vector(y)
          model <- "binary"
        } else  model <- "multinomial"

      } else {
        if ( la == 2 ) {
          y <- as.vector(y)
          model <- "binary"
        } else  model <- "ordinal"    
      }
    }
    
    ## count data
    if ( sum( is.vector(y) ) == 1 ) {
      if ( ( sum( floor(y) - y ) == 0  &  la > 2 ) )  model <- "poisson"
    }
  
  }  
  ##### model checking
 
     ## univariate gaussian model
  if ( model == "gaussian"  &  is.vector(y) ) {
     if ( !robust ) {
       mod <- lm(y ~ ., data = as.data.frame(x), weights = wei )
     } else {
       # cont = robust::lmRob.control(mxr = 2000, mxf = 2000, mxs = 2000 )
       # mod <- robust::lmRob(y ~ ., data = as.data.frame(x), control = cont) 
       mod <- MASS::rlm(y ~., data=as.data.frame(x), maxit = 2000, weights = wei)
     }

  } else if ( model == "gaussian"  &  sum( class(y) == "matrix" ) == 1 ) {
    mod <- lm(y ~ ., data = as.data.frame(x) )
 
  } else if ( model == "binomial" &  sum( class(y) == "matrix" ) == 1 ) {
    mod <- glm(y[, 1] / y[, 2] ~ ., data = as.data.frame(x), weights = y[, 2], family = binomial )
    
    ## median (quantile) regression 
  } else if ( model == "median" ) {
    mod <- quantreg::rq(y ~ ., data = as.data.frame(x), weights = wei ) 

    ## binary logistic regression
  } else if ( model == "binary" ) {
    #if ( robust == FALSE ) {
      mod <- glm(y ~ ., data = as.data.frame(x), binomial, weights = wei)
    #} else {
    #  mod <- robust::glmRob(y ~ ., data = as.data.frame(x), binomial, maxit = 100)
    #}

    ## multinomial logistic regression
  } else if ( model == "multinomial" ) {
    mod <- nnet::multinom(y ~ ., data = as.data.frame(x), trace = FALSE, weights = wei)

    ## ordinal logistic regression
  } else if ( model == "ordinal" ) {
    mod <- ordinal::clm(target ~ ., data = as.data.frame(x), weights = wei )

    ## poisson regression
  } else if ( model == "poisson" ) {
    #if ( robust == FALSE ) {
      mod <- glm(y ~ ., data = as.data.frame(x), poisson, weights = wei)
    #} else {
    #  mod <- robust::glmRob(y ~ ., data = as.data.frame(x), poisson, maxit = 100)
    #}
 
    ## negative binomial regression
  } else if ( model == "negbin" ) {
    mod <- MASS::glm.nb(y ~ ., data = as.data.frame(x), weights = wei )

    ## zero inflated poisson regression with constant zero part
  } else if ( model == "zip" ) {
    mod <- zip.mod(y, x, wei = wei)

    ## beta regression
  } else if ( model == "beta" ) {
    mod <- beta.mod(y, x, wei = wei )
  
    ## Cox proportional hazards
  } else if ( model == "Cox" ) {
    model <- survival::coxph(target ~ ., data = as.data.frame(x), weights = wei )

    ## Weibull regression
  } else if ( model == "Weibull" ) {
    model <- survival::survreg(target ~ ., data = as.data.frame(x), weights = wei )

    ## Exponential regression
  } else if ( model == "Exponential" ) {
    model <- survival::survreg(target ~ ., data = as.data.frame(x), weights = wei, dist = "exponential" )

    ## (generalised) linear mixed models for longitudinal data
  } else if ( model == "longitudinal" ) {
    if ( la > 2  &  sum( round(y) - y ) != 0  ) {
      if ( slopes ) {
        mod <- lme4::lmer( target ~ . -group + (reps|group), REML = reml, data = as.data.frame( cbind(reps, x) ), weights = wei ) 
      } else  mod <- lme4::lmer( target ~ . -group + (1|group), REML = reml, data = as.data.frame( cbind(reps, x) ), weights = wei )
	  
    } else if ( la > 2  &  sum( round(y) - y ) == 0 ) {
      if ( slopes ) {
        mod <- lme4::glmer( target ~ . - group + (reps|group), REML = reml, family = poisson,  data = ( cbind(reps, x) ), weights = wei ) 
      } else  mod <- lme4::glmer( target ~ . -group  + (1|group), REML = reml, family = poisson, data = ( cbind(reps, x) ), weights = wei )
	  
    } else  if ( la == 2 ) {
      y <- as.vector(y)   
      if ( slopes ) {
        mod <- lme4::glmer( target ~ . - group + (reps|group), REML = reml, family = binomial, data = as.data.frame( cbind(reps, x) ), weights = wei ) 
      } else  mod <- lme4::glmer( target ~ . -group + (1|group), REML = reml, family = binomial, data = as.data.frame( cbind(reps, x) ), weights = wei )
    }

    ## (generalised) linear mixed models for grouped data
  } else if ( model == "grouped" ) {
    if ( sum( round(y) - y ) != 0  &  la > 2 ) {
      mod <- lme4::lmer( target ~ . -group + (1|group), REML = reml, data = as.data.frame(x), weights = wei )
    } else if ( sum( round(y) - y ) == 0  &  la > 2) {
      mod <- lme4::glmer( target ~ . -group + (1|group), REML = reml, family = poisson, data = as.data.frame(x), weights = wei )
    } else if ( la == 2 ) {
      mod <- lme4::glmer( target ~ . -group + (1|group), REML = reml, family = binomial, data = as.data.frame(x), weights = wei )
    }
  }
  
  if ( !is.null(xnew) )  {
    xnew <- as.data.frame(xnew)  
    colnames(xnew) <- colnames(x)
    pred <- predict(mod, xnew) 	
    result <- list(mod = mod, pred = pred)
  } else result <- list(mod = mod) 

  result 

}













    
      

