##########
### A generic regression model accepting many regression models
###
##########

reg.fit <- function(y, dataset, event = NULL, reps = NULL, group = NULL, slopes = FALSE, 
                    reml = TRUE, model = NULL, xnew = NULL) {
  ## possible models are "gaussian" (default), "binary", "multinomial", "poisson",
  ## "ordinal", "survival", "zip0", "zipx", "beta", "median", "negbin",
  ## "robust", "longitudinal" or "grouped".
  ## y is the target variable, can be a numerical variable, a matrix, a factor, ordinal factor, percentages, or time to event
  ## dataset is the indendent variable(s). It can be a vector, a matrix or a dataframe with continuous only variables, 
  ## a data frame with mixed or only categorical variables
  ## event is NULL unless you have time to event data (survival regression)
  ## reps is NULL unless you have time measurements (longitudinal data)
  ## group is NULL unless you have grouped (or clustered) data or longitudinal data (needs reps)
  ## slopes is for the longitudinal data only, TRUE or FALSE 
  ## reml is for mixed models. If TRUE, REML will be used, otherwise ML will be used

  x <- as.data.frame(dataset)  ## just in case
  if ( is.null( colnames(x) ) )  {  ## checks for column names
    colnames(x) <- paste("X", 1:ncol(x), sep = "")
  }	

  ## dependent (target) variable checking if no model was given, 
  ## but other arguments are given. For some cases, these are default cases

  if ( is.null(model) ) {

    ## linear regression 
    if ( class(y) == "numeric" & is.null(event) & is.null(reps) & is.null(group) ) {
      model <- "gaussian"  
    }
    
    ## multivariate data
    if ( class(y) == "matrix" ) {
      a <- rowSums(y)
      if ( all(y > 0) & round( sd(a, 16) ) == 0 )  ## are they compositional data?
        model <- "gaussian"
    }
    
    ## percentages
    if ( all( y > 0 & y < 1) ) {  ## are they percentages?
      y <- log( y / (1 - y) ) 
      model <- "gaussian"
    }
    
    ## surival data
    if ( !is.null(event) ) {
      target <- survival::Surv(time = y, event = event)
      model <- "survival"
    }

    ## longitudinal data
    if ( !is.null(reps) & !is.null(group) ) {
      model <- "longitudinal"
    }
    
    ## grouped data
    if ( is.null(reps) & !is.null(group) ) {
      model <- "grouped"
    }
    
    ## binary data
    if ( length( unique(y) ) == 2 ) {
      model <- "binary"   
    }
    
    ## ordinal, multinomial or perhaps binary data
    if ( is.factor(y) ) {
      if ( !is.ordered(y) ) {
        if ( length(unique(y) ) == 2 ) {
          y <- as.vector(y)
          model <- "binary"
        } else {
          model <- "multinomial"
        }  
      } else {
        if ( length(unique(y) ) == 2 ) {
          y <- as.vector(y)
          model <- "binary"
        } else {
          model <- "ordinal"    
        }
      }
    }
    
    ## count data
    if ( is.vector(y) & sum( round(y) - y ) == 0 ) {
      model <- "poisson"
    }
  
  }

  ##### model checking
 
  ## univariate or multivariate gaussian model
  if (model == "gaussian" ) {
     mod <- lm(y ~ ., data = as.data.frame(x) )

  ## median (quantile) regression 
  } else if (model == "median") {
    mod <- quantreg::rq(y ~ ., data = as.data.frame(x) ) 

  ## linear regression via M estimation  
  } else if (model == "robust") {
    mod <- MASS::rlm(y ~ ., data = as.data.frame(x), maxit = 2000) 
  
  ## binary logistic regression
  } else if (model == "binary") {
    mod <- glm(y ~ ., data = as.data.frame(x), binomial)

  ## multinomial logistic regression
  } else if (model == "multinomial"){
    mod <- nnet::multinom(y ~ ., data = as.data.frame(x), trace = FALSE)

  ## ordinal logistic regression
  } else if (model == "ordinal") {
    mod <- ordinal::clm(target ~ ., data = as.data.frame(x) )

  ## poisson regression
  } else if (model == "poisson") {
    mod <- glm(y ~ ., data = as.data.frame(x), poisson)
 
  ## negative binomial regression
  } else if (model == "negbin") {
    mod <- MASS::glm.nb(y ~ ., data = as.data.frame(x) )

  ## zero inflated poisson regression with constant zero part
  } else if (model == "zip0") {
    mod <- pscl::zeroinfl( y ~ .| 1, data = as.data.frame(x) )

  ## zero inflated poisson regression with variable zero part
  } else if (model == "zipx") {
    mod <- pscl::zeroinfl( y ~ . | ., data = as.data.frame(x) )

  ## beta regression
  } else if (model == "beta") {
    mod <- betareg::betareg(y ~ ., data = as.data.frame(x) )
  
  ## survival regression
  } else if (model == "survival") {
    model <- survival::coxph (target ~ ., data = as.data.frame(x))

  ## (generalised) linear mixed models for longitudinal data
  } else if (model == "longitudinal") {
    if ( length( unique(y) ) != 2 ) {
      if (slopes == TRUE ) {
        mod <- lme4::lmer( target ~ . -group + (reps|group), REML = reml, data = as.data.frame( cbind(reps, x) ) ) 
      } else {
        mod <- lme4::lmer( target ~ . -group + (1|group), REML = reml, data = as.data.frame( cbind(reps, x) ) )
      }
    }
    if ( sum( round(y) - y ) != 0 ) {
      if (slopes == TRUE ) {
        mod <- lme4::glmer( target ~ . - group + (reps|group), REML = reml, family = poisson,  data = ( cbind(reps, x) ) ) 
      } else {
        mod <- lme4::glmer( target ~ . -group  + (1|group), REML = reml, family = poisson, data = ( cbind(reps, x) ) )
      }
    }
    if ( length( unique(y) ) == 2 ) {
      y <- as.vector(y)   
      if (slopes == TRUE ) {
        mod <- lme4::glmer( target ~ . - group + (reps|group), REML = reml, family = binomial, data = as.data.frame( cbind(reps, x) ) ) 
      } else {
        mod <- lme4::glmer( target ~ . -group + (1|group), REML = reml, family = binomial, data = as.data.frame( cbind(reps, x) ) )
      }
    }

    ## (generalised) linear mixed models for grouped data
  } else if (model == "longitudinal") {
    if ( sum( round(y) - y ) != 0 & length( unique(y) ) != 2 ) {
      mod <- lme4::lmer( target ~ . -group + (1|group), REML = reml, data = as.data.frame(x) )
    }
    if ( sum( round(y) - y ) != 0 ) {
      mod <- lme4::glmer( target ~ . -group + (1|group), REML = reml, family = poisson, data = as.data.frame(x) )
    }
    if ( length( unique(y) ) == 2 ) {
      mod <- lme4::glmer( target ~ . -group + (1|group), REML = reml, family = binomial, data = as.data.frame(x) )
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













    
      

