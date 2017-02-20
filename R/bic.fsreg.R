bic.fsreg <- function( target, dataset, test = NULL, wei = NULL, tol = 2, robust = FALSE, ncores = 1 ) {

  p <- ncol(dataset)  ## number of variables
  bico <- numeric( p )
  moda <- list()
  k <- 1   ## counter
  n <- length(target)  ## sample size
  con <- log(n)
  tool <- NULL
  info <- matrix( 0, ncol = 2 )
  result <- NULL
  sela <- NULL
  
  #check for NA values in the dataset and replace them with the variable median or the mode
  if( any(is.na(dataset)) ) {
    #dataset = as.matrix(dataset);
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    if (class(dataset) == "matrix")  {
      dataset <- apply( dataset, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
    }else{
      poia <- which( is.na(dataset), arr.ind = TRUE )[2]
      for( i in poia )  {
        xi <- dataset[, i]
        if(class(xi) == "numeric")
        {                    
          xi[ which( is.na(xi) ) ] <- median(xi, na.rm = TRUE) 
        } else if ( is.factor( xi ) ) {
          xi[ which( is.na(xi) ) ] <- levels(xi)[ which.max( as.vector( table(xi) ) )]
        }
        dataset[, i] <- xi
      }
    }
  }
  
  ##################################
  # target checking and initialize #
  ##################################
  

  if ( is.null( colnames(dataset) ) )   colnames(dataset) <- paste("X", 1:p, sep = "")

  if ( is.null(test) ) {

    ## linear regression 
    if ( sum( class(target) == "numeric" ) == 1  ||  sum( class(target) == "vector" ) == 1  ) {
     
	 la <- length( unique(target) )
	 
	 if ( la > 2 ) {
	   if ( sum( round(target) - target ) == 0 )   test <- "testIndPois" 
		
	 } else if ( la == 2 ) {
	    test <- "testIndLogistic"  
		
     } else if ( min( target ) > 0  &  max( target ) < 1 )   target <- log( target / (1 - target) ) 
      test <- "testIndReg"  
      
    }
      
    ## surival data
    if ( sum( class(target) == "Surv" ) == 1 ) {
      test <- "censIndCR"
    }

    ## ## ordinal, multinomial or perhaps binary data 
    if ( length( unique(target) ) == 2  &  is.factor(target) ) {
      test <- "testIndLogistic"   
    }
 
  }

    #available conditional independence tests
    av_models = c("testIndReg", "testIndRQ", "testIndBeta", "testIndCR", "testIndWR", "testIndLogistic", "testIndPois", "testIndNB", "testIndZIP", "testIndSpeedglm");
    
    #cat(test)
    
  if ( ( test == "testIndLogistic"  &  length( unique(target) ) == 2 )  ||  test == "testIndPois"  ||  test == "testIndReg" ) {
   
    result <- bic.glm.fsreg( target, dataset, wei = wei, tol = tol, heavy = FALSE, robust = robust, ncores = ncores ) 
  
  } else if ( test == "testIndspeedglm" ) {
    
    result <- bic.glm.fsreg( target, dataset, wei = wei, tol = tol, heavy = TRUE, robust = robust, ncores = ncores ) 
    
  } else if ( test == "testIndBeta" ) {
    
    result <- bic.betafsreg(target, dataset, wei = wei, tol = tol, ncores = ncores )
  
  } else if ( test == "testIndZip" ) {
    
    result <- bic.zipfsreg(target, dataset, wei = wei, tol = tol, ncores = ncores )
    
  } else {
 
    ci_test <- test <- match.arg(test, av_models ,TRUE);
    #convert to closure type

  if ( test == "censIndCR" ) {
      test <- survival::coxph 
      robust <- FALSE

	} else if ( test == "censIndWR" ) {
      test <- survival::survreg 
      robust <- FALSE

  } else if ( test == "testIndLogistic" ) {
	  if ( is.ordered(target) ) {      
        test <- ordinal::clm
        robust <- FALSE
      } else {
		    test <- nnet::multinom  
	      robust <- FALSE
	    }
      
  } else if ( test == "testIndNB" ) {
      test <- MASS::glm.nb
      robust <- FALSE

	} else if ( test == "testIndRQ" ) {
      test <- quantreg::rq 
      robust <- FALSE

  } else if ( test == "testIndReg"  &  !is.matrix(target) ) {
	    if ( !robust ) {
        test <- lm 
      } else  test <- MASS::rlm
  }
    
    runtime <- proc.time()
      
    ini <- test( target ~ 1 )
    ini <-  - 2 * as.numeric( logLik(ini) ) + con  ## initial BIC
    
    if (ncores <= 1) {
        for (i in 1:p) {
          mi <- test( target ~ dataset[, i], weights = wei )
		      la <- logLik(mi)
          bico[i] <-  - 2 * as.numeric( la ) + attr(la, "df") * con
        }

      mat <- cbind(1:p, bico)

      if( any( is.na(mat) ) )     mat[ which( is.na(mat) ) ] <- ini

    } else {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      bico <- numeric(p)
      mod <- foreach( i = 1:p, .combine = rbind) %dopar% {
        ww <- test( target ~ dataset[, i], weights = wei )
		    la <- logLik(ww)
        bico[i] <-  - 2 * as.numeric( la ) +  attr(la, "df") * con
      }
      stopCluster(cl)

      mat <- cbind(1:p, mod)

      if ( any( is.na(mat) ) )  mat[ which( is.na(mat) ) ] <- ini
      
    }

    colnames(mat) <- c("variable", "BIC")
    rownames(mat) <- 1:p
    sel <- which.min( mat[, 2] )
    
    if ( ini - mat[sel, 2] > tol ) {
      info[1, ] <- mat[sel, ]
      mat <- mat[-sel, ]
      if ( !is.matrix(mat) )   mat <- matrix(mat, ncol = 2) 
      sela <- sel
      mi <- test( target ~ dataset[, sel], weights = wei )
      la <- logLik(mi)
      tool[1] <-  - 2 * as.numeric( la ) +  attr(la, "df") * con
      moda[[ 1 ]] <- mi
    } else  {
      info <- info  
      sela <- NULL
    }
    
    ######
    ###     k equals 2
    ######

    if ( length(moda) > 0  &  nrow(mat) > 0 ) {

      k <- 2
      pn <- p - k  + 1
      mod <- list()

      if ( ncores <= 1 ) {
        bico <- numeric( pn )
        for ( i in 1:pn ) {
          ma <- test( target ~., data = as.data.frame( dataset[, c(sel, mat[i, 1]) ] ), weights = wei )
		      la <- logLik(ma)
          bico[i] <-  - 2 * as.numeric( la ) +  attr(la, "df") * con
        }

        mat[, 2] <- bico

      } else {

        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        bico <- numeric(pn)
        mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
          ww <- test( target ~., data = as.data.frame( dataset[, c(sel, mat[i, 1]) ] ), weights = wei )
          la <- logLik(ww)
          bico[i] <-  - 2 * as.numeric( la ) +  attr(la, "df") * con
        }
        stopCluster(cl)

        mat[, 2] <- mod

      }

      ina <- which.min( mat[, 2] )
      sel <- mat[ina, 1]

      if ( tool[1] - mat[ina, 2] <= tol ) {
        info <- info

      } else {
        tool[2] <- mat[ina, 2]
        info <- rbind(info, mat[ina, ] )
        sela <- info[, 1]
        mat <- mat[-ina , ]
        if ( !is.matrix(mat) )  mat <- matrix(mat, ncol = 2) 
        mi <- test( target ~., data = as.data.frame( dataset[, sela] ), weights = wei )
	      la <- logLik(mi)
        tool[2] <-  - 2 * as.numeric( la ) +  attr(la, "df") * con
        moda[[ 2 ]] <- mi
      }
   }

      #########
      ####      k is greater than 2
      #########

    if ( nrow(info) > 1  &  nrow(mat) > 0 ) {
      while ( ( k < n - 10 ) & ( tool[ k - 1 ] - tool[ k ] > tol ) & ( nrow(mat) > 0 ) ) {

        k <- k + 1
        pn <- p - k + 1

        if (ncores <= 1) {
          for ( i in 1:pn ) {
            ma <- test( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), weights = wei )
		        la <- logLik(ma)
            mat[i, 2] <-  - 2 * as.numeric( la ) +  attr(la, "df") * con
          }

        } else {
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          bico <- numeric(pn)
          mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
            ww <- test( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), weights = wei )
			      la <- logLik(ww)
            bico[i] <-  - 2 * as.numeric( la ) +  attr(la, "df") * con
          }
          stopCluster(cl)

          mat[, 2] <- mod

        }

        ina <- which.min( mat[, 2] )
        sel <- mat[ina, 1]

        if ( tool[k - 1] - mat[ina, 2]  <= tol ) {
          info <- rbind( info,  c( -10, Inf ) )
          tool[k] <- Inf

        } else {

          tool[k] <- mat[ina, 2]
          info <- rbind(info, mat[ina, ] )
          sela <- info[, 1]
          mat <- mat[-ina , ]
          if ( !is.matrix(mat) )  mat <- matrix(mat, ncol = 2) 
          ma <- test( target ~., data =as.data.frame( dataset[, sela] ), weights = wei )
		      la <- logLik(ma)
          tool[k] <-  - 2 * as.numeric( la ) +  attr(la, "df") * con
          moda[[ k ]] <- ma

        }

      }

    }

    runtime <- proc.time() - runtime

    d <- length(sela)
    final <- NULL
    models <- NULL

    if ( d >= 1 ) {
      models <- NULL
      xx <- as.data.frame( dataset[, sela] )
      colnames(xx) <- paste("V", sela, sep = "")

      if ( d == 1 ) {
        xx <- as.data.frame( dataset[, sela] )
        colnames(xx) <- paste("V", sela, sep = "")
        models <- final <- test( target ~., data = as.data.frame( xx ), weights = wei )

      } else    for (i in 1:d)  models[[ i ]] <- test( target ~., data = as.data.frame( xx[, 1:i] ), weights = wei )
      
      final <- summary( models[[ d ]] )
      info <- info[1:d, ]
      if ( d == 1 )  info <- matrix(info, nrow = 1)
      colnames(info) <- c( "variables", "BIC" )
      rownames(info) <- info[, 1]
    }

    result = list(mat = t(mat), info = info, models = models, final = final, ci_test = ci_test, runtime = runtime )
  } 

  result

}
