bs.reg <- function(target, dataset, threshold = 0.05, wei = NULL, test = NULL, user_test = NULL, robust = FALSE) {
  
  threshold <- log(threshold)
  dm <- dim(dataset)
  n <- dm[1]  ## sample size 
  p <- dm[2]  ## number of variables

  if ( p > n ) {
    res <- paste("The number of variables is hiher than the sample size. No backward procedure was attempted")
  
  } else {
  
   tic <- proc.time()

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

 if ( is.null( colnames(dataset) ) )  colnames(dataset) <- paste("X", 1:p, sep = "")
  
  la <- length( unique(target) )
 
  ## dependent (target) variable checking if no test was given, 
  ## but other arguments are given. For some cases, these are default cases
  
  if ( is.null(test)  &  is.null(user_test) ) {
    
    ## surival data
    if ( sum( class(target) == "Surv" ) == 1 ) {
      ci_test <- test <- "censIndCR"
      ## ordinal, multinomial or perhaps binary data
    } else if ( is.factor(target) ||  is.ordered(target) || length( unique(target) ) == 2 ) {
      ci_test <- test <- "testIndLogistic"
      ## count data
    } else if ( length( unique(target) ) > 2  &  !is.factor(target) ) {
      if ( sum( round(target) - target ) == 0 ) {
        ci_test <- test <- "testIndPois"
      } else  ci_test <- test <- "testIndReg"  
    }
  }
   
  av_models = c("testIndReg", "testIndBeta", "censIndCR", "testIndRQ", "censIndWR", "testIndLogistic", "testIndPois", "testIndNB", "testIndZIP", "testIndSpeedglm", "testIndBinom");
  
  ci_test <- test
  test <- match.arg(test, av_models, TRUE);

   ############ 
   ###  GLMs 
   ############
    
    heavy = FALSE
    if (test == "testIndSpeedglm")   heavy <- TRUE
    
    if ( test == "testIndPois"  ||  test == "testIndReg"  ||  ( test == "testIndLogistic"  &  la == 2 )  ||  test == "testIndSpeedglm" || test == "testIndBinom" ) {

      res <- glm.bsreg(target = target, dataset = dataset, wei = wei, threshold = exp( threshold ), heavy = heavy, robust = robust ) 
	   
    } else if ( test == "testIndBeta" ) {

      res <- beta.bsreg(target = target, dataset = dataset, wei = wei, threshold = exp( threshold ) ) 
	  
    } else if ( test == "testIndZIP" ) {

      res <- zip.bsreg(target = target, dataset = dataset, wei = wei, threshold = exp( threshold ) ) 	
 
    } else {
	   
     if ( test == "censIndCR" ) {
        test = survival::coxph 
        robust = FALSE
		
	  } else if ( test == "censIndWR" ) {
        test = survival::survreg 
        robust = FALSE
            
      } else if ( test == "testIndLogistic" ) {
	  
	    if ( is.ordered(target) )  {
          test = ordinal::clm
          robust = FALSE
       
	    } else {
          test = nnet::multinom
          robust = FALSE
        }
	  
      } else if ( test == "testIndNB" ) {
        test = MASS::glm.nb
        robust = FALSE

      } else if ( test == "testIndRQ" ) {
        test = quantreg::rq
        robust = FALSE
	
	  } else if (test == "testIndRQ") {
	    test = quantreg::rq
	    robust = FALSE
	  }

     ###################
     ###################
  
    ini <- test( target ~.,  data = as.data.frame(dataset), weights = wei )
	  dofini <- length( coef(ini) )
	  likini <- logLik(ini) 
	  stat <- dof <- numeric(p)
		
	  for (j in 1:p) {
		  mod <- test( target ~.,  data = as.data.frame(dataset[, -j]), weights = wei )
		  stat[j] <- 2 * abs( likini - logLik(mod) )
		  dof[j] <- dofini - length( coef(mod) ) 
	  }
	    
	    if ( ci_test == "censIndCR")   dof <- dof + 1
      mat <- cbind(1:p, pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
      colnames(mat) <- c("variable", "log.p-values", "statistic" )
      rownames(mat) <- 1:p 

      sel <- which.max( mat[, 2] )
      info <- matrix( c(0, -10, -10) , ncol = 3 )
      sela <- sel 

      if ( mat[sel, 2] < threshold ) {
        final <- ini 
    
      } else {
       
        info[1, ] <- mat[sel, ]
        mat <- mat[-sel, , drop = FALSE] 
        dat <- as.data.frame( dataset[, -sel] ) 

      } 
    
      i <- 1  

      if ( info[1, 2] > threshold ) {
         
         while ( info[i, 2] > threshold  &  NCOL(dat) > 0 )  {   

           ini <- test( target ~., data = dat, weights = wei )
           likini <- logLik(ini) 
           dofini <- length( coef(ini) )

           i <- i + 1        
           k <- p - i + 1
		      
		      if ( k == 1 ) {
		        mod <- test(target ~ 1, weights = wei)
		        if ( ci_test == "censIndCR")  {
		          dof <- dof + 1
		          stat <- 2 * abs( likini - mod$loglik )
		        } else {
		          stat <- 2 * abs( likini - logLik(mod) )
		          dof <- dofini - length( coef(mod) ) 
		        }  
		        pval <- pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE)
		      
		        if (pval > threshold ) {
		          final <- "No variables were selected"
		          info <- rbind(info, c(mat[, 1], stat, pval) )
		          dat <- as.data.frame( dataset[, -info[, 1] ] )
		          mat <- NULL
		        } else {
		          info <- rbind(info, c(0, -10, -10)) 
		          final <- ini
		        }  
		              
		      } else { 
            stat <- dof <- numeric(k)
             
	          for (j in 1:k) {
		          mod <- test( target ~.,  data = as.data.frame(dat[, -j]), weights = wei )
		          stat[j] <- 2 * abs( likini - logLik(mod) )
		          dof[j] <- dofini - length( coef(mod) ) 
	          }
            
            mat[, 2:3] <- cbind( pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
            sel <- which.max( mat[, 2] )
           
             if ( mat[sel, 2] < threshold ) {
               final <- ini
               info <- rbind(info, c(0, -10, -10) )

             } else {
      
               info <- rbind(info, mat[sel, ] )
               mat <- mat[-sel, ,drop = FALSE] 
               dat <- as.data.frame( dataset[, -info[, 1] ] )
             }
 
           }
           
         }
	   
      }
      
      runtime <- proc.time() - tic		
      info <- info[ info[, 1] > 0, ]
      colnames(mat) <- c("Variables", "log.p-values", "statistic")
      res <- list(runtime = runtime, info = info, mat = mat, ci_test = ci_test, final = final ) 
    }
  
  }
  
  res
  
}  
     


 
  
    