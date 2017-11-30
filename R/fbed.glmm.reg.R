fbed.glmm.reg <- function(y, x, id, alpha = 0.05, wei = NULL, K = 0, method = "LR", gam = NULL, backward = TRUE, type = "gaussian") {

  #check for NA values in the dataset and replace them with the variable median or the mode
  if ( any( is.na(x) ) ) {  
    #dataset = as.matrix(dataset);
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    x <- apply( x, 2, function(x){ x[which(is.na(x))] = median(x, na.rm = TRUE) ; return(x) } ) 
  }
  
  if ( method =="LR" ) {
    if (type == "gaussian") {
      result <- fbed.lmm(y = y, x = x, id = id, alpha = alpha, wei = wei, K = K)
    } else {
      result <- fbed.glmm(y = y, x = x, id = id, alpha = alpha, wei = wei, K = K, type = type)
    }
    
    result$back.rem <- 0
    result$back.n.tests <- 0
    
    if ( backward ) {
      
      if (result$info[1, 1] > 0) {
        a <- glmm.bsreg(y, x[, result$res[, 1], drop = FALSE], id, threshold = alpha, wei = wei, type = type)
        
        if ( typeof(a) == "list" ) {
          result$back.rem <- result$res[a$info[, 1], 1]
          back.n.tests <- sum( dim(result$res)[1] : dim(a$mat)[1] )
          sel <- result$res[a$mat[, 1], 1] 
          stat <- a$mat[, 3]
          pval <- a$mat[, 2]
          result$res <- cbind(sel, stat, pval)
          result$back.n.tests <- back.n.tests
          result$runtime <- result$runtime + a$runtime
        } else {
          back.rem <- 0
          back.n.tests <- 0
          result$back.rem <- back.rem
          result$back.n.tests <- back.n.tests
          result$runtime <- result$runtime 
        }  ## end if ( typeof(a) == "list" ) 
      }   ## end if (result$info[1, 1] > 0)
    }  ## end if ( backward )
    
  } else {
    
    if (type == "gaussian") {
      result <- ebic.fbed.lmm(y, x, id, gam = gam, wei = wei, K = K)
    } else {
      result <- ebic.fbed.glmm(y, x, id, gam = gam, wei = wei, K = K, type = type) 
    } 
      
    result$back.rem <- 0
    result$back.n.tests <- 0
    
    if ( backward ) {
      
      if (result$info[1, 1] > 0) {
        a <- ebic.glmm.bsreg(y, x[, result$res[, 1], drop = FALSE], id, wei = wei, gam = gam, type = type)

        if ( typeof(a) == "list" ) {
          back.n.tests <- sum( dim(result$res)[1] : length(a$mat[, 1]) )
          
          result$back.rem <- result$res[a$info[, 1], 1]
          sel <- result$res[ a$mat[, 1], 1]
          val <- a$mat[, 2]
          result$res <- cbind(sel, val)
          colnames(result$res) <- c("Vars", "eBIC")
          result$back.n.tests <- back.n.tests
          result$runtime <- result$runtime + a$runtime
        } else {
          back.rem <- 0
          back.n.tests <- 0
          result$back.rem <- back.rem
          result$back.n.tests <- back.n.tests
          result$runtime <- result$runtime 
        }  
      }   ## end if (result$info[1, 1] > 0) 
      
    }  ## end if ( backward )
    
  }  ## end if ( method == "LR" ) 
  
  result
}

