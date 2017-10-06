######## Internal for Poisson regression
internaliamb.poisbs <- function(target, dataset, threshold, wei, p, heavy = FALSE, robust = FALSE) {
  if ( !is.null(dataset) |  p > 0 ) {
    if ( p > 1 ) {
      if ( !heavy ) {
        ini <- glm( target ~.,  data = data.frame(dataset), family = poisson(log), weights = wei, y = FALSE, model = FALSE )
        tab <- drop1( ini, test = "Chisq" )
        dof <- tab[-1, 1]
        stat <- tab[-1, 4]
        
      } else {
        ini <- speedglm::speedglm( target ~.,  data = data.frame(dataset), family = poisson(log), weights = wei )
        dofini <- length( coef(ini) )
        stat <- dof <- numeric(p)
        for (i in 1:p) {
          mod <- speedglm::speedglm( target ~.,  data = data.frame(dataset[, -i]), family = poisson(log), weights = wei )
          stat[i] <- mod$deviance - ini$deviance
          dof[i] <- dofini - length( coef(mod) ) 
        }
      }	
    } else {
      if ( !heavy ) {
        ini <- glm( target ~.,  data = data.frame(dataset), family = poisson(log), weights = wei, y = FALSE, model = FALSE )
        mod0 <- glm(target ~ 1, family = poisson(log), weights = wei, y = FALSE, model = FALSE)
      } else {
        ini <- speedglm::speedglm( target ~.,  data = data.frame(dataset), family = poisson(log), weights = wei )
        mod0 <- speedglm::speedglm( target ~ 1,  data = data.frame(dataset), family = poisson(log), weights = wei )
      }  
      stat <- mod0$deviance - ini$deviance
      dof <- length( coef(ini) ) - 1
    }
    
    mat <- cbind(1:p, pchisq( stat, dof, lower.tail = FALSE, log.p = TRUE), stat )
    colnames(mat) <- c("variable", "log.p-values", "statistic" )
    rownames(mat) <- 1:p  
    
    info <- matrix( c(0, -10, -10) , ncol = 3 )
    sel <- which( mat[, 2] > threshold ) 
    
    if ( length(sel) == 0 ) {
      final <- ini 
      
    } else {
      
      info <- mat[sel, , drop = FALSE]
      mat <- mat[-sel, , drop = FALSE] 
      dat <- as.data.frame( dataset[, -sel] ) 
      
      if ( p - length(sel) == 0 ) {
        final <- "No variables were selected"
        mat <- NULL
      } else if ( p - length(sel) == 1 ) {
        if ( !heavy ) {
          mod1 <- glm(target ~., data = data.frame(dat), family = poisson(log), weights = wei, y = FALSE, model = FALSE)
          mod0 <- glm(target ~ 1, family = poisson(log), weights = wei, y = FALSE, model = FALSE)
        } else  {
          mod1 <- speedglm::speedglm( target ~ .,  data = data.frame(dat), family = poisson(log), weights = wei )
          mod0 <- speedglm::speedglm( target ~ 1,  data = data.frame(dat), family = poisson(log), weights = wei )
        }  
        stat <- abs( mod1$deviance - mod0$deviance )
        pval <- pchisq( stat, length( coef(mod1) ) - 1, lower.tail = FALSE, log.p = TRUE)
        if (pval > threshold ) {
          final <- "No variables were selected"
          mat <- NULL
        } else final <- mod1
      } else {
        if ( !heavy ) {
          final <- glm(target ~., data = data.frame(dat), family = poisson(log), weights = wei, y = FALSE, model = FALSE)
        } else  final <- speedglm::speedglm( target ~ .,  data = data.frame(dat), family = poisson(log), weights = wei )
      }
    }
    info <- info[ info[, 1] > 0, , drop = FALSE]
    
  } else { 
    info <- NULL  
    mat <- NULL
    final <- "No variables were selected"
  } 
  
  list(info = info, mat = mat, final = final) 
}  