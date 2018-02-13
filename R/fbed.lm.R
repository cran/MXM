fbed.lm <- function(y, x, alpha = 0.05, univ = NULL, wei = NULL, K = 0) { 
  dm <- dim(x)
  n <- dm[1]
  p <- dm[2]
  ind <- 1:p
  sig <- log(alpha)
  lik1 <- Rfast::Var(y) * (n - 1)
  lik2 <- numeric(p)
  dof <- numeric(p)
  sela <- NULL
  card <- 0
  sa <- NULL
  pva <- NULL
  
  runtime <- proc.time()
  
  if ( is.null(univ) ) {
    if ( is.null(wei) ) {
      for ( i in ind ) {
        X <- model.matrix(y ~ x[, i])
        fit2 <- .lm.fit( X, y)
        lik2[i] <- sum(fit2$residuals^2)
        dof[i] <- length(fit2$coefficients)
      }
    } else {  
      for ( i in ind ) {
        X <- model.matrix(y ~ x[, i])
        fit2 <- lm.wfit(X, y, w = wei)
        lik2[i] <- sum(fit2$residuals^2)
        dof[i] <- length(fit2$coefficients)
      }
    }
    n.tests <- p
    stat <- (lik1 - lik2) * (n - dof) / ( lik2 * (dof - 1) )
    pval <- pf(stat, dof - 1, n - dof, lower.tail = FALSE, log.p = TRUE)
    univ$stat <- stat
    univ$pvalue <- pval
  } else {
    stat <- univ$stat
    pval <- univ$pvalue
    n.tests <- 0
    sel <- which.min(pval)
  }  
  s <- which(pval < sig)
  
  if ( length(s) > 0 ) {
    sel <- which.min(pval)
    sela <- sel
    s <- s[ - which(s == sel) ]
    sa <- stat[sel]
    pva <- pval[sel]
    pval <- numeric(p)
    #########
    while ( sum(s>0) > 0 ) {
      for ( i in ind[s] )  {
        fit2 <- lm( y ~., data = x[, c(sela, i)], weights = wei )
        if ( any( is.na(fit2$coefficients) ) )  {
          stat[i] <- 0
          pval <- 0
        } else {
          mod <- anova(fit2)
          pr <- dim(mod)[1] - 1
          df1 <- mod[pr, 1]
          df2 <- mod[pr + 1, 1]
          stat[i] <- mod[pr, 4]
          pval[i] <- pf(mod[pr, 4], df1, df2, lower.tail = FALSE, log.p = TRUE)
        }
      }
      n.tests <- n.tests + length( ind[s] ) 
      s <- which(pval < sig) 
      sel <- which.min(pval) * ( length(s)>0 )
      sa <- c(sa, stat[sel]) 
      pva <- c(pva, pval[sel])
      sela <- c(sela, sel[sel > 0] )
      s <- s[ - which(s == sel) ]
      if (sel > 0) {
        stat <- numeric(p)
        pval <- numeric(p)
      } 
    } ## end while ( sum(s > 0) > 0 )
    
  card <- sum(sela > 0)
  
  if (K == 1) {
    for ( i in ind[-sela] )  {
      fit2 <- lm( y ~., data = x[, c(sela, i)], weights = wei )
      if ( any( is.na(fit2$coefficients) ) )  {
        stat[i] <- 0
        pval <- 0
      } else {
        mod <- anova(fit2)
        pr <- dim(mod)[1] - 1
        df1 <- mod[pr, 1]
        df2 <- mod[pr + 1, 1]
        stat[i] <- mod[pr, 4]
        pval[i] <- pf(mod[pr, 4], df1, df2, lower.tail = FALSE, log.p = TRUE)
      }
    }
    n.tests[2] <- length( ind[-sela] )
    s <- which(pval < sig)
    sel <- which.min(pval) * ( length(s)>0 )
    sa <- c(sa, stat[sel]) 
    pva <- c(pva, pval[sel])
    sela <- c(sela, sel[sel > 0] )
    s <- s[ - which(s == sel) ]
    pval <- numeric(p)
    while ( sum(s>0) > 0 ) {
      for ( i in ind[s] )  {
        fit2 <- lm( y ~., data = x[, c(sela, i)], weights = wei )
        if ( any( is.na(fit2$coefficients) ) )  {
          stat[i] <- 0
          pval <- 0
        } else {
          mod <- anova(fit2)
          pr <- dim(mod)[1] - 1
          df1 <- mod[pr, 1]
          df2 <- mod[pr + 1, 1]
          stat[i] <- mod[pr, 4]
          pval[i] <- pf(mod[pr, 4], df1, df2, lower.tail = FALSE, log.p = TRUE)
        }
      }
      n.tests[2] <- n.tests[2] + length( ind[s] )
      s <- which(pval < sig)
      sel <- which.min(pval) * ( length(s)>0 )
      sa <- c(sa, stat[sel]) 
      pva <- c(pva, pval[sel])
      sela <- c(sela, sel[sel>0])
      s <- s[ - which(s == sel) ]
      if (sel > 0) {
        stat <- numeric(p)
        pval <- numeric(p)
      } 
    } ## end while ( sum(s>0) > 0 ) 
    card <- c(card, sum(sela>0) )
  }  ## end if ( K == 1 ) 
  
  if ( K > 1) {
    
    for ( i in ind[-sela] )  {
      fit2 <- lm( y ~., data = x[, c(sela, i)], weights = wei )
      if ( any( is.na(fit2$coefficients) ) )  {
        stat[i] <- 0
        pval <- 0
      } else {
        mod <- anova(fit2)
        pr <- dim(mod)[1] - 1
        df1 <- mod[pr, 1]
        df2 <- mod[pr + 1, 1]
        stat[i] <- mod[pr, 4]
        pval[i] <- pf(mod[pr, 4], df1, df2, lower.tail = FALSE, log.p = TRUE)
      }
    }
    n.tests[2] <- length( ind[-sela] ) 
    s <- which(pval < sig)
    sel <- which.min(pval) * ( length(s)>0 )
    sa <- c(sa, stat[sel]) 
    pva <- c(pva, pval[sel])
    sela <- c(sela, sel[sel>0])
    s <- s[ - which(s == sel) ]
    pval <- numeric(p)
    while ( sum(s > 0) > 0 ) {
      for ( i in ind[s] )  {
        fit2 <- lm( y ~., data = x[, c(sela, i)], weights = wei )
        if ( any( is.na(fit2$coefficients) ) )  {
          stat[i] <- 0
          pval <- 0
        } else {
          mod <- anova(fit2)
          pr <- dim(mod)[1] - 1
          df1 <- mod[pr, 1]
          df2 <- mod[pr + 1, 1]
          stat[i] <- mod[pr, 4]
          pval[i] <- pf(mod[pr, 4], df1, df2, lower.tail = FALSE, log.p = TRUE)
        }
      }
      n.tests[2] <- n.tests[2] + length( ind[s] )  
      s <- which(pval < sig)
      sel <- which.min(pval) * ( length(s)>0 )
      sa <- c(sa, stat[sel]) 
      pva <- c(pva, pval[sel])
      sela <- c(sela, sel[sel>0])
      s <- s[ - which(s == sel) ]
      if (sel > 0) {
        stat <- numeric(p)
        pval <- numeric(p)
      } 
    } ## end while ( sum(s>0) > 0 ) 
    
    card <- c(card, sum(sela > 0) )
    vim <- 1
    while ( vim < K  & card[vim + 1] - card[vim] > 0 ) {
      vim <- vim + 1
      for ( i in ind[-sela] )  {
        fit2 <- lm( y ~., data = x[, c(sela, i)], weights = wei )
        if ( any( is.na(fit2$coefficients) ) )  {
          stat[i] <- 0
          pval <- 0
        } else {
          mod <- anova(fit2)
          pr <- dim(mod)[1] - 1
          df1 <- mod[pr, 1]
          df2 <- mod[pr + 1, 1]
          stat[i] <- mod[pr, 4]
          pval[i] <- pf(mod[pr, 4], df1, df2, lower.tail = FALSE, log.p = TRUE)
        }
      }
      n.tests[vim + 1] <- length( ind[-sela] )
      s <- which(pval < sig)
      sel <- which.min(pval) * ( length(s)>0 )
      sa <- c(sa, stat[sel]) 
      pva <- c(pva, pval[sel])
      sela <- c(sela, sel[sel>0])
      s <- s[ - which(s == sel) ]
      pval <- numeric(p)
      while ( sum(s > 0) > 0 ) {
        for ( i in ind[s] )  {
          fit2 <- lm( y ~., data = x[, c(sela, i)], weights = wei )
          if ( any( is.na(fit2$coefficients) ) )  {
            stat[i] <- 0
            pval <- 0
          } else {
            mod <- anova(fit2)
            pr <- dim(mod)[1] - 1
            df1 <- mod[pr, 1]
            df2 <- mod[pr + 1, 1]
            stat[i] <- mod[pr, 4]
            pval[i] <- pf(mod[pr, 4], df1, df2, lower.tail = FALSE, log.p = TRUE)
          }
        }
        n.tests[vim + 1] <- n.tests[vim + 1] + length( ind[s] )
        s <- which(pval < sig)
        sel <- which.min(pval) * ( length(s)>0 )
        sa <- c(sa, stat[sel]) 
        pva <- c(pva, pval[sel])
        sela <- c(sela, sel[sel>0])
        s <- s[ - which(s == sel) ]
        if (sel > 0) {
          stat <- numeric(p)
          pval <- numeric(p)
        } 
      } ## end while ( sum(s > 0) > 0 ) 
      card <- c(card, sum(sela>0) )
    }  ## end while ( vim < K )
  } ## end if ( K > 1)
  } ## end if ( length(s) > 0 )

  runtime <- proc.time() - runtime
  len <- sum( sela > 0 )
  if (len > 0) {
    res <- cbind(sela[1:len], sa[1:len], pva[1:len] )
    info <- matrix(nrow = length(card), ncol = 2)
    info[, 1] <- card
    info[, 2] <- n.tests
  } else {
    res <- matrix(c(0, 0, 0), ncol = 3)
    info <- matrix(c(0, p), ncol = 2)
  }  
  colnames(res) <- c("Vars", "stat", "log p-value")
  rownames(info) <- paste("K=", 1:length(card)- 1, sep = "")
  colnames(info) <- c("Number of vars", "Number of tests")
  list(univ = univ, res = res, info = info, runtime = runtime)
}
