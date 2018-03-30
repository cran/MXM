fbed.clogit <- function(y, x, alpha = 0.05, univ = NULL, wei = NULL, K = 0) { 
  sig <- log(alpha)
  dm <- dim(x)
  p <- dm[2]
  ind <- 1:p
  lik2 <- numeric(p)
  dof <- numeric(p)
  sela <- NULL
  card <- 0
  sa <- NULL
  pva <- NULL
  runtime <- proc.time()
  id <- y[, 2]
  case <- y[, 1]
  
  lik1 <- survival::clogit( case ~ x[, 1] + strata(id) )$loglik[1]
  
  if ( is.null(univ) ) {
    for ( i in ind ) {
      fit2 <- survival::clogit( case ~ x[, i] + strata(id) )
      lik2[i] <- logLik(fit2)
      dof[i] <- length(fit2$coefficients)
    }
    stat <- 2 * (lik2 - lik1)
    pval <- pchisq(stat, dof, lower.tail = FALSE, log.p = TRUE)
    univ$stat <- stat
    univ$pvalue <- pval
    n.tests <- p
  } else {
    n.tests <- 0
    stat <- univ$stat
    pval <- univ$pvalue
    lik2 <- univ$stat + 2 * lik1
  } 
  
  s <- which(pval < sig)
  
  if ( length(s) > 0 ) {
    sel <- which.min(pval)
    sela <- sel
    s <- s[ - which(s == sel) ]
    lik1 <- lik2[sel] 
    d1 <- dim( model.matrix( y~., data.frame(x[, sel]) ) )[2] - 1
    sa <- stat[sel]
    pva <- pval[sel]
    lik2 <- rep( lik1, p )
    dof <- rep(10000, p)
    #########
    while ( sum(s>0) > 0 ) {
      for ( i in ind[s] )  {
        fit2 <- try( survival::clogit( case ~. + strata(id), data = x[, c(sela, i)] ), silent = TRUE)
          if ( identical( class(fit2), "try-error" ) ) {
          lik2[i] <- lik1
        } else {
          lik2[i] <- logLik(fit2)  
          dof[i] <- length( fit2$coefficients )          
        } 
      } ## end  for
      n.tests <- n.tests + length( ind[s] ) 
      stat <- 2 * (lik2 - lik1)
      pval <- pchisq(stat, dof - d1, lower.tail = FALSE, log.p = TRUE)
      s <- which(pval < sig) 
      sel <- which.min(pval) * ( length(s)>0 )
      sa <- c(sa, stat[sel]) 
      pva <- c(pva, pval[sel])
      sela <- c(sela, sel[sel>0])
      s <- s[ - which(s == sel) ]
      if (sel > 0) {
        lik1 <- lik2[sel] 
        d1 <- dof[sel]
        lik2 <- rep(lik1, p)
        dof <- rep(10000, p)
      }  
    } ## end while ( sum(s > 0) > 0 )
    
    card <- sum(sela > 0)
    
    if (K == 1) {
      for ( i in ind[-sela] )  {
        fit2 <- try( survival::clogit( case ~. + strata(id), data = x[, c(sela, i)] ), silent = TRUE)
          if ( identical( class(fit2), "try-error" ) ) {
          lik2[i] <- lik1
        } else {
          lik2[i] <- logLik(fit2)  
          dof[i] <- length( fit2$coefficients )          
        } 
      }
      n.tests[2] <- length( ind[-sela] )
      stat <- 2 * (lik2 - lik1)
      pval <- pchisq(stat, dof - d1, lower.tail = FALSE, log.p = TRUE)
      s <- which(pval < sig)
      sel <- which.min(pval) * ( length(s)>0 )
      sa <- c(sa, stat[sel]) 
      pva <- c(pva, pval[sel])
      sela <- c(sela, sel[sel>0])
      s <- s[ - which(s == sel) ]
      if (sel > 0) {
        lik1 <- lik2[sel] 
        d1 <- dof[sel]
        lik2 <- rep(lik1, p)
        dof <- rep(10000, p)
      }  
      while ( sum(s>0) > 0 ) {
        for ( i in ind[s] )  {
          fit2 <- try( survival::clogit( case ~. + strata(id), data = x[, c(sela, i)] ), silent = TRUE)
          if ( identical( class(fit2), "try-error" ) ) {
            lik2[i] <- lik1
          } else {
            lik2[i] <- logLik(fit2)  
            dof[i] <- length( fit2$coefficients )          
          } 
        }
        n.tests[2] <- n.tests[2] + length( ind[s] )
        stat <- 2 * (lik2 - lik1)
        pval <- pchisq(stat, dof - d1, lower.tail = FALSE, log.p = TRUE)
        s <- which(pval < sig)
        sel <- which.min(pval) * ( length(s)>0 )
        sa <- c(sa, stat[sel]) 
        pva <- c(pva, pval[sel])
        sela <- c(sela, sel[sel>0])
        s <- s[ - which(s == sel) ]
        if (sel > 0) {
          lik1 <- lik2[sel] 
          d1 <- dof[sel]
          lik2 <- rep(lik1, p)
          dof <- rep(10000, p)
        }  
      } ## end while ( sum(s>0) > 0 ) 
      card <- c(card, sum(sela>0) )
    }  ## end if ( K == 1 ) 
    
    if ( K > 1) {
      
      for ( i in ind[-sela] )  {
        fit2 <- try( survival::clogit( case ~. + strata(id), data = x[, c(sela, i)] ), silent = TRUE)
          if ( identical( class(fit2), "try-error" ) ) {
          lik2[i] <- lik1
        } else {
          lik2[i] <- logLik(fit2)  
          dof[i] <- length( fit2$coefficients )          
        } 
      }
      n.tests[2] <- length( ind[-sela] ) 
      stat <- 2 * (lik2 - lik1)
      pval <- pchisq(stat, dof - d1, lower.tail = FALSE, log.p = TRUE)
      s <- which(pval < sig)
      sel <- which.min(pval) * ( length(s)>0 )
      sa <- c(sa, stat[sel]) 
      pva <- c(pva, pval[sel])
      sela <- c(sela, sel[sel>0])
      s <- s[ - which(s == sel) ]
      if (sel > 0) {
        lik1 <- lik2[sel] 
        d1 <- dof[sel]
        lik2 <- rep(lik1, p)
        dof <- rep(10000, p)
      }  
      while ( sum(s > 0) > 0 ) {
        for ( i in ind[s] )  {
          fit2 <- try( survival::clogit( case ~. + strata(id), data = x[, c(sela, i)] ), silent = TRUE)
          if ( identical( class(fit2), "try-error" ) ) {
            lik2[i] <- lik1
          } else {
            lik2[i] <- logLik(fit2)  
            dof[i] <- length( fit2$coefficients )          
          } 
        }
        n.tests[2] <- n.tests[2] + length( ind[s] )  
        stat <- 2 * (lik2 - lik1)
        pval <- pchisq(stat, dof - d1, lower.tail = FALSE, log.p = TRUE)
        s <- which(pval < sig)
        sel <- which.min(pval) * ( length(s)>0 )
        sa <- c(sa, stat[sel]) 
        pva <- c(pva, pval[sel])
        sela <- c(sela, sel[sel>0])
        s <- s[ - which(s == sel) ]
        if (sel > 0) {
          lik1 <- lik2[sel] 
          d1 <- dof[sel]
          lik2 <- rep(lik1, p)
          dof <- rep(10000, p)
        }  
      } ## end while ( sum(s>0) > 0 ) 
      
      card <- c(card, sum(sela > 0) )
      vim <- 1
      while ( vim < K  & card[vim + 1] - card[vim] > 0 ) {
        vim <- vim + 1
        for ( i in ind[-sela] )  {
          fit2 <- try( survival::clogit( case ~. + strata(id), data = x[, c(sela, i)] ), silent = TRUE)
          if ( identical( class(fit2), "try-error" ) ) {
            lik2[i] <- lik1
          } else {
            lik2[i] <- logLik(fit2)  
            dof[i] <- length( fit2$coefficients )          
          } 
        }
        n.tests[vim + 1] <- length( ind[-sela] )
        stat <- 2 * (lik2 - lik1)
        pval <- pchisq(stat, dof - d1, lower.tail = FALSE, log.p = TRUE)
        s <- which(pval < sig)
        sel <- which.min(pval) * ( length(s)>0 )
        sa <- c(sa, stat[sel]) 
        pva <- c(pva, pval[sel])
        sela <- c(sela, sel[sel>0])
        s <- s[ - which(s == sel) ]
        if (sel > 0) {
          lik1 <- lik2[sel] 
          d1 <- dof[sel]
          lik2 <- rep(lik1, p)
          dof <- rep(10000, p)
        }     
        while ( sum(s > 0) > 0 ) {
          for ( i in ind[s] )  {
            fit2 <- try( survival::clogit( case ~. + strata(id), data = x[, c(sela, i)] ), silent = TRUE)
            if ( identical( class(fit2), "try-error" ) ) {
              lik2[i] <- lik1
            } else {
              lik2[i] <- logLik(fit2)  
              dof[i] <- length( fit2$coefficients )          
            } 
          }
          n.tests[vim + 1] <- n.tests[vim + 1] + length( ind[s] )
          stat <- 2 * (lik2 - lik1)
          pval <- pchisq(stat, dof - d1, lower.tail = FALSE, log.p = TRUE)
          s <- which(pval < sig)
          sel <- which.min(pval) * ( length(s)>0 )
          sa <- c(sa, stat[sel]) 
          pva <- c(pva, pval[sel])
          sela <- c(sela, sel[sel>0])
          s <- s[ - which(s == sel) ]
          if (sel > 0) {
            lik1 <- lik2[sel] 
            d1 <- dof[sel]
            lik2 <- rep(lik1, p)
            dof <- rep(10000, p)
          }  
        } ## end while ( sum(s > 0) > 0 ) 
        card <- c(card, sum(sela>0) )
      }  ## end while ( vim < K )
    } ## end if ( K > 1)
  } ## end if ( length(s) > 0 )
  
  runtime <- proc.time() - runtime
  len <- sum( sela > 0 )
  if (len > 0) {
    res <- cbind(sela[1:len], 2 * sa[1:len], pva[1:len] )
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
