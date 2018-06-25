big.fbed.g2 <- function(z, threshold = 0.01, univ = NULL, K = 0, backward = TRUE) { 
  dm <- dim(z)
  n <- dm[1]    ;    p <- dim(z)[2]  - 1
  ind <- 1:p
  sig <- log(threshold)
  sela <- NULL
  card <- 0
  sa <- NULL
  pva <- NULL
  runtime <- proc.time()
  
  m <- colSums(z[]) 
  x2 <- colSums(z[]^2)
  s <- (x2 - m^2/n)
  s <- s[-(p + 1)]
  ind[s == 0] <- 0

  all.dc <- numeric(p)
  for (i in 1:(p + 1) )  all.dc[i] <- max(z[, i]) + 1

  if ( is.null(univ) ) {
    y <- z[, p + 1]
    stat <- pval <- numeric(p)
    for (i in ind) {
      a <- Rfast::Table(z[, i], y, names = FALSE)
      dof <- prod(dim(a) - 1)
      rs <- Rfast::rowsums(a)
      cs <- Rfast::colsums(a)
      est <- outer(rs, cs, "*")/sum(a)
      stat[i] <- 2 * sum(a * log(a/est), na.rm = TRUE)
      pval[i] <- pchisq(stat[i], dof, lower.tail = FALSE, log.p = TRUE)
    }
    univ$stat <- stat
    univ$pvalue <- pval
    n.tests <- p
  } else {
    stat <- univ$stat
    pval <- univ$pvalue
    n.tests <- 0
  }
  
  s <- which(pval < sig)
  
  if ( length(s) > 0 ) {
    sel <- which.min(pval)
    sela <- sel
    s <- s[ - which(s == sel) ]
    sa <- stat[sel]
    pva <- pval[sel]
    stat <- numeric(p)
    pval <- numeric(p)
    #########
    while ( sum(s>0) > 0 ) {
      for ( i in ind[s] )  {
        k <- length(sela)
        dc <- all.dc[c(p + 1, i, sela)]
        mod <- Rfast::g2Test(z[, c(p + 1, i, sela)], x = 1, y = 2, cs = 3:(2 + k), dc = dc)           
        stat[i] <- mod$statistic
        pval[i] <- pchisq(mod$statistic, mod$df, lower.tail = FALSE, log.p = TRUE)
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
      stat <- numeric(p)
      pval <- numeric(p)
      for ( i in ind[-sela] )  {
        k <- length(sela)
        dc <- all.dc[c(p + 1, i, sela)]
        mod <- Rfast::g2Test(z[, c(p + 1, i, sela)], x = 1, y = 2, cs = 3:(2 + k), dc = dc)           
        stat[i] <- mod$statistic
        pval[i] <- pchisq(mod$statistic, mod$df, lower.tail = FALSE, log.p = TRUE)
      }
      n.tests[2] <- length( ind[-sela] )
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
      while ( sum(s>0) > 0 ) {
        for ( i in ind[s] )  {
          k <- length(sela)
          dc <- all.dc[c(p + 1, i, sela)]
          mod <- Rfast::g2Test(z[, c(p + 1, i, sela)], x = 1, y = 2, cs = 3:(2 + k), dc = dc)           
          stat[i] <- mod$statistic
          pval[i] <- pchisq(mod$statistic, mod$df, lower.tail = FALSE, log.p = TRUE)
        }
        n.tests[2] <- n.tests[2] + length( ind[s] )
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
      } ## end while ( sum(s>0) > 0 ) 
      card <- c(card, sum(sela>0) )
    }  ## end if ( K == 1 ) 
    
    if ( K > 1) {
      
      for ( i in ind[-sela] )  {
        k <- length(sela)
        dc <- all.dc[c(p + 1, i, sela)]
        mod <- Rfast::g2Test(z[, c(p + 1, i, sela)], x = 1, y = 2, cs = 3:(2 + k), dc = dc)           
        stat[i] <- mod$statistic
        pval[i] <- pchisq(mod$statistic, mod$df, lower.tail = FALSE, log.p = TRUE)
      }
      n.tests[2] <- length( ind[-sela] ) 
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
      while ( sum(s > 0) > 0 ) {
        for ( i in ind[s] )  {
          k <- length(sela)
          dc <- all.dc[c(p + 1, i, sela)]
          mod <- Rfast::g2Test(z[, c(p + 1, i, sela)], x = 1, y = 2, cs = 3:(2 + k), dc = dc)           
          stat[i] <- mod$statistic
          pval[i] <- pchisq(mod$statistic, mod$df, lower.tail = FALSE, log.p = TRUE)
        }
        
        n.tests[2] <- n.tests[2] + length( ind[s] )  
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
      } ## end while ( sum(s>0) > 0 ) 
      
      card <- c(card, sum(sela > 0) )
      vim <- 1
      while ( vim < K  & card[vim + 1] - card[vim] > 0 ) {
        vim <- vim + 1
        for ( i in ind[-sela] )  {
          k <- length(sela)
          dc <- all.dc[c(p + 1, i, sela)]
          mod <- Rfast::g2Test(z[, c(p + 1, i, sela)], x = 1, y = 2, cs = 3:(2 + k), dc = dc)           
          stat[i] <- mod$statistic
          pval[i] <- pchisq(mod$statistic, mod$df, lower.tail = FALSE, log.p = TRUE)
        }
        n.tests[vim + 1] <- length( ind[-sela] )
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
        while ( sum(s > 0) > 0 ) {
          for ( i in ind[s] )  {
            k <- length(sela)
            dc <- all.dc[c(p + 1, i, sela)]
            mod <- Rfast::g2Test(z[, c(p + 1, i, sela)], x = 1, y = 2, cs = 3:(2 + k), dc = dc)           
            stat[i] <- mod$statistic
            pval[i] <- pchisq(mod$statistic, mod$df, lower.tail = FALSE, log.p = TRUE)
          }
          n.tests[vim + 1] <- n.tests[vim + 1] + length( ind[s] )
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
  result <- list(univ = univ, res = res, info = info, runtime = runtime)
  
  result$back.rem <- 0
  result$back.n.tests <- 0

  result
}
