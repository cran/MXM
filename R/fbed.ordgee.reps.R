fbed.ordgee.reps <- function(y, x, id, reps, univ = NULL, alpha = 0.05, wei = NULL, K = 0) { 
  
  dm <- dim(x)
  p <- dm[2]
  n <- dm[1]
  ind <- 1:p
  sig <- log(alpha)
  stat <- numeric(p)
  sela <- NULL
  card <- 0
  sa <- NULL
  pva <- NULL
  k <- length( unique(y) ) + 1
  
  ep <- Rfast::check_data(x)
  if ( sum(ep>0) > 0 )  x[, ep] <- rnorm( n * ep )
  
  runtime <- proc.time()
  
  if ( is.null(univ) ) {
    for ( i in ind ) {
      fit2 <- geepack::ordgee( y ~ reps + x[, i], id = id, weights = wei )
      mod <- summary(fit2)
      stat[i] <- mod[k, 3]
    }
    n.tests <- p
    pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
    univ <- list()
    univ$stat <- stat
    univ$pvalue <- pval
  } else {
    n.tests <- 0
    stat <- univ$stat
    pval <- univ$pvalue
  }  
  s <- which(pval < sig)
  
  if ( length(s) > 0 ) {
    sel <- which.min(pval)
    sela <- sel
    s <- s[ - which(s == sel) ]
    sa <- stat[sel]
    pva <- pval[sel]
    stat <- numeric(p )
    #########
    while ( sum(s>0) > 0 ) {
      nr <- k + length(sela)
      for ( i in ind[s] )  {
        fit2 <- geepack::ordgee( y ~ reps + x[, sela] + x[, i], id = id, weights = wei )
        mod <- summary(fit2)
        stat[i] <- mod[nr, 3]
      }
      n.tests <- n.tests + length( ind[s] ) 
      pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
      s <- which(pval < sig) 
      sel <- which.min(pval) * ( length(s) > 0 )
      sa <- c(sa, stat[sel]) 
      pva <- c(pva, pval[sel])
      sela <- c(sela, sel[sel>0])
      s <- s[ - which(s == sel) ]
      if (sel > 0)  stat <- rep(0, p)
    } ## end while ( sum(s > 0) > 0 )
    
    card <- sum(sela > 0)
    
    if (K == 1) {
      nr <- k + length(sela)
      for ( i in ind[-sela] )  {
        fit2 <- geepack::ordgee( y ~ reps + x[, sela] + x[, i], id = id, weights = wei )
        mod <- summary(fit2)
        stat[i] <- mod[nr, 3]
      }
      n.tests[2] <- length( ind[-sela] )
      pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
      s <- which(pval < sig)
      sel <- which.min(pval) * ( length(s)>0 )
      sa <- c(sa, stat[sel]) 
      pva <- c(pva, pval[sel])
      sela <- c(sela, sel[sel>0])
      s <- s[ - which(s == sel) ]
      if (sel > 0)   stat <- numeric(p)
      while ( sum(s>0) > 0 ) {
        nr <- k + length(sela)
        for ( i in ind[s] )  {
          fit2 <- geepack::ordgee( y ~ reps + x[, sela] + x[, i], id = id, weights = wei )
          mod <- summary(fit2)
          stat[i] <- mod[nr, 3]
        }
        n.tests[2] <- n.tests[2] + length( ind[s] )
        pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
        s <- which(pval < sig)
        sel <- which.min(pval) * ( length(s)>0 )
        sa <- c(sa, stat[sel]) 
        pva <- c(pva, pval[sel])
        sela <- c(sela, sel[sel>0])
        s <- s[ - which(s == sel) ]
        if (sel > 0)  stat <- numeric(p)
      } ## end while ( sum(s>0) > 0 ) 
      card <- c(card, sum(sela>0) )
    }  ## end if ( K == 1 ) 
    
    if ( K > 1) {
      nr <- k + length(sela)
      for ( i in ind[-sela] )  {
        fit2 <- geepack::ordgee( y ~ reps + x[, sela] + x[, i], id = id, weights = wei )
        mod <- summary(fit2)
        stat[i] <- mod[nr, 3]
      }
      n.tests[2] <- length( ind[-sela] )
      pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
      s <- which(pval < sig)
      sel <- which.min(pval) * ( length(s)>0 )
      sa <- c(sa, stat[sel]) 
      pva <- c(pva, pval[sel])
      sela <- c(sela, sel[sel>0])
      s <- s[ - which(s == sel) ]
      if (sel > 0)  stat <- numeric(p)
      
      while ( sum(s > 0) > 0 ) {
        nr <- k + length(sela)
        for ( i in ind[s] )  {
          fit2 <- geepack::ordgee( y ~ reps + x[, sela] + x[, i], id = id, weights = wei )
          mod <- summary(fit2)
          stat[i] <- mod[nr, 3]
        }
        n.tests[2] <- n.tests[2] + length( ind[s] )  
        pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
        s <- which(pval < sig)
        sel <- which.min(pval) * ( length(s)>0 )
        sa <- c(sa, stat[sel]) 
        pva <- c(pva, pval[sel])
        sela <- c(sela, sel[sel>0])
        s <- s[ - which(s == sel) ]
        if (sel > 0)  stat <- numeric(p)
      } ## end while ( sum(s>0) > 0 ) 
      
      card <- c(card, sum(sela > 0) )
      vim <- 1
      while ( vim < K  & card[vim + 1] - card[vim] > 0 ) {
        vim <- vim + 1
        nr <- k + length(sela)
        for ( i in ind[-sela] )  {
          fit2 <- geepack::ordgee( y ~ reps + x[, sela] + x[, i], id = id, weights = wei )
          mod <- summary(fit2)
          stat[i] <- mod[nr, 3]
        }
        n.tests[vim + 1] <- length( ind[-sela] )
        pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
        s <- which(pval < sig)
        sel <- which.min(pval) * ( length(s)>0 )
        sa <- c(sa, stat[sel]) 
        pva <- c(pva, pval[sel])
        sela <- c(sela, sel[sel>0])
        s <- s[ - which(s == sel) ]
        if (sel > 0)  stat <- numeric(p)
        while ( sum(s > 0) > 0 ) {
          nr <- k + length(sela)
          for ( i in ind[s] )  {
            fit2 <- geepack::ordgee( y ~ reps + x[, sela] + x[, i], id = id, weights = wei )
            mod <- summary(fit2)
            stat[i] <- mod[nr, 3]
          }
          n.tests[vim + 1] <- n.tests[vim + 1] + length( ind[s] )
          pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = TRUE)
          s <- which(pval < sig)
          sel <- which.min(pval) * ( length(s)>0 )
          sa <- c(sa, stat[sel]) 
          pva <- c(pva, pval[sel])
          sela <- c(sela, sel[sel>0])
          s <- s[ - which(s == sel) ]
          if (sel > 0)  stat <- numeric(p)
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
  list(res = res, info = info, runtime = runtime)
}
