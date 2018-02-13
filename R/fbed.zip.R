fbed.zip <- function(y, x, alpha = 0.05, univ = NULL, wei = NULL, K = 0) { 
  dm <- dim(x)
  n <- dm[1]
  p <- dm[2]
  ind <- 1:p
  sig <- log(alpha)
  lik2 <- numeric(p)
  sela <- NULL
  card <- 0
  sa <- NULL
  pva <- NULL
  
  if ( is.null(wei) ) {
    mo <- 2 * Rfast::zip.mle(y)$loglik
    lgy <- sum( lgamma(y + 1) )  
  } else  {
    mo <- 2 * zipmle.wei(y, wei)$loglik    
    lgy <- sum( wei * lgamma(y + 1) )  
  }
  runtime <- proc.time()
  
  if ( is.null(univ) ) {
    mod <- zip.regs(y, x, wei, logged = TRUE)
    n.tests <- p
    stat <- mod[, 1]
    pval <- mod[, 2]
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
      lik1 <- stat[sel] + mo 
      d1 <- dim( model.matrix( y~., data.frame(x[, sel]) ) )[2] - 1
      sa <- stat[sel]
      pva <- pval[sel]
      lik2 <- rep( lik1, p )
      dof <- rep(10000, p)
      #########
      while ( sum(s>0) > 0 ) {
        for ( i in ind[s] )  {
          fit2 <- zip.reg( y, x[, c(sela, i)], wei = wei, lgy = lgy )
          lik2[i] <- 2 * fit2$loglik
          dof[i] <- length( fit2$be )
        }
        n.tests <- n.tests + length( ind[s] ) 
        stat <- lik2 - lik1
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
        fit2 <- zip.reg( y, x[, c(sela, i)], wei = wei, lgy = lgy )
        lik2[i] <- 2 * fit2$loglik
        dof[i] <- length( fit2$be )
      }
      n.tests[2] <- length( ind[-sela] )
      stat <- lik2 - lik1
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
          fit2 <- zip.reg( y, x[, c(sela, i)], wei = wei, lgy = lgy )
          lik2[i] <- 2 * fit2$loglik
          dof[i] <- length( fit2$be )
        }
        n.tests[2] <- n.tests[2] + length( ind[s] )
        stat <- lik2 - lik1
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
        fit2 <- zip.reg( y, x[, c(sela, i)], wei = wei, lgy = lgy )
        lik2[i] <- 2 * fit2$loglik
        dof[i] <- length( fit2$be )
      }
      n.tests[2] <- length( ind[-sela] ) 
      stat <- lik2 - lik1
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
          fit2 <- zip.reg( y, x[, c(sela, i)], wei = wei, lgy = lgy )
          lik2[i] <- 2 * fit2$loglik
          dof[i] <- length( fit2$be )
        }
        n.tests[2] <- n.tests[2] + length( ind[s] )  
        stat <- lik2 - lik1
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
        fit2 = zip.reg( y, x[, c(sela, i)], wei = wei, lgy = lgy )
        lik2[i] <- 2 * fit2$loglik
        dof[i] <- length( fit2$be )
      }
      n.tests[vim + 1] <- length( ind[-sela] )
      stat <- lik2 - lik1
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
          fit2 <- zip.reg( y, x[, c(sela, i)], wei = wei, lgy = lgy )
          lik2[i] <- 2 * fit2$loglik
          dof[i] <- length( fit2$be )
        }
        n.tests[vim + 1] <- n.tests[vim + 1] + length( ind[s] )
        stat <- lik2 - lik1
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
 