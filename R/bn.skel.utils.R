bn.skel.utils <- function(mod, G = NULL, roc = TRUE, alpha = 0.01) {
  area <- NULL
  preds <- mod$pvalue
  preds <- preds[ upper.tri(preds) ]
  if ( !is.null(G) ) {
    group <- G[ upper.tri(G) ] 
    area <- auc(group,  -preds, roc = roc)
  }  
  p <- exp(preds)
  sig.p <- p[ which(p <= alpha) ]
  fdr <- length(p) / length(sig.p) * max(sig.p)
  theseis <- which( mod$pvalue < log(alpha), arr.ind = TRUE )
  sig.p <- cbind(theseis, preds[ preds < log(alpha) ] )
  sig.p <- sig.p[order(sig.p[, 3]), ]
  sig.p[, 3] <- exp(sig.p[, 3])
  theseis <- NULL
  list(area = area, fdr = fdr,sig.pvales = sig.p)
}