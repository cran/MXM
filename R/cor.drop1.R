cor.drop1 <- function(y, x, logged = FALSE) {
  if ( is.matrix(x) )    x <- data.frame(x)
  mod <- lm(y ~., x)
  a <- drop1(mod, test = "F")
  dof <- summary(mod)[[10]][3]
  r <- sqrt( a[-1, 5] / (a[-1, 5] + dof) )
  stat <- abs( 0.5 * log( (1 + r) / (1 - r) ) ) * sqrt(dof - 1) 
  if ( logged ) {
    pval <- log(2) + pt(stat, dof - 1, lower.tail = FALSE, log.p = TRUE)
  } else pval <- 2 * pt(stat, dof - 1, lower.tail = FALSE)
  cbind(stat, pval)
}