\name{Correlation based tests}
\alias{testIndFisher}
\alias{testIndSpearman}
\title{
Fisher and Spearman conditional independence test for continuous class variables
}

\description{
The main task of this test is to provide a p-value PVALUE for the null hypothesis: feature 'X' is independent from 'TARGET' given a conditioning set CS.
}

\usage{
testIndFisher(target, dataset, xIndex, csIndex, dataInfo = NULL, univariateModels = NULL,
hash = FALSE, stat_hash = NULL, pvalue_hash = NULL, robust = FALSE)

testIndSpearman(target, dataset, xIndex, csIndex,dataInfo = NULL,univariateModels = NULL,
hash = FALSE, stat_hash = NULL, pvalue_hash = NULL, robust = FALSE)
}

\arguments{
  \item{target}{
  A numeric vector containing the values of the target variable. If the values are proportions or percentages, i.e. strictly within 0 and 1 they are mapped into R using log( target/(1 - target) ).
}
  \item{dataset}{
  A numeric matrix containing the variables for performing the test. Rows as samples and columns as features.
}
  \item{xIndex}{
  The index of the variable whose association with the target we want to test.
}
  \item{csIndex}{
  The indices of the variables to condition on.
}
  \item{dataInfo}{
  list object with information on the structure of the data. Default value is NULL.
}
  \item{univariateModels}{
  Fast alternative to the hash object for univariate test. List with vectors "pvalues" (p-values), "stats" (statistics) and "flags" (flag = TRUE if the test was succesful) representing the univariate association of each variable with the target. Default value is NULL.
}
  \item{hash}{
A boolean variable which indicates whether (TRUE) or not (FALSE) to use the hash-based implementation of the statistics of SES. Default value is FALSE. If TRUE you have to specify the stat_hash argument and the pvalue_hash argument.
}
  \item{stat_hash}{
A hash object (hash package required) which contains the cached generated statistics of a SES run in the current dataset, using the current test.
}
  \item{pvalue_hash}{
A hash object (hash package required) which contains the cached generated p-values of a SES run in the current dataset, using the current test.
}
  \item{robust}{
A boolean variable which indicates whether (TRUE) or not (FALSE) to use a robustified version of Fisher's correlation coefficient via M-estimation available from the "rlm" in MASS. Two regressions are fitted and the square root ot the absolute value of the beta coefficients is used to calculate the correlation coefficient (Shevlyakov and Smirnov, 2011). For the conditional correlation the correlation of the residuals of the two robust regressions is calcualted. For more ways of calculating the correlation coefficient see the references. It takes more time than non robust version but it is suggested in case of outliers. Default value is FALSE.

In the case of testIndSpearman, this is not used, as Spearman correlation is robust by default.
}
}

\details{
If hash = TRUE, testIndFisher requires the arguments 'stat_hash' and 'pvalue_hash' for the hash-based implementation of the statistic test. These hash Objects are produced or updated by each run of SES (if hash == TRUE) and they can be reused in order to speed up next runs of the current statistic test. If "SESoutput" is the output of a SES run, then these objects can be retrieved by SESoutput@hashObject$stat_hash and the SESoutput@hashObject$pvalue_hash.

Important: Use these arguments only with the same dataset that was used at initialization.

For all the available conditional independence tests that are currently included on the package, please see "?CondIndTests".

Note that if the \code{\link{testIndReg}} is used instead the results will not be be the same, unless the sample size is very large. This is because the Fisher test uses the t distribution stemming from the Fisher's z transform and not the t distribution of the correlation coefficient.

BE CAREFUL with testIndSpearman. The Pearson's correlation coefficient is actually calculated. So, you must have transformed the data into their ranks before plugging them here. The reason for this is to speed up the computation time, as this test can be used in SES, MMPC and mmhc.skel. The variance of the Fisher transformed Spearman's correlation is \eqn{\frac{1.06}{n-3}} and the variance of the Fisher transformed Pearson's correlation coefficient is \eqn{\frac{1}{n-3}}.
}

\value{
A list including:
\item{pvalue}{
A numeric value that represents the logarithm of the generated p-value due to Fisher's method (see reference below).
}
\item{stat}{
A numeric value that represents the generated statistic due to Fisher's method (see reference below).
}
\item{flag}{
A numeric value (control flag) which indicates whether the test was succesful (0) or not (1).
}
\item{stat_hash}{
The current hash object used for the statistics. See argument stat_hash and details. If argument hash = FALSE this is NULL.
}
\item{pvalue_hash}{
The current hash object used for the p-values. See argument stat_hash and details. If argument hash = FALSE this is NULL.
}
}

\references{
Peter Spirtes, Clark Glymour, and Richard Scheines. Causation, Prediction, and Search. The MIT Press, Cambridge, MA, USA, second edition, January 2001.

Lee Rodgers J., and Nicewander W.A. (1988). "Thirteen ways to look at the correlation coefficient." The American Statistician 42(1): 59-66.

Shevlyakov G. and Smirnov P. (2011). Robust Estimation of the Correlation Coefficient: An Attempt of Survey. Austrian Journal of Statistics, 40(1 & 2): 147-156.
}

\author{
Vincenzo Lagani and Ioannis Tsamardinos (Copyright 2012)

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> Vincenzo Lagani <vlagani@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{testIndSpearman}, \link{testIndReg}, \link{SES}, \link{testIndLogistic}, \link{gSquare}, \link{CondIndTests}}
}

\examples{

#simulate a dataset with continuous data
dataset <- matrix(runif(1000 * 200, 1, 1000), nrow = 1000 )
#the target feature is the last column of the dataset as a vector
target <- dataset[, 200]
results <- testIndFisher(target, dataset, xIndex = 44, csIndex = 100)

#require(gRbase)  #for faster computations in the internal functions
#define class variable (here tha last column of the dataset)
target <- 200;
#run the SES algorithm using the testIndFisher conditional independence test
sesObject <- SES(target, dataset, max_k = 3, threshold = 0.05, test = "testIndFisher");

#print summary of the SES output
summary(sesObject);
#plot the SES output
plot(sesObject, mode = "all");

}

\keyword{ Fisher's Test}
\keyword{ Conditional Independence Test }
