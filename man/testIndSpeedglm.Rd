\name{Conditional independence test for continuous, binary and count data with thousands of samples}
\alias{testIndSpeedglm}
\title{
Conditional independence test for continuous, binary and discrete (counts) variables with thousands of observations
}

\description{
The main task of this test is to provide a p-value PVALUE for the null hypothesis: feature 'X' is independent from 'TARGET' given a conditioning set CS. The pvalue is calculated by comparing a logistic model based on the conditioning set CS against a model whose regressor are both X and CS. The comparison is performed through a chi-square test with the aproprirate degrees of freedom on the difference between the deviances of the two models.
}

\usage{
testIndSpeedglm(target,dataset,xIndex,csIndex,dataInfo = NULL,univariateModels = NULL, 
hash = FALSE, stat_hash = NULL, pvalue_hash = NULL, target_type = 0, robust = FALSE)
}

\arguments{
\item{target}{
A numeric vector containing the values of the target variable. It can be either continuous or percentages (values within 0 and 1), binary or discrete (counts). 
}
\item{dataset}{
A numeric matrix or data frame, in case of categorical predictors (factors), containing the variables for performing the test. Rows as samples and columns as features.
}
\item{xIndex}{
The index of the variable whose association with the target we want to test.
}
\item{csIndex}{
The indices of the variables to condition on.
}
\item{dataInfo}{
A list object with information on the structure of the data. Default value is NULL.
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
\item{target_type}{
A numeric vector that represents the type of the target. Default value is 0. See details for more.
\itemize{
\item target_type = 1 (binary target)
\item target_type = 2 (nominal target)
\item target_type = 3 (discrete target)
}
}
\item{robust}{
A boolean variable which indicates whether (TRUE) or not (FALSE) to use a robustified version of the logistic regressions available here. Currently it is not available for these cases.
}
}

\details{
If argument target_type=0 then testIndSpeedglm requires the dataInfo argument to indicate the type of the current target:
\itemize{
\item dataInfo$target_type = "normal" (continuous target)
\item dataInfo$target_type = "binary" (binary target)
\item dataInfo$target_type = "discrete" (discrete target)
}

If hash = TRUE, testIndSpeedglm requires the arguments 'stat_hash' and 'pvalue_hash' for the hash-based implementation of the statistic test. These hash Objects are produced or updated by each run of SES (if hash == TRUE) and they can be reused in order to speed up next runs of the current statistic test. If "SESoutput" is the output of a SES run, then these objects can be retrieved by SESoutput@hashObject$stat_hash and the SESoutput@hashObject$pvalue_hash.

Important: Use these arguments only with the same dataset that was used at initialization.

This test is designed for large sample sized data, tens and hundreds of thousands and it works for linear, logistic and poisson regression. The classical \code{\link{lm}} and \code{\link{glm}} functions will use too much memory when many observations are available. The package "speedglm" handles such data more efficiently. You can try and see, in the first case the computer will jam, whereas in the second it will not. Hence, this test is to be used in these cases only. We have not set a threshold on the sample size, so that the algorithm decides whether to shift to speedglm or not, because this depends upon the user's computing fascilities. When there are up to \eqn{20,000} observations, the built-in function \code{\link{lm}} is faster, but when \eqn{n=30,000}, the \code{\link{speedlm}} is more than twice as fast. 

For all the available conditional independence tests that are currently included on the package, please see "?CondIndTests".
}

\value{
A list including:
\item{pvalue}{
A numeric value that represents the logarithm of the generated p-value.
}
\item{stat}{
A numeric value that represents the generated statistic.
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
McCullagh, Peter, and John A. Nelder. Generalized linear models.  CRC press, USA, 2nd edition, 1989.
}

\author{
Vincenzo Lagani and Ioannis Tsamardinos 

R implementation and documentation: Vincenzo Lagani <vlagani@csd.uoc.gr>, Giorgos Athineou <athineou@csd.uoc.gr> and Michail Tsagris <mtsagris@yahoo.gr>
}

%\note{
%
% }

\seealso{
\code{\link{SES}, \link{testIndLogistic}, \link{testIndReg}, \link{testIndPois}, \link{CondIndTests}}
}

\examples{
#require(gRbase) #for faster computations in the internal functions
dataset <- matrix(runif(40000 * 10, 1, 50), ncol = 10 ) 
#the target feature is the last column of the dataset as a vector
target <- rpois(40000, 10)
system.time( testIndPois(target, dataset, xIndex = 1, csIndex = 2) )
system.time( testIndSpeedglm(target, dataset, xIndex = 1, csIndex = 2) )

}

\keyword{ Linear regression }
\keyword{ Binary logistic regression }
\keyword{ Poisson regression }
\keyword{ Conditional independence test }
