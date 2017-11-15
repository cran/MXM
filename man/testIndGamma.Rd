\name{Conditional independence tests for positive data}
\alias{testIndGamma}
\alias{testIndNormLog}
\alias{permGamma}
\alias{permNormLog}
\alias{waldGamma}
\alias{waldNormLog}
\title{
Regression conditional independence test for positive response variables.
}

\description{
The main task of this test is to provide a p-value PVALUE for the null hypothesis: feature 'X' is independent from 'TARGET' given a conditioning set CS. The pvalue is calculated by comparing a Poisson regression model based on the conditioning set CS against a model whose regressor are both X and CS. The comparison is performed through a chi-square test with the appropriate degrees of freedom on the difference between the deviances of the two models. The models supported here are poisson, zero inlftaed poisson and negative binomial.
}

\usage{
testIndGamma(target, dataset, xIndex, csIndex, wei = NULL, dataInfo = NULL, 
univariateModels = NULL, hash = FALSE, stat_hash = NULL, pvalue_hash = NULL, 
robust = FALSE)
  
testIndNormLog(target, dataset, xIndex, csIndex, wei = NULL, dataInfo = NULL, 
univariateModels = NULL, hash = FALSE, stat_hash = NULL, pvalue_hash = NULL, 
robust = FALSE)

permGamma(target, dataset, xIndex, csIndex, wei = NULL, dataInfo = NULL, 
univariateModels = NULL, hash = FALSE, stat_hash = NULL, pvalue_hash = NULL, 
robust = FALSE, threshold = 0.05, R = 999)

permNormLog(target, dataset, xIndex, csIndex, wei = NULL, dataInfo = NULL, 
univariateModels = NULL, hash = FALSE, stat_hash = NULL, pvalue_hash = NULL, 
robust = FALSE, threshold = 0.05, R = 999)

waldGamma(target, dataset, xIndex, csIndex, wei = NULL, dataInfo = NULL, 
univariateModels = NULL, hash = FALSE, stat_hash = NULL, pvalue_hash = NULL, 
robust = FALSE) 

waldNormLog(target, dataset, xIndex, csIndex, wei = NULL, dataInfo = NULL, 
univariateModels = NULL, hash = FALSE, stat_hash = NULL, pvalue_hash = NULL,
robust = FALSE) 
}

\arguments{
\item{target}{
A numeric vector containing the values of the target variable. For the Gamma based tests, the values must be strictly greater than zero. For the NormLog case, zeros can be included.
}
\item{dataset}{
A numeric matrix or data frame, in case of categorical predictors (factors), containing the variables for performing the test. Rows as samples and columns as features. In the cases of "waldPois", "waldNB" and "waldZIP" this is strictly a matrix. 
}
\item{xIndex}{
The index of the variable whose association with the target we want to test.
}
\item{csIndex}{
The indices of the variables to condition on. If you have no variables set this equal to 0.
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. 
}
\item{dataInfo}{
A list object with information on the structure of the data. Default value is NULL.
}
\item{univariateModels}{
Fast alternative to the hash object for univariate test. List with vectors "pvalues" (p-values), "stats" (statistics) and "flags" (flag = TRUE if the test was succesful) representing the univariate association of each variable with the target. Default value is NULL.
}
\item{hash}{
A boolean variable which indicates whether (TRUE) or not (FALSE) to use tha hash-based implementation of the statistics of SES. Default value is FALSE. If TRUE you have to specify the stat_hash argument and the pvalue_hash argument.
}
\item{stat_hash}{
A hash object (hash package required) which contains the cached generated statistics of a SES run in the current dataset, using the current test.
}
\item{pvalue_hash}{
A hash object (hash package required) which contains the cached generated p-values of a SES run in the current dataset, using the current test.
}
\item{robust}{
A boolean variable which indicates whether (TRUE) or not (FALSE) to use a robust version of the statistical test if it is available. It takes more time than non robust version but it is suggested in case of outliers. Default value is FALSE as it is currently nor supported.
}
\item{threshold}{
Threshold (suitable values in [0,1]) for assessing p-values significance. 
}
\item{R}{
The number of permutations, set to 999 by default. There is a trick to avoind doing all permutations. As soon as the number of times the permuted test statistic is more than the observed test statistic is more than 50 (if threshold = 0.05 and R = 999), the p-value has exceeded the signifiance level (threshold value) and hence the predictor variable is not significant. There is no need to continue do the extra permutations, as a decision has already been made. 
}
}

\details{
If hash = TRUE, all three tests require the arguments 'stat_hash' and 'pvalue_hash' for the hash-based implementation of the statistic test. These hash Objects are produced or updated by each run of SES (if hash = TRUE) and they can be reused in order to speed up next runs of the current statistic test. If "SESoutput" is the output of a SES run, then these objects can be retrieved by SESoutput@hashObject$stat_hash and the SESoutput@hashObject$pvalue_hash.
  
Important: Use these arguments only with the same dataset that was used at initialization. For all the available conditional independence tests that are currently included on the package, please see "?CondIndTests".
  
For the testIndGamma and testIndNormLog the F test is used and not the log-likelihood ratio test because both of these regression models have a nusiance parameter. 
}

\value{
  A list including:
    \item{pvalue}{
      A numeric value that represents the logarithm of the generated p-value due to the count data regression (see references below).
    }
  \item{stat}{
    A numeric value that represents the generated statistic due to Poisson regression(see reference below).
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
McCullagh P., and Nelder J.A. (1989). Generalized linear models.  CRC press, USA, 2nd edition.
}

\author{
Michail Tsagris
  
R implementation and documentation: Michail Tsagris <mtsagris@csd.uoc.gr>
}

%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{\link{testIndSpeedglm}, \link{testIndNB}, \link{testIndZIP}, \link{gSquare}, \link{CondIndTests}}
}

\examples{
#simulate a dataset with continuous data
dataset <- matrix(rnorm(400 * 50, 1, 5), ncol = 50 ) 
#the target feature is the last column of the dataset as a vector
target <- rgamma(400, 1, 3)
testIndGamma(target, dataset, xIndex = 24, csIndex = 10)
testIndNormLog(target, dataset, xIndex = 24, csIndex = 10)
#run the SES algorithm using the testIndPois conditional independence test
sesObject <- SES(target, dataset, max_k = 3, threshold = 0.05, test = "testIndGamma");
sesObject2 <- SES(target, dataset, max_k = 3, threshold = 0.05, test = "testIndNormLog");
}

\keyword{ Gamma regression }
\keyword{ Log link }

