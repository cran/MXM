\name{ conditional independence test for longitudinal and clustered data}
\alias{testIndGLMM}
\title{
Linear mixed models conditional independence test for longitudinal class variables
}

\description{
The main task of this test is to provide a p-value PVALUE for the null hypothesis: feature 'X' is independent from 'TARGET' given a conditioning set CS. The pvalue is calculated by comparing a linear model based on the conditioning set CS against a model whose GLMMressor are both X and CS. The comparison is performed through an F test the appropriate degrees of freedom on the difference between the deviances of the two models. This test accepts a longitudinal target and longitudinal, categorical, continuous or mixed data as predictor variables. 
}

\usage{
testIndGLMM(target, reps = NULL, group, dataset, xIndex, csIndex, dataInfo=NULL, 
univariateModels=NULL, hash = FALSE, stat_hash=NULL, pvalue_hash=NULL, 
target_type=0, slopes=FALSE)
}

\arguments{
  \item{target}{
  A numeric vector containing the values of the target variable. If the values are proportions or percentages, i.e. strictly within 0 and 1 they are mapped into R using log( target/(1 - target) ). In both cases a linear mixed model is applied. It can also be a binary variable (binary logistic regression) or a discrete, counts (Poisson regression), thus fitting generalised linear mixed models.
  }
  \item{reps}{
  A numeric vector containing the time points of the subjects. It's length is equal to the length of the target variable. If you have clustered data, leave this NULL. 
}
  \item{group}{
  A numeric vector containing the subjects or groups. It must be of the same legnth as target. 
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
A boolean variable which indicates whether (TRUE) or not (FALSE) to use tha hash-based implementation of the statistics of SES. Default value is FALSE. If TRUE you have to specify the stat_hash argument and the pvalue_hash argument.
}
  \item{stat_hash}{
A hash object (hash package required) which contains the cached generated statistics of a SES run in the current dataset, using the current test.
}
  \item{pvalue_hash}{
A hash object (hash package required) which contains the cached generated p-values of a SES run in the current dataset, using the current test.
}
  \item{target_type}{
The type of the target variable. It is set to 0 by default. 1 for continuous variable (normal), 2 for binary variable (binomial) and 3 for discrete variable (Poisson).
}
  \item{slopes}{
A boolean variable which indicates whether (TRUE) to or not (FALSE) random slopes in the time effect as well. By deault random intercepts are considered. 
}
}

\details{
If hash = TRUE, testIndGLMM requires the arguments 'stat_hash' and 'pvalue_hash' for the hash-based implementation of the statistic test. These hash Objects are produced or updated by each run of SES (if hash == TRUE) and they can be reused in order to speed up next runs of the current statistic test. If "SESoutput" is the output of a SES.temp run, then these objects can be retrieved by SESoutput@hashObject$stat_hash and the SESoutput@hashObject$pvalue_hash.

The time effect is linear. It could be of higer order as well, but this would be a hyper-parameter, increasing the complexity of the models to be tested. 

Important: Use these arguments only with the same dataset that was used at initialization.

For all the available conditional independence tests that are currently included on the package, please see "?CondIndTests".
}

\value{
A list including:
\item{pvalue}{
A numeric value that represents the logarithm of the generated p-value due to the (generalised) linear mixed model (see reference below).
}
\item{stat}{
A numeric value that represents the generated statistic due to the (generalised) linear mixed model (see reference below).
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
Pinheiro, Jose, and Douglas Bates. Mixed-effects models in S and S-PLUS. Springer Science \& Business Media, 2006.
}

\author{
Vincenzo Lagani, Ioannis Tsamardinos, Michail Tsagris and Giorgos Athineou

R implementation and documentation: Giorgos Athineou <athineou@csd.gr>, Vincenzo Lagani <vlagani@csd.uoc.gr> and Michail Tsagris <mtsagris@yahoo.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{SES.temporal}, \link{MMPC.temporal}, \link{CondIndTests}}
}

\examples{
#data(sleepstudy)
#attach(sleepstudy)
#target <- Reaction
#x <- matrix(rnorm(180 * 10),ncol = 10) ## unrelated preidctor variables
#testIndGLMM(target, Days, Subject, x, 1,0,target_type = 1)
}

\keyword{ Linear mixed model }
\keyword{ Conditional independence test }
