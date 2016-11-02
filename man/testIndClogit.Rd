\name{Conditional independence test for case control data}
\alias{testIndClogit}
\title{
Conditional independence test based on conditional logistic regression for case control studies 
}

\description{
The main task of this test is to provide a p-value PVALUE for the null hypothesis: feature 'X' is independent from 'TARGET' given a conditioning set CS. The pvalue is calculated by comparing a conditional logistic regression model based on the conditioning set CS against a model whose regressor are both X and CS. The comparison is performed through a chi-square test with the appropriate degrees of freedom on the difference between the deviances of the two models. This is suitable for a case control design
}

\usage{
testIndClogit(target, dataset, xIndex, csIndex, wei = NULL, dataInfo = NULL, 
univariateModels = NULL, hash = FALSE, stat_hash = NULL, pvalue_hash = NULL, 
robust = FALSE)
}

\arguments{
\item{target}{
A matrix with two columns, the first one must be 0 and 1, standing for 0 = control and 1 = case. The second column is the id of the patients. A numerical variable, for example c(1,2,3,4,5,6,7,1,2,3,4,5,6,7). 
}
\item{dataset}{
A numeric matrix or a data.frame in case of categorical predictors (factors), containing the variables for performing the test. Rows as samples and columns as features.
}
\item{xIndex}{
The index of the variable whose association with the target we want to test.
}
\item{csIndex}{
The indices of the variables to condition on.
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL and it should stay NULL as weights are ignored the conditional logistic regression. See the \bold{survival} package for more information about this type of regression.
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
A boolean variable which indicates whether (TRUE) or not (FALSE) to use a robustified version of Beta regression. Currently it is not available for this test.
}
}

\details{
If hash = TRUE, testIndClogit requires the arguments 'stat_hash' and 'pvalue_hash' for the hash-based implementation of the statistic test. These hash Objects are produced or updated by each run of SES (if hash == TRUE) and they can be reused in order to speed up next runs of the current statistic test. If "SESoutput" is the output of a SES run, then these objects can be retrieved by SESoutput@hashObject$stat_hash and the SESoutput@hashObject$pvalue_hash.

Important: Use these arguments only with the same dataset that was used at initialization.

For all the available conditional independence tests that are currently included on the package, please see "?CondIndTests".
}

\value{
A list including:
\item{pvalue}{
A numeric value that represents the logarithm of the generated p-value due to Beta regression (see reference below).
}
\item{stat}{
A numeric value that represents the generated statistic due to Beta regression (see reference below).
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
Ferrari S.L.P. and Cribari-Neto F. (2004). Beta Regression for Modelling Rates and Proportions. Journal of Applied Statistics, 31(7): 799-815.
}

\author{
Vincenzo Lagani, Ioannis Tsamardinos, Michail Tsagris and Giorgos Athineou

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr>, Vincenzo Lagani <vlagani@csd.uoc.gr> and Michail Tsagris <mtsagris@yahoo.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{SES}, \link{testIndLogistic}, \link{censIndCR}, \link{censIndWR}}
}

\examples{
#simulate a dataset with continuous data
dataset <- matrix(rnorm(300 * 100), nrow = 300 ) 
#the target feature is the last column of the dataset as a vector
case <- rbinom(300, 1, 0.6)
ina <- which(case==1)
ina <- sample(ina, 100)
case[-ina] = 0 
id <- rep(1:100,3)
target <- cbind(case, id)

results <- testIndClogit(target, dataset, xIndex = 44, csIndex = 60)
results

#run the SES algorithm using the testIndClogit conditional independence test
a1 <- SES(target, dataset, max_k = 3, threshold = 0.05, test = "testIndClogit");
a2 <- MMPC(target, dataset, max_k = 3, threshold = 0.05, test = "testIndClogit");
# print summary of the SES output
summary(a1);
# plot the SES output
# plot(a1, mode = "all");
}

\keyword{ Conditional logistic regression}
\keyword{ Case-control studies }
\keyword{ Conditional independence test }
