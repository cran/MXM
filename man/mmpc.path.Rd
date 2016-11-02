\name{MMPC solution paths for many combinations of hyper-parameters}
\alias{mmpc.path}
\title{
MMPC solution paths for many combinations of hyper-parameters
}

\description{
MMPC solution paths for many combinations of hyper-parameters.
}

\usage{
mmpc.path(target, dataset, wei = NULL, max_ks = NULL, thresholds = NULL, test = NULL,
user_test = NULL, robust = FALSE, ncores = 1)
}

\arguments{
  \item{target}{
The class variable. Provide either a string, an integer, a numeric value, a vector, a factor, an ordered factor or a Surv object. See also Details.
}
\item{dataset}{
The data-set; provide either a data frame or a matrix (columns = variables , rows = samples).
Alternatively, provide an ExpressionSet (in which case rows are samples and columns are features, see bioconductor for details).
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. 
}
\item{max_ks}{
A vector of possible max_k values. Can be a number as well, but this does not really make sense to do. If nothing is given, the values max_k=3 and max_k=2 are used by default. 
}
\item{thresholds}{
A vector of possible threshold values. Can be a number as well, but this does not really make sense to do. If nothing is given, the values (0.1, 0.05, 0.01) are used by default.  
}
\item{test}{
The conditional independence test to use. Default value is NULL. See also \code{\link{CondIndTests}}.
}
\item{user_test}{
A user-defined conditional independence test (provide a closure type object). Default value is NULL. If this is defined, the "test" argument is ignored.
}
\item{robust}{
A boolean variable which indicates whether (TRUE) or not (FALSE) to use a robust version of the statistical test if it is available. It takes more time than a non robust version but it is suggested in case of outliers. Default value is FALSE.
}
\item{ncores}{
How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes and tens of thousands of variables and a regression based test which requires numerical optimisation. In other cases it will not make a difference in the overall time (in fact it can be slower). The parallel computation is used in the first step of the algorithm, where univariate associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in comparison to 1 core. Note also, that the amount of reduction is not linear in the number of cores. This argument is used only in the first run of MMPC and for the univariate associations only and the results are stored (hashed). In the enxt runs of MMPC the results are used (cashed) and so the process is faster. 
}
}

\details{
For different combinations of the hyper-parameters, max_k and the significance level (threshold or alpha) the MMPC algorith is run. 
}

\value{
The output of the algorithm is an object of the class 'SESoutput' for SES or 'MMPCoutput' for MMPC including:
\item{bic}{
A matrix with the BIC values of the final fitted model based on the selected variables identified by each configuration, combination of the hyper-parameters.
}
\item{size}{
A matrix with the legnth of the selected variables identified by each configuration of MMPC.
}
\item{variables}{
A list containing the variables from each configuration of MMPC
}
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
}

\references{
I. Tsamardinos, V. Lagani and D. Pappas (2012). Discovering multiple, equivalent biomarker signatures. In proceedings of the 7th conference of the Hellenic Society for Computational Biology & Bioinformatics - HSCBB12.

Tsamardinos, Brown and Aliferis (2006). The max-min hill-climbing Bayesian network structure learning algorithm. Machine learning, 65(1), 31-78.
}

\author{
Ioannis Tsamardinos, Vincenzo Lagani

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> Vincenzo Lagani <vlagani@csd.uoc.gr>  
}

\seealso{
\code{\link{CondIndTests}, \link{cv.ses}}
}

\examples{
set.seed(123)
require(hash)
# simulate a dataset with continuous data
dataset <- matrix(runif(1000 * 101, 1, 100), nrow = 1000 ) 
#the target feature is the last column of the dataset as a vector
target <- dataset[, 101]
dataset <- dataset[, -101]

a <- mmpc.path(target, dataset, max_ks = NULL, thresholds = NULL, test = NULL, 
user_test = NULL, robust = FALSE, ncores = 1)
}

\keyword{ SES }
\keyword{ Multiple Feature Signatures }
\keyword{ Feature Selection }
\keyword{ Variable Selection }
