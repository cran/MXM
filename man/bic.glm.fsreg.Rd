\name{BIC based forward selection with generalised linear models}
\alias{bic.glm.fsreg}
\title{
Variable selection in generalised linear models with forward selection based on BIC

}

\description{
Variable selection in generalised linear models with forward selection based on BIC
}

\usage{
bic.glm.fsreg( target, dataset, test = NULL, robust = FALSE, tol = 0, ncores = 1 ) 
}

\arguments{
\item{target}{
The class variable. It can be either a vector with binary data (binomial regression), counts (poisson regression). If none of these is identified, linear regression is used.
}
\item{dataset}{
The dataset; provide either a data frame or a matrix (columns = variables, rows = samples). These can be continous and or categorical.
}
\item{test}{
This can be either "gaussian" for a continuous (univariate and multivariate) target, "binary" for a binary target, "poisson" for counts or NULL (default value).
}
\item{robust}{
A boolean variable which indicates whether (TRUE) or not (FALSE) to use a robust version of the statistical test if it is available. It takes more time than a non robust version but it is suggested in case of outliers. Default value is FALSE and is currently supported only by linear regression 
}
\item{tol}{
The difference bewtween two successive values of BIC. By default this is is set to 2. If for example, the BIC difference between two succesive models is less than 2, the process stops and the last variable, even though significant does not enter the model. 
}
\item{ncores}{
How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes and tens of thousands of variables and a regression based test which requires numerical optimisation. In other cammmb it will not make a difference in the overall time (in fact it can be slower). The parallel computation is used in the first step of the algorithm, where univariate associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in comparison to 1 core. Note also, that the amount of reduction is not linear in the number of cores.
}
%\item{maxit}{
%The maximum number of iterations required by the MM estimators in the robust binary or poisson regrression. This argument is taken into %account only if the argument robust is set to TRUE.
%}
}

\details{
Forward selection via the BIC is implemented. A variable which results in a reduction of BIC will be included, until the reduction is below a threshold set by the user (argument "tol").
}

\value{
The output of the algorithm is S3 object including:
\item{mat}{
A matrix with the variables and their latest test statistics and p-values.
}
\item{info}{
A matrix with the selected variables, their p-values and test statistics. Each row corresponds to a model which contains the variables up to that line. The BIC in the last column is the BIC of that model.
}
\item{models}{
The regression models, one at each step.
}
\item{final}{
The final regression model.
}
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
}

%\references{
%Tsamardinos I., Aliferis C. F. and Statnikov, A. (2003). Time and sample efficient discovery of Markov blankets and direct causal %relations. In Proceedings of the 9th ACM SIGKDD international conference on Knowledge discovery and data mining (pp. 673-678). 
%}


\author{
Michail Tsagris

R implementation and documentation: Giorgos Aathineou <athineou@csd.uoc.gr> Michail Tsagris <mtsagris@yahoo.gr>
}

\seealso{
\code{\link{fs.reg}, \link{lm.fsreg}, \link{bic.fsreg}, \link{CondIndTests}, \link{MMPC}, \link{SES}}
}

\examples{
set.seed(123)
#require(gRbase) #for faster computations in the internal functions
require(hash)

#simulate a dataset with continuous data
dataset <- matrix( runif(1000 * 50, 1, 100), ncol = 50 )

#define a simulated class variable 
target <- 3 * dataset[, 10] + 2 * dataset[, 20] + 3 * dataset[, 30] + rnorm(1000, 0, 5)
a <- bic.glm.fsreg(target, dataset, robust = FALSE, tol = 2, ncores = 1 ) 
}

\keyword{ Markov Blanket }
\keyword{ Variable Selection }
