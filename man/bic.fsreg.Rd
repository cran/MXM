\name{BIC based forward selection}
\alias{bic.fsreg}
\title{
Variable selection in regression models with forward selection using BIC
}

\description{
Variable selection in regression models with forward selection using BIC
}

\usage{
bic.fsreg( target, dataset, test = NULL, robust = FALSE, tol = 0, ncores = 1 )
}

\arguments{
\item{target}{
The class variable. Provide either a string, an integer, a numeric value, a vector, a factor, an ordered factor or a Surv object. See also Details.
}
\item{dataset}{
The dataset; provide either a data frame or a matrix (columns = variables, rows = samples). The data can be either euclidean, categorical or both.
}
\item{test}{
The regression model to use. Available options are 
"gaussian" for normal linear regression, "beta" for beta regression, "Cox" for Cox proportional hazards, "Weibull" for Weibull regression, "binary" for binomial regression, "multinomial" for multinomial regression, "ordinal" for ordinal regression, "poisson" for poisson regression, "nb" for negative binomial regression, "zip" for zero inflated poisson regression and "speedglm" for linear, binary or poisson regression with large datasets (tens of thousands of observations).
}
\item{robust}{
A boolean variable which indicates whether (TRUE) or not (FALSE) to use a robust version of the statistical test if it is available. It takes more time than a non robust version but it is suggested in case of outliers. Default value is FALSE.
}
\item{tol}{
The difference bewtween two successive values of the stopping rule. By default this is is set to 2. If for example, the BIC difference between two succesive models is less than 2, the process stops and the last variable, even though significant does not enter the model.
}
\item{ncores}{
How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes and tens of thousands of variables and a regression based test which requires numerical optimisation. In other cammmb it will not make a difference in the overall time (in fact it can be slower). The parallel computation is used in the first step of the algorithm, where univariate associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in comparison to 1 core. Note also, that the amount of reduction is not linear in the number of cores.
}
}

\details{
If the current 'test' argument is defined as NULL or "auto" and the user_test argument is NULL then the algorithm automatically selects the best test based on the type of the data. Particularly:
\itemize{
	\item if target is a factor, the multinomial or the binary logistic regression is used. If the target has two values only, binary logistic regression will be used.
	\item if target is a ordered factor, the ordinal regression is used.
	\item if target is a numerical vector or a matrix with at least two columns (multivariate) linear regression is used. 
	\item if target is discrete numerical (counts), the poisson regression conditional independence test is used. If there are only two values, the binary logistic regression is to be used.
	\item if target is a Surv object, the Survival conditional independence test is used.
}
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

\references{
Tsamardinos I., Aliferis C. F. and Statnikov, A. (2003). Time and sample efficient discovery of Markov blankets and direct causal relations. In Proceedings of the 9th ACM SIGKDD international conference on Knowledge discovery and data mining (pp. 673-678). 
}


\author{
Michail Tsagris

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> Michail Tsagris <mtsagris@yahoo.gr>
}

\seealso{
\code{\link{glm.fsreg}, \link{lm.fsreg}, \link{bic.glm.fsreg}, \link{CondIndTests}, \link{MMPC}, \link{SES}}
}

\examples{
set.seed(123)
#require(gRbase) #for faster computations in the internal functions
require(hash)

#simulate a dataset with continuous data
dataset <- matrix( runif(1000 * 50, 1, 100), ncol = 50 )

#define a simulated class variable 
target <- 3 * dataset[, 10] + 2 * dataset[, 20] + 3 * dataset[, 30] + rnorm(1000, 0, 5)
a1 <- bic.fsreg(target, dataset, robust = FALSE, tol = 2, ncores = 1 ) 
a2 <- bic.fsreg(target, dataset, robust = TRUE, tol = 2, ncores = 1 ) 
a3 <- MMPC(target, dataset, robust= FALSE, ncores = 1)

target <- round(target)
b1 <- bic.fsreg(target, dataset, robust = FALSE, tol = 2, ncores = 1 ) 
b2 <- bic.fsreg(target, dataset, robust = TRUE, tol = 2, ncores = 1 ) 
#b3 <- MMPC(target, dataset, robust= FALSE, ncores = 1)
}

\keyword{ Markov Blanket }
\keyword{ Variable Selection }
