\name{glm.bsreg}
\alias{glm.bsreg}
\title{
Variable selection in generalised linear regression models with backward selection

}

\description{
Variable selection in generalised linear regression models with backward selection
}

\usage{
glm.bsreg(target, dataset, threshold = 0.05)
}

\arguments{
\item{target}{
The class variable. Provide either an integer, a numeric value, or a factor. See also Details.
}
\item{dataset}{
The dataset; provide either a data frame or a matrix (columns = variables, rows = samples). In either case, only two cases are avaialble, either all data are continuous, or categorical. 
}
\item{threshold}{
Threshold (suitable values in [0,1]) for asmmmbsing p-values significance. Default value is 0.05.
}
}

\details{
This functions currently implements only linear, binary logistic and Poisson regression. If the sample size is less than the number of variables a meesage will appear and no backward regression is performed. 
}

\value{
The output of the algorithm is S3 object including:
\item{mat}{
A matrix with the variables and their latest test statistics and p-values.
}
\item{final}{
The final regression model.
}
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
}

%\references{
%Tsamardinos I., Aliferis C. F. and Statnikov, A. (2003). Time and sample efficient discovery of Markov blankets and direct causal relations. In Proceedings of the 9th ACM SIGKDD international conference on Knowledge discovery and data mining (pp. 673-678). 
%}


\author{
Michail Tsagris

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> Michail Tsagris <mtsagris@yahoo.gr>
}

\seealso{
\code{\link{fs.reg}, \link{lm.fsreg}, \link{bic.fsreg}, \link{bic.glm.fsreg}, \link{CondIndTests}, \link{MMPC}, \link{SES}}
}

\examples{
set.seed(123)
#require(gRbase) #for faster computations in the internal functions
require(hash)

#simulate a dataset with continuous data
dataset <- matrix( runif(1000 * 50, 1, 100), ncol = 50 )

#define a simulated class variable 
target <- rpois(1000, 10)

a <- glm.bsreg(target, dataset, threshold = 0.05) 
}

\keyword{ Markov Blanket }
\keyword{ Variable Selection }
