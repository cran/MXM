\name{Backward selection with generalised linear regression models}
\alias{glm.bsreg}
\title{
Variable selection in generalised linear regression models with backward selection

}

\description{
Variable selection in generalised linear regression models with backward selection
}

\usage{
glm.bsreg(target, dataset, threshold = 0.05, wei = NULL, heavy = FALSE, robust = FALSE)
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
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. It is not suggested when robust is set to TRUE.
}
\item{heavy}{
A boolean variable specifying whether heavy computations are required or not. If for exmaple the dataset contains tens of thousands of rows, it is advised to used memory efficient GLMs and hence set this to TRUE.
}
\item{robust}{
A boolean variable which indicates whether (TRUE) or not (FALSE) to use a robust version of the statistical test if it is available. It takes more time than a non robust version but it is suggested in case of outliers. 
Default value is FALSE. Currently only the linear regression option is supported.
}
}

\details{
This functions currently implements only linear, binary logistic and Poisson regression. If the sample size is less than the number of variables a meesage will appear and no backward regression is performed. 
}

\value{
The output of the algorithm is S3 object including:
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
\item{info}{
A matrix with the variables and their latest test statistics and p-values.
}
\item{ci_test}{
The conditional independence test used. 
}
\item{final}{
The final regression model.
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

#simulate a dataset with continuous data
dataset <- matrix( runif(1000 * 10, 1, 100), ncol = 10 )

#define a simulated class variable 
target <- rpois(1000, 10)
a <- glm.bsreg(target, dataset, threshold = 0.05) 

target <- rnorm(1000)
b <- glm.bsreg(target, dataset, threshold = 0.05)

target <- rbinom(1000, 1, 0.6)
d <- glm.bsreg(target, dataset, threshold = 0.05)

}

\keyword{ Markov Blanket }
\keyword{ Variable Selection }
