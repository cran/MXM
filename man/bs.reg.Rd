\name{Backward selection regression}
\alias{bs.reg}
\title{
  Variable selection in regression models with forward selection
}

\description{
  Variable selection in regression models with forward selection
}

\usage{
bs.reg(target, dataset, threshold = 0.05, test = NULL, wei = NULL, user_test = NULL, 
heavy = FALSE, robust = FALSE)
}

\arguments{
\item{target}{
The class variable. Provide either a string, an integer, a numeric value, a vector, a factor, an ordered factor or a Surv object. See also Details.
}
\item{dataset}{
The dataset; provide either a data frame or a matrix (columns = variables, rows = samples). In either case, only two cases are avaialble, either all data are continuous, or categorical. 
  }
\item{threshold}{
Threshold (suitable values in [0,1]) for asmmmbsing p-values significance. Default value is 0.05.
}
\item{test}{
The regression model to use. Available options are most of the tests for SES and MMPC. The ones NOT available are "gSquare", "censIndER", "testIndMVreg", "testIndClogit", "testIndBinom", "testIndSpearman" and "testIndFisher".
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. It is not suggested when robust is set to TRUE.
}
\item{user_test}{
A user-defined conditional independence test (provide a closure type object). Default value is NULL. If this is defined, the "test" argument is ignored.
}
\item{heavy}{
A boolean variable specifying whether heavy computations are required or not. If for exmaple the dataset contains tens of thousands of rows, it is advised to used memory efficient GLMs and hence set this to TRUE.
}
\item{robust}{
A boolean variable which indicates whether (TRUE) or not (FALSE) to use a robust version of the statistical test if it is available.  It takes more time than a non robust version but it is suggested in case of outliers. Default value is FALSE.
}
  %\item{maxit}{
    %The maximum number of iterations required by the MM estimators in the robust binary or poisson regrression. This argument is %taken into account only if the argument robust is set to TRUE.
    %}
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

\author{
  Michail Tsagris
  
  R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> Michail Tsagris <mtsagris@yahoo.gr>
}

\seealso{
  \code{\link{glm.fsreg}, \link{lm.fsreg}, \link{bic.fsreg}, \link{bic.glm.fsreg}, \link{CondIndTests}, \link{MMPC}, \link{SES}}
}

\examples{
set.seed(123)
#simulate a dataset with continuous data
dataset <- matrix( runif(1000 * 10, 1, 100), ncol = 10 )

#define a simulated class variable 
target <- rnorm(1000)

a <- bs.reg(target, dataset, threshold = 0.05, test = "testIndRQ") 
b <- bs.reg(target, dataset, threshold = 0.05, test = "testIndReg") 
}

\keyword{ Backward regression }
\keyword{ Markov Blanket }
\keyword{ Variable Selection }
