\name{Beta regression}
\alias{beta.mod}
\alias{beta.reg}
\title{
Beta regression
}

\description{
Beta regression.
}

\usage{
beta.mod(target, dataset, wei = NULL, xnew= NULL) 
beta.reg(target, dataset, wei = NULL) 
}

\arguments{
\item{target}{
The target (dependent) variable. It must be a numerical vector with integers. 
}
\item{dataset}{
The indendent variable(s). It can be a vector, a matrix or a dataframe with continuous only variables, a data frame with mixed or only categorical variables.
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. 
}
\item{xnew}{
If you have new values for the predictor variables (dataset) whose target variable you want to predict insert them here. If you put the "dataset" or leave it NULL it will calculate the regression fitted values.
}
}

\details{
The beta regression is fitted. The "beta.reg" is an internal wrapper function and is used for speed up purposes. It is not to be called directly by the user unless they know what they are doing. 
}

\value{
A list including:
\item{iters}{
The iterations required until convergence.
}
\item{phi}{
The estimated precision parameter. 
}
\item{be}{
The estimated coefficients of the model.
}
\item{loglik}{
The log-likelihood of the regression model. 
}
}

\references{
Ferrari S.L.P. and Cribari-Neto F. (2004). Beta Regression for Modelling Rates and Proportions. Journal of Applied Statistics, 31(7): 799-815.
}

\author{
 Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{beta.regs}, \link{testIndBeta}, \link{reg.fit}, \link{ridge.reg} }
}

\examples{
y <- rbeta(500, 3, 5)
x <- matrix( rnorm(500 * 2), ncol = 2)
a1 <- beta.mod(y, x)
w <- runif(500)
a2 <- beta.mod(y, x, w)
}

\keyword{ Regression model }
\keyword{ beta distribution }