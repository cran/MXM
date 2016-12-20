\name{Many simple beta regressions}
\alias{beta.regs}
\title{
Many simple beta regressions.
}

\description{
Many simple beta regressions.
}

\usage{
beta.regs(target, dataset, wei = NULL, logged = FALSE, ncores = 1) 
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
\item{logged}{
A boolean variable; it will return the logarithm of the pvalue if set to TRUE.
}
\item{ncores}{
The number of cores to use. The default value is 1.
}
}

\details{
Many simple beta regressions are fitted. 
}

\value{
A matrix with the test statistic values, their relevant (logged) p-values and the BIC values. 
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
\code{ \link{beta.mod}, \link{testIndBeta}, \link{reg.fit}, \link{ridge.reg} }
}

\examples{
y <- rbeta(500, 5, 3)
x <- matrix( rnorm(500 * 20), ncol = 20)
a <- beta.regs(y, x)
}

\keyword{ Regression models }
\keyword{ beta distribution }