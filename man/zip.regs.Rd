\name{Many simple zero inflated Poisson regressions}
\alias{zip.regs}
\title{
Many simple zero inflated Poisson regressions.
}

\description{
Many simple zero inflated Poisson regressions.
}

\usage{
zip.regs(target, dataset, wei = NULL, logged = FALSE, ncores = 1) 
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
Many simple zero inflated Poisson regressions are fitted. 
}

\value{
A matrix with the test statistic values, their relevant (logged) p-values and the BIC values. 
}

\references{
Lambert D. (1992). Zero-inflated Poisson regression, with an application to defects in manufacturing. Technometrics, 34(1):1-14.
}

\author{
 Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{zip.mod}, \link{testIndZIP}, \link{reg.fit}, \link{ridge.reg} }
}

\examples{
y <- rpois(1000, 3)
x <- matrix( rnorm(1000 * 20), ncol = 20)
y[1:100] <- 0
a <- zip.regs(y, x)
}

\keyword{ Regression models }
\keyword{ zero inflated Poisson }