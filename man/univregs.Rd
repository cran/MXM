\name{Univariate regression based tests}
\alias{univregs}
\title{
Univariate regression based tests
}

\description{
Univariate regression based tests.
}

\usage{
univregs(target, dataset, test, wei = NULL, robust = FALSE, ncores = 1) 
}

\arguments{
\item{target}{
The target (dependent) variable. It must be a numerical vector with integers. 
}
\item{dataset}{
The indendent variable(s). It can be a vector, a matrix or a dataframe with continuous only variables, a data frame with mixed or only categorical variables.
}
\item{test}{
This argument should be one of the following strictly. \bold{testIndFisher}, \bold{testIndSpearman}, \bold{gSquare}, \bold{testIndBeta}, \bold{testIndReg}, \bold{testIndSpeedglm}, \bold{testIndLogistic}, \bold{testIndPois}, \bold{testIndZip}, \bold{censIndCR} or \bold{censIndWR}. Note that you must give the name of the test, without "". 
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. 
}
\item{robust}{
A boolean variable which indicates whether (TRUE) or not (FALSE) to use a robust version of the statistical test if it is available. It takes more time than non robust version but it is suggested in case of outliers. Default value is FALSE as it is currently nor supported.
}
\item{ncores}{
How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes and tens of thousands of variables and a regression based test which requires numerical optimisation. In other cases it will not make a difference in the overall time (in fact it can be slower). The parallel computation is used in the first step of the algorithm, where univariate associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in comparison to 1 core. Note also, that the amount of reduction is not linear in the number of cores.
}
}

\details{
This function is more as a help function for SES and MMPC, but it can also be called directly by the user. In some, one should specify the regression model to use and the function will perform all simple regressions, i.e. all regression models between the target and each of the variables in the dataset. 
}

\value{
A list including:
\item{stat}{
The value of the test statistic.
}
\item{pvalue}{
The logged p-value of the test. 
}
\item{flag}{
A number, either 1 (test performed) or 0 (test not performed).
}
}

%\references{

%}

\author{
 Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{SES}, \link{MMPC}, \link{CondIndTests}, \link{reg.fit}, \link{ridge.reg} }
}

\examples{
y <- rpois(100, 15)
x <- matrix( rnorm(100 * 20), ncol = 20)
a1 <- univregs(y, x, testIndPois)
}

\keyword{ Regression models }
