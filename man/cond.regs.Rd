\name{Conditional independence regression based tests}
\alias{cond.regs}
\title{
Conditional independence regression based tests
}

\description{
Conditional independence regression based tests.
}

\usage{
cond.regs(target, dataset, xIndex, csIndex, test = NULL, wei = NULL, ncores = 1)
}

\arguments{
\item{target}{
The target (dependent) variable. It must be a numerical vector. 
}
\item{dataset}{
The indendent variable(s). For the "univregs" this can be a matrix or a dataframe with continuous only variables, a data frame with mixed or only categorical variables. For the "wald.univregs", "perm.univregs", "score.univregs" and "rint.regs" this can only by a numerical matrix.
}
\item{xIndex}{
The indices of the variables whose association with the target you want to test.
}
\item{csIndex}{
The index or indices of the variable(s) to condition on. If this is 0, the the function \code{\link{univregs}} will be called.
}
\item{test}{
One of the following: \bold{testIndBeta}, \bold{testIndReg}, \bold{testIndLogistic}, \bold{testIndMultinom}, 
\bold{testIndOrdinal}, \bold{testIndPois}, \bold{testIndZIP}, \bold{testIndNB}, \bold{testIndClogit}, 
\bold{testIndBinom}, \bold{testIndIGreg}, \bold{censIndCR}, \bold{censIndWR}, \bold{censIndWR}, 
\bold{testIndGamma} or \bold{testIndNormLog}. 

Note that in all cases you must give the name of the test, without " ". 
}
\item{wei}{
A vector of weights to be used for weighted regression. The default value is NULL. 
}
\item{ncores}{
How many cores to use. This plays an important role if you have tens of thousands of variables or really large sample sizes and tens of thousands of variables and a regression based test which requires numerical optimisation. In other cases it will not make a difference in the overall time (in fact it can be slower). The parallel computation is used in the first step of the algorithm, where univariate associations are examined, those take place in parallel. We have seen a reduction in time of 50\% with 4 cores in comparison to 1 core. Note also, that the amount of reduction is not linear in the number of cores.
}
}

\details{
This function is more as a help function for SES and MMPC, but it can also be called directly by the user. In some, one should specify the regression model to use and the function will perform all simple regressions, i.e. all regression models between the target and each of the variables in the dataset. }

\value{
A list including:
\item{stat}{
The value of the test statistic.
}
\item{pvalue}{
The logarithm of the p-value of the test. 
}
}

\references{
Chen J. and Chen Z. (2008). Extended Bayesian information criteria for model selection with large model spaces. Biometrika, 95(3): 759-771.
}

\author{
 Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{univregs}, \link{SES}, \link{MMPC}, \link{CondIndTests}, \link{reg.fit}, \link{ridge.reg} }
}

\examples{
y <- rpois(50, 15)
x <- matrix( rnorm(50 * 10), ncol = 10)
a1 <- univregs(y, x, test = testIndPois)
a2 <- perm.univregs(y, x, test = permPois)
a3 <- wald.univregs(y, x, test = waldPois)
}

\keyword{ Regression models }
