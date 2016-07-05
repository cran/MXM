\name{Correlation based tests with and without permutation p-value}
\alias{condi}
\title{
Fisher conditional independence test for continuous class variables with and without permutation based p-value
}

\description{
The main task of this test is to provide a permutation based p-value PVALUE for the null hypothesis: feature 'X' is independent from 'TARGET' given a conditioning set CS.
}

\usage{
condi(ind1, ind2, cs, dat, type = "pearson", rob = FALSE, R = 1) 
}

\arguments{
\item{ind1}{
The index of the one variable to be considered. 
}
\item{ind2}{
The index of the other variable to be considered. 
}
\item{cs}{
The index or indices of the conditioning set of variable(s).
}
\item{dat}{
A matrix with the data. 
}
\item{type}{
Do you want the Pearson (type = "pearson") or the Spearman (type = "spearman") correlation to be used. 
}
\item{rob}{
If you choose type="pearson" then you can sapecify whether you want a robust version of it (TRUE) or not (FALSE).
}
\item{R}{
If R = 1 then the asymptotic p-value is calculated. If R > 1 a permutation based p-value is returned. 
}
}

\details{
This test is currently designed for usage by the PC algorithm. 
}

\value{
A vector including the test statistic, it's associated p-value and the relevant degrees of freedom. In the case of a permutation based p-value, the returned test statistic is the observed test statistic divided by the relevant degrees of freedom. This is for the case of ties between many permutation based p-values. The PC algorithm choose a pair of variables based on the p-values. If they are equal it will use the test statistic. 
}

\references{
Hampel F. R., Ronchetti E. M., Rousseeuw P. J., and Stahel W. A. (1986). Robust statistics: the approach based on influence functions. John Wiley & Sons. 

Peter Spirtes, Clark Glymour, and Richard Scheines. Causation, Prediction, and Search. The MIT Press, Cambridge, MA, USA, second edition, January 2001.

Lee Rodgers J., and Nicewander W.A. (1988). "Thirteen ways to look at the correlation coefficient." The American Statistician 42(1): 59-66.

Shevlyakov G. and Smirnov P. (2011). Robust Estimation of the Correlation Coefficient: An Attempt of Survey. Austrian Journal of Statistics, 40(1 & 2): 147-156.
}

\author{
Michail Tsagris

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> and Michail Tsagris <mtsagris@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{testIndFisher}, \link{testIndSpearman}, \link{pc.skel}, \link{gSquare}, \link{CondIndTests}
}
}

\examples{
#simulate a dataset with continuous data
dataset <- matrix(runif(1000 * 10, 1, 1000), nrow = 10 )
testIndFisher(dataset[, 1], dataset[, -1], xIndex = 1, csIndex = 2)
condi(ind1 = 1, ind2 = 2, cs = 3, dataset, R = 1)
condi(ind1 = 1, ind2 = 2, cs = 3, dataset, R = 999)
}

\keyword{ Fisher's Test}
\keyword{ Conditional Independence Test }
\keyword{ Permutation based p-value }
