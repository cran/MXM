\name{Permutation based p-value for the Pearson correlation coefficient}
\alias{permcor}
\alias{permcorrels}
\title{
Permutation based p-value for the Pearson correlation coefficient
}

\description{
The main task of this test is to provide a p-value PVALUE for the null hypothesis: feature 'X' is independent from 'TARGET' given a conditioning set CS.
}

\usage{
permcor(x, R = 999) 
permcorrels(y, x, R = 999)
}

\arguments{
\item{x}{
For the case of "permcor" this is a matrix with two columns, two continuous variables. In the case of "permcorrels" this is a matrix with many variables. 
}
\item{y}{
A vector whose length is equal to the number of rows of x.
}
\item{R}{
The number of permutations to be conducted; set to 999 by default. 
}
}

\details{
This is a computational non parametric correlation coefficient test and is advised to be used when a small sample size is available. If you want to use the Spearman correlation instead, simply provide the ranks of x or of y and x. 
}

\value{
For the case of "permcor" a vector consisting of two values, the Pearson correlation and the permutation based p-value. 
For the "permcorrels" a vector with three values, the Pearson correlation, the test statistic value and the permutation based p-value.
}

\references{
Legendre Pierre (2000). Comparison of permutation methods for the partial correlation and partial Mantel tests. Journal of Statistical Computation and Simulation 67(1):37-73.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{pc.skel}, \link{testIndSpearman}, \link{testIndFisher}, \link{SES}, \link{CondIndTests}}
}

\examples{
permcor(iris[, 1:2])
permcor(iris[, 1:2], R = 9999)
x <- matrix(rnorm(50 * 2000), ncol = 2000)
a <- permcorrels(iris[1:50, 1], x)
}

\keyword{ Permutation test }
\keyword{ conditional independence test }
