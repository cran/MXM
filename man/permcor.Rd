\name{Permutation based p-value for the Pearson correlation coefficient}
\alias{permcor}
\title{
Permutation based p-value for the Pearson correlation coefficient
}

\description{
The main task of this test is to provide a p-value PVALUE for the null hypothesis: feature 'X' is independent from 'TARGET' given a conditioning set CS.
}

\usage{
permcor(x, R = 999) 
}

\arguments{
  \item{x}{
  A matrix with two columns, two continuous variables.
}
  \item{R}{
The number of permutations to be conducted; set to 999 by default. 
}
}

\details{
This is a computational non parametric correlation coefficient test and is advised to be used when a small sample size is available. This is internally called by \code{\link{condi.perm}}.
}

\value{
A vector consisting of two values, the Pearson correlation and the permutation based p-value. 
}

%\references{
%
%}

\author{
Vincenzo Lagani and Ioannis Tsamardinos

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> Vincenzo Lagani <vlagani@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{pc.skel}, \link{testIndSpearman}, \link{testIndFisher}, \link{SES}, \link{CondIndTests}}
}

\examples{
permcor(iris[,1:2])
permcor(iris[,1:2], R = 9999)
}

\keyword{ Permutation test }
\keyword{ Conditional Independence Test }
