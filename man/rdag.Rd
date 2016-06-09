\name{rdag}
\alias{rdag}
\title{
G square conditional independence test for discrete data
}

\description{
The main task of this test is to provide a p-value PVALUE for the null hypothesis: feature 'X' is independent 
from 'TARGET' given a conditioning set CS. This test is based on the log likelihood ratio test.
}

\usage{
rdag(n, p, s, a = 0, m, A = NULL, seed = FALSE) 
}

\arguments{
  \item{n}{
  A number indicating the sample size. 
}
  \item{p}{
  A number indicating the number of nodes (or vectices, or variables).
}
  \item{s}{
  A number in \eqn{(0, 1)}. This defines somehow the sparseness of the model. It is the probability that a node has an edge. 
}
  \item{a}{
  A number in \eqn{(0, 1)}. The defines the percentage of outliers to be included in the simulated data. If \eqn{a=0}, no outliers are generated. 
}
  \item{m}{
  A vector equal to the number of nodes. This is the mean vector of the normal distribution from which the data are to be generated. This is used only when \eqn{a>0} 
 so as to define the mena vector of the multivariate normal from which the outliers will be generated. 
}
  \item{A}{
If you already have an an adjacency matrix in mind, plug it in here, otherwise, leave it NULL. 
}
  \item{seed}{
If seed is TRUE, the simulated data will always be the same. 
}
}

\details{
In the case where no adjacency matrix is given, an \eqn{p \times p} matrix with zeros everywhere is created. 
very element below the diagonal is is replaced by random values from a Bernoulli distribution with probability of success equal to s. 
This is the matrix B. Every value of 1 is replaced by a uniform value in \eqn{0.1, 1}. This final matrix is called A. 
The data are generated from a multivariate normal distribution with a zero mean vector and covariance matrix equal to 
\eqn{\left({\bf I}_p- A\right)^{-1}\left({\bf I}_p- A\right)}, where \eqn{{\bf I}_p} is the \eqn{p \times p} identiy matrix. 
If a is greater than zero, the outliers are generated from a multivariate normal with the same covariance matrix and mean vector the one 
specified by the user, the argument "m". The flexibility of the outliers is that you cna specifiy outliers in some variables only or in all of them. 
For example, m=c(0,0,5) introduces outliers in the third variable only, whereas m=c(5,5,5) introduces outliers in all variables. 
The user is free to decide on the type of outliers to include in the data.
}


\value{
A list including:
\item{nout}{
The number of outliers.
}
\item{G}{
The adcacency matrix used. If G[i, j] = 1, then G[j, i] = 2 and this means that there is an arrow from j to i.
}
\item{A}{
The matrix with the with the uniform values in the interval \eqn{0.1, 1}. 
}
\item{x}{
The simulated data. 
}
}

\references{
Colombo, Diego, and Marloes H. Maathuis (2014). Order-independent constraint-based causal structure learning." The Journal of Machine Learning Research 15(1): 3741--3782.
}

\author{
R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%}

\seealso{
\code{\link{SES}, \link{testIndFisher}, \link{testIndLogistic}, \link{censIndCR}, \link{CondIndTests}}
}

\examples{
y = rdag(100, 20, 0.2)
x = y$x
tru = y$G 

mod = pc.con(x)
b = pc.or(mod)
plota(tru) 
dev.new()
plota(b$G)

}

\keyword{ \eqn{G^2} test }
\keyword{ Log likelihood ratio }
\keyword{ Conditional independence test }
