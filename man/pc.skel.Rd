\name{Skeleton of the PC algorithm}
\alias{pc.skel}
\alias{pc.con}
\title{
The skeleton of a Bayesian network produced by the PC algorithm
}

\description{
The skeleton of a Bayesian network produced by the PC algorithm. No orientations are involved. The pc.con is a faster implementation for continuous datasets only. pc.skel is more general. 
}

\usage{
pc.skel(dataset, method = "pearson", alpha = 0.05, rob = FALSE, graph = FALSE) 

pc.con(dataset, method = "pearson", alpha = 0.05, graph = FALSE) 
}

\arguments{
  \item{dataset}{
  A matrix with the variables. The user must know if they are continuous or if they are categorical. data.frame or matrix are both supported, as the dataset is converted into a matrix.   
}
  \item{method}{
If you have continuous data, you can either choose "pearson" or "spearman". If you have categorical data though, this must be "cat".
}
  \item{alpha}{
The significance level ( suitable values in (0, 1) ) for assessing the p-values. Default value is 0.05.
}
\item{rob}{
A boolean variable which indicates whether (TRUE) or not (FALSE) to use a robust version of the statistical test if it is available. It takes more time than a non robust version but it is suggested in case of outliers. Default value is FALSE. This will only be taken into account if method is "pearson". 
}
\item{graph}{
Boolean that indicates whether or not to generate a plot with the graph. Package RgraphViz is required.
}
}

\details{
The PC algorithm as proposed by Spirtes et al. (2000) is implemented. The variables must be either continuous or categorical, only. The skeleton of the PC algorithm is order independent, since we are using the third heuristic (Spirte et al., 2000, pg. 90). At every ste of the alogirithm use the pairs which are least statistically associated. The conditioning set consists of variables which are most statistically associated with each either of the pair of variables. 

For example, for the pair (X, Y) there can be two coniditoning sets for example (Z1, Z2) and (W1, W2). All p-values and test statistics and degrees of freedom have been computed at the first step of the algorithm. Take the p-values between (Z1, Z2) and (X, Y) and between (Z1, Z2) and (X, Y). The conditioning set with the minimum p-value is used first. If the minimum p-values are the same, use the second lowest p-value. If the unlikely, but not impossible event of all p-values being the same, the test statistic divided by the degrees of freedom is used as a means of choosing which conditioning set is to be used first. 

If two or more p-values are below the machine epsilon (.Machine$double.eps which is equal to 2.220446e-16), all of them are set to 0. To make the comparison or the ordering feasible we use the logarithm of the p-value. Hence, the logarithm of the p-values is always calcualte dand used.

In the case of the \eqn{G^2} test of independence we have incorporated a rule of thumb. I the number of samples is at least 5 times the number of the parameters to be estimated, the test is performed, otherwise, independence is not rejected (see Tsmardinos et al., 2006, pg. 43).

The pc.con is a faster implementation of the PC algorithm but for continuous data only, without the robust option, unlike pc.skel which is more general and even for the continuous datasets slower. pc.con accepts only "pearson" and "spearman" as correlations. 
If in addition, you have more than 1000 variables, a trick to calculate the correlation matrix is implemented which can reduce the time required by this first step by nearly 50%. 

}

\value{
A list including:
\item{stat}{
The test statistics of the univariate associations. 
}
\item{pvalue}{
The logarithm of the p-values of the univariate associations. 
}
\item{info}{
Some summary statistics about the edges, minimum, maximum, mean, median number of edges.
}
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
\item{kappa}{
The maximum value of k, the maximum cardinality of the conditioning set at which the algorithm stopped.
}
\item{G}{
The adjancency matrix. A value of 1 in G[i, j] appears in G[j, i] also, indicating that i and j have an edge between them.
}
\item{sepset}{
A list with the separating sets for every value of k.
}
\item{title}{
The name of the dataset. 
}

Bear in mind that the values can be extracted with the $ symbol, i.e. this is an S3 class output. 

}
\references{
Spirtes P.,  Glymour C. and Scheines R. (2001). Causation, Prediction, and Search. The MIT Press, Cambridge, MA, USA, 3nd edition.
}

\author{
 Michail Tsagris

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> and Michail Tsagris <mtsagris@yahoo.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{SES}, \link{MMPC}, \link{mmhc.skel} }
}

\examples{
# simulate a dataset with continuous data
dataset <- matrix( runif(1000 * 50, 1, 100), nrow = 1000 ) 
a <- mmhc.skel(dataset, max_k = 3, threshold = 0.05, test = "testIndFisher" ) 
b <- pc.skel( dataset, method = "pearson", alpha = 0.05 ) 
b2 <- pc.con( dataset, method = "pearson" ) 
a$runtime ## 
b$runtime ## 
b2$rntime ## check the diffrerences in the runtimes
}

\keyword{ pc algorithm }