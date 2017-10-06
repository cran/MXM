\name{Skeleton of the PC algorithm}
\alias{pc.skel}
\alias{pc.con}
\title{
The skeleton of a Bayesian network produced by the PC algorithm
}

\description{
The skeleton of a Bayesian network produced by the PC algorithm. No orientations are involved. The pc.con is for continuous data only and it calls the same C++ as pc.skel. Hence, you are advised to use pc.skel. 
}

\usage{
pc.skel(dataset, method = "pearson", alpha = 0.01, rob = FALSE, R = 1) 

pc.con(dataset, method = "pearson", alpha = 0.01) 
}

\arguments{
\item{dataset}{
A matrix with the variables. The user must know if they are continuous or if they are categorical. \bold{If you have categorical data though, the user must transform the data.frame into a matrix. In addition, the numerical matrix must have values starting from 0. For example, 0, 1, 2, instead of "A", "B" and "C"}. In the case of mixed variables, continuous, binary and ordinal this must a data.frame and the non continuous variables must be ordered factors, even the binary variables.
}
\item{method}{
If you have continuous data, you can choose either "pearson", "spearman" or "distcor". The latter uses the distance correlation and should not be used with lots of observations as it is by default really slow. If you have categorical data, this must be "cat". If you have a mix of continuous, binary and ordinal data (we will expand the available dataset in the future) then choose "comb.fast" or "comb.mm". These two methods perform the symmetric test for mixed data (Tsagris et al., 2017). See details for more information on this.
}
\item{alpha}{
The significance level ( suitable values in (0, 1) ) for assessing the p-values. Default value is 0.01.
}
\item{rob}{
A boolean variable which indicates whether (TRUE) or not (FALSE) to use a robust version of the statistical test if it is available. It takes more time than a 
non robust version but it is suggested in case of outliers. Default value is FALSE. This will only be taken into account if method is "pearson". 
}
\item{R}{
The number of permutations to be conducted. This is taken into consideration for the "pc.skel" only. The Pearson correlation coefficient is calculated and the 
p-value is assessed via permutations. There was an extra argument for plotting the skeleton but it does not work with the current visualisation 
packages, hence we removed the argument. Use \code{\link{plotnetwork}} to plot the skeleton.   
}
}

\details{
The PC algorithm as proposed by Spirtes et al. (2000) is implemented. The variables must be either continuous or categorical, only. The skeleton of the PC algorithm is order independent, since we are using the third heuristic (Spirte et al., 2000, pg. 90). At every ste of the alogirithm use the pairs which are least statistically associated. The conditioning set consists of variables which are most statistically associated with each either of the pair of variables. 

For example, for the pair (X, Y) there can be two coniditoning sets for example (Z1, Z2) and (W1, W2). All p-values and test statistics and degrees of freedom have been computed at the first step of the algorithm. Take the p-values between (Z1, Z2) and (X, Y) and between (Z1, Z2) and (X, Y). The conditioning set with the minimum p-value is used first. If the minimum p-values are the same, use the second lowest p-value. In the event of 2 or more p-values being the same (with permutations for example), the test statistic divided by the degrees of freedom is used as a means of choosing which conditioning set is to be used first. If two or more p-values are below the machine epsilon (.Machine$double.eps which is equal to 2.220446e-16), all of them are set to 0. To make the comparison or the ordering feasible we use the logarithm of the p-value. Hence, the logarithm of the p-values is always calculated and used.

In the case of the \eqn{G^2} test of independence (for categorical data) we have incorporated a rule of thumb. I the number of samples is at least 5 times the number of the parameters to be estimated, the test is performed, otherwise, independence is not rejected (see Tsamardinos et al., 2006).

The "comb.fast" and "comb.mm" methods are used with mixed variables, continuous, binary and ordinal. The "comb.mm" performs two log-likleihood ratio tests. For every pair of variables each of the two variables is treated as response and the suitable regression model is fitted. Then, two likelihood ratio tests are performed and the 2 p-values are combined in a meta-analytic way. In the case of "comb.fast" one regression model is fitted, the easiest (between the two) to implement. The ordering of the "easiness" is as follows: linear regression > logistic regression > ordinal regression.  

The pc.con is a faster implementation of the PC algorithm but for continuous data only, without the robust option, unlike pc.skel which is more general and even for the continuous datasets slower. pc.con accepts only "pearson" and "spearman" as correlations. 

If there are missing values they are placed by their median in case of continuous data and by their mode (most frequent value) if they are categorical. 
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
\item{density}{
The number of edges divided by the total possible number of edges, that is #edges / \eqn{n(n-1)/2}, where \eqn{n} is the number of variables.
}
\item{info}{
Some summary statistics about the edges, minimum, maximum, mean, median number of edges.
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

Tsagris M., Borboudakis G., Lagani V. and Tsamardinos I. (2017). Constraint-based Causal Discovery with Mixed Data. 
The 2017 ACM SIGKDD Workshop on Causal Discovery, 14/8/2017, Halifax, Nova Scotia, Canada. 
http://nugget.unisa.edu.au/CD2017/papersonly/constraint-based-causal-r1.pdf

Sedgewick, A. J., Ramsey, J. D., Spirtes, P., Glymour, C., & Benos, P. V. (2017). Mixed Graphical Models for Causal Analysis of Multi-modal Variables. arXiv preprint arXiv:1704.02621.

Szekely G.J. and Rizzo, M.L. (2014). Partial distance correlation with methods for dissimilarities. The Annals of Statistics, 42(6): 2382--2412.

Szekely G.J. and Rizzo M.L. (2013). Energy statistics: A class of statistics based on distances. Journal of Statistical Planning and Inference 143(8): 1249--1272.
}

\author{
 Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{MMPC}, \link{mmhc.skel}, \link{corfs.network}  }
}

\examples{
# simulate a dataset with continuous data
dataset <- matrix( runif(1000 * 50, 1, 100), nrow = 1000 ) 
a <- mmhc.skel(dataset, max_k = 3, threshold = 0.05, test = "testIndFisher" ) 
b <- pc.skel( dataset, method = "pearson", alpha = 0.05 ) 
b2 <- pc.con( dataset, method = "pearson" ) 
a$runtime ## 
b$runtime ## 
b2$runtime ## check the diffrerences in the runtimes
}

\keyword{ pc algorithm }