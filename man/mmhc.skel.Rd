\name{Skeleton of the max-min hill-climbing (MMHC) algorithm}
\alias{mmhc.skel}
\title{
The skeleton of a Bayesian network as produced by MMHC
}

\description{
The skeleton of a Bayesian network produced by MMHC. No orientations are involved.
}

\usage{
mmhc.skel(dataset, max_k = 3, threshold = 0.05, test = "testIndFisher", type = "MMPC", 
backward = TRUE, symmetry = TRUE, nc = 1, ini.pvalue = NULL, hash = FALSE) 
}

\arguments{
\item{dataset}{
A matrix with the variables. The user must know if they are continuous or if they are categorical. If you have a matrix with categorical data, i.e. 0, 1, 2, 3 where each number indicates a category, 
the minimum number for each variable must be 0. data.frame is also supported, as the dataset in this case is converted into a matrix.   
}
\item{max_k}{
The maximum conditioning set to use in the conditional indepedence test (see Details of SES or MMPC). 
}
\item{threshold}{
Threshold ( suitable values in (0, 1) ) for assessing p-values significance. Default value is 0.05.
}
\item{test}{
The conditional independence test to use. Default value is "testIndFisher". This procedure allows for "testIndFisher", "testIndSPearman" for continuous variables and "gSquare" for categorical variables. In case the dataset is a dta.frame with mixed types of data leave this NULL and an appropriate test will be selected. See \code{link{MMPC}} for the automatic choice of tests.
}
\item{type}{
The type of variable selection to take place for each variable (or node). The default (and standard) in "MMPC". You can also choose to run it via "SES" and thus allow for multiple signatures of variables to be connected 
to a variable.
}
\item{backward}{
If TRUE, the backward (or symmetry correction) phase will be implemented. This removes any falsely included variables in the parents and children set of the target variable. It call the \code{link{mmpcbackphase}} for this purpose. For perm.ses and wald.ses this is not yet applicable.
}
\item{symmetry}{
In order for an edge to be added, a statistical relationship must have been found from both directions. If you want this symmetry correction to take place, leave this boolean variable to TRUE. If you set it to FALSE, then if 
a relationship between Y and X is detected but not between X and Y, the edge is still added. 
}
\item{nc}{
How many cores to use. This plays an important role if you have many variables, say thousands or so. You can try with nc = 1 and with nc = 4 for example to see the differences. If you have a multicore machine, this is a must 
option. There was an extra argument for plotting the skeleton but it does not work with the current visualisation 
packages, hence we removed the argument. Use \code{\link{plotnetwork}} to plot the skeleton.   
}
\item{ini.pvalue}{
This is a list with the matrix of the univariate p-values. If you want to run mmhc.skel again, the univariate associations need not be calculated again.
}
\item{hash}{
A boolean variable which indicates whether (TRUE) or not (FALSE) to store the statistics calculated during MMPC or SES execution in a hash-type object. Default value is FALSE. If TRUE a hashObject is produced. You should put it to TRUE is you want an accurate count of the the number of tests performed. The only reason I put FLASE is because I want it to pass the CRAN check.
}
}

\details{
The MMPC is run on every variable. The backward phase (see Tsamardinos et al., 2006) takes place automatically. After all variables have been used, the matrix is checked for inconsistencies and they are corrected. 

A trick mentioned in that paper to make the procedure faster is the following. In the k-th variable, the algorithm checks how many previously scanned variables have an edge with the this variable and keeps them (it discards the other variables with no edge) along with the next (unscanned) variables. 

This trick reduces time, but can lead to different results. For example, if the i-th variable is removed, the k-th node might not remove an edge between the j-th variable, simply because the i-th variable that could d-sepate them is missing. 

The user is given this option via the argument "fast", which can be either TRUE or FALSE. Parallel computation is also available. 

}

\value{
A list including:
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
\item{density}{
The number of edges divided by the total possible number of edges, that is #edges / \eqn{n(n-1)/2}, where \eqn{n} is the number of variables.
}
\item{info}{
Some summary statistics about the edges, minimum, maximum, mean, median number of edges.
}
\item{ms}{
If you run "MMPC" for each variable this is NULL. If you run "SES" is a vector denoting which variables had more than one signature, i.e. more than one set of variables associated with them.
}
\item{ntests}{
The number of tests MMPC (or SES) performed at each variable.
}
\item{ini.pvalue}{
A matrix with the p-values of all pairwise univariate assocations.
}
\item{pvalue}{
A matrix with the final p-values. These are the maximum p-values calculated during the process. When the process finishes, the matrix is not symmetric. It becomes symmetric though by keeping the maximum between any two off-diagonal elements. These p-values now can be used to sort the strength of the edges. If you know the true adjacency matrix you can use them and create a ROC curve (see \code{link{bn.skel.utils}} for more information).
}
\item{G}{
The adjancency matrix. A value of 1 in G[i, j] appears in G[j, i] also, indicating that i and j have an edge between them.
}

Bear in mind that the values can be extracted with the $ symbol, i.e. this is an S3 class output. 
}

\references{
Tsamardinos, Brown and Aliferis (2006). The max-min hill-climbing Bayesian network structure learning algorithm. Machine learning, 65(1), 31-78.

Brown L. E., Tsamardinos I., and Aliferis C. F. (2004). A novel algorithm for scalable and accurate Bayesian network learning. Medinfo, 711-715.
}

\author{
 Michail Tsagris

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> and Michail Tsagris <mtsagris@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{pc.skel}, \link{pc.or}, \link{corfs.network}, \link{bn.skel.utils} }
}

\examples{
# simulate a dataset with continuous data
y <- rdag2(500, p = 20, nei = 3)
x <- y$x
a <- mmhc.skel(x, max_k = 5, threshold = 0.01, test = "testIndFisher" ) 
b <- pc.skel( x, alpha = 0.01 ) 
a$runtime ## 
b$runtime ## 
 
}

\keyword{ MMHC algorithm }