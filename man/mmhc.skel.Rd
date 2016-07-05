\name{Skeleton of the max-min hill-climbing (MMHC) algorithm}
\alias{mmhc.skel}
\title{
The skeleton of a Bayesian network as produced by MMHC
}

\description{
The skeleton of a Bayesian network produced by MMHC. No orientations are involved.
}

\usage{
mmhc.skel(dataset, max_k = 3, threshold = 0.05, test = NULL, rob = FALSE, fast = FALSE,
nc = 1, graph = FALSE) 
}

\arguments{
\item{dataset}{
A matrix with the variables. The user must know if they are continuous or if they are categorical. data.frame or matrix are both supported, as the dataset is converted into a matrix.   
}
\item{max_k}{
The maximum conditioning set to use in the conditional indepedence test (see Details of SES or MMPC). 
}
  \item{threshold}{
Threshold ( suitable values in (0, 1) ) for assessing p-values significance. Default value is 0.05.
}
  \item{test}{
The conditional independence test to use. Default value is "testIndFisher". This procedure allows for "testIndFisher", "testIndSPearman" for continuous variables and "gSquare" for categorical variables. 
}
\item{rob}{
A boolean variable which indicates whether (TRUE) or not (FALSE) to use a robust version of the statistical test if it is available. It takes more time than a non robust version but it is suggested in case of outliers. Default value is FALSE. This will only be taken into account if test is "testIndFisher".
}
\item{fast}{
A bollean variable indicating a faster procedure to take place. By default this is set to FALSE. See details about this. 
}
\item{nc}{
How many cores to use. This plays an important role if you have many variables, say thousands or so. You can try with nc = 1 and with nc = 4 for example to see the differences. If you have a multicore machine, this is a must option.
}
\item{graph}{
Boolean that indicates whether or not to generate a plot with the graph. Package RgraphViz is required.
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
\item{G}{
The adjancency matrix. A value of 1 in G[i, j] appears in G[j, i] also, indicating that i and j have an edge between them.
}

Bear in mind that the values can be extracted with the $ symbol, i.e. this is an S3 class output. 
}
\references{
Tsamardinos, Brown and Aliferis (2006). The max-min hill-climbing Bayesian network structure learning algorithm. Machine learning, 65(1), 31-78.
}

\author{
 Michail Tsagris

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> and Michail Tsagris <mtsagris@yahoo.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{SES}, \link{MMPC}, \link{pc.skel} }
}

\examples{
# simulate a dataset with continuous data
dataset <- matrix(runif(1000 * 50, 1, 100), nrow = 1000 ) 
a1 <- mmhc.skel(dataset, max_k = 3, threshold = 0.05, test = "testIndFisher", 
rob = FALSE, nc = 1) 
a2 <- mmhc.skel(dataset, max_k = 3, threshold = 0.05, test = "testIndSpearman", 
rob = FALSE, nc = 1)
a3 <- pc.con(dataset)
a4 <- pc.skel(dataset, R = 1) 
## check the diffrerences in the runtimes
a1$runtime  
a2$runtime 
a3$runtime 


}

\keyword{ MMHC algorithm }