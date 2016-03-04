\name{nei}
\alias{nei}
\title{
Returns and plots, if asked, the node(s) and their neighbour(s), if there are any.
}

\description{
Returns and plots, if asked, the node(s) and their neighbour(s) of one or more nodes (if there are any).
}

\usage{
nei(G, ver, graph = TRUE)
}

\arguments{
  \item{G}{
  The adjacency matrix as produced from \code{\link{mmhc.skel}} or any other algorithm. 
}
  \item{ver}{
  A vector with one or more numbers indicating the seleted node(s) (or variable(s)). 
}
  \item{graph}{
  A boolean variable. If TRUE the relevane graph will appear (if there are neighbours). 
}
}

\details{
 This is a way to see the network for some given nodes. It is useful if you have many nodes and the whole network is a bit difficult to see clearly. 
}

\value{
A list object called "geit" containing the neighbours of the node(s). If there are no neighbours a message appears and no plot is presented. If the "graph" argument is set to TRUE and there are neighbours, a plot will appear. 

Bear in mind that the values can be extracted with the $ symbol, i.e. this is an S3 class output. 

}

%\references{
% 
%}

\author{
 Michail Tsagris

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> and Michail Tsagris <mtsagris@yahoo.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{mmhc.skel}, \link{SES}, \link{MMPC} }
}

\examples{
# simulate a dataset with continuous data
dataset <- matrix(runif(1000 * 50, 1, 100), nrow = 1000 ) 
a <- mmhc.skel(dataset, max_k = 3, threshold = 0.05, test = "testIndFisher", 
rob = FALSE, nc = 1) 
nei(a$G, 2)
nei(a$G, 2:4)
}

\keyword{ Network plot }
\keyword{ Neighbour nodes }