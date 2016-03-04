\name{plota}
\alias{plota}
\title{
Plot of an undirected graph
}

\description{
Plot of an undirected graph.
}

\usage{
plota(G)
}

\arguments{
  \item{G}{
  The adjacency matrix as produced from \code{\link{mmhc.skel}} or any other algorithm. Currently this only supports undirected graphs.
}
}

\details{
 This visualises an undirected graph.
}

\value{
The plot of the undirected graph.
}

\references{
\% 
}

\author{
 Michail Tsagris

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> and Michail Tsagris <mtsagris@yahoo.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{mmhc.skel}, \link{nei}, \link{MMPC} }
}

\examples{
# simulate a dataset with continuous data
dataset <- matrix(runif(1000 * 20, 1, 100), nrow = 1000 ) 
a <- mmhc.skel(dataset, max_k = 3, threshold = 0.05, test = "testIndFisher", 
rob = FALSE, nc = 1) 
plota(a$G)
}

\keyword{ Network plot }
\keyword{ Undirected graph }