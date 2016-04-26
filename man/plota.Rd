\name{plota}
\alias{plota}
\title{
Plot of an (un)directed graph
}

\description{
Plot of an (un)directed graph.
}

\usage{
plota(G)
}

\arguments{
  \item{G}{
  The adjacency matrix as produced from \code{\link{mmhc.skel}} or any other algorithm. This can correspond to an undirected, partially directed or a completely directed graph. 
}
}

\details{
 This visualises the directed graph.
}

\value{
The plot of the directed graph.
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
\keyword{ Directed graph }