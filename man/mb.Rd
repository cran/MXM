\name{mb}
\alias{mb}
\title{
Returns and plots, if asked, the Markov blanket of a node (or variable).
}

\description{
Returns and plots, if asked, the Markov blanket of a node (or variable).
}

\usage{
mb(G, node, graph = FALSE)
}

\arguments{
  \item{G}{
  The graph matrix as produced from \code{\link{pc.or}} or any other algorithm which produces directed graphs. 
}
  \item{node}{
  A vector with one or more numbers indicating the seleted node(s) (or variable(s)). 
}
  \item{graph}{
  A boolean variable. If TRUE the relevant graph will appear (if the Markov blanket is the non emty set). 
}
}

\details{
This is a way to see the network for some given nodes. It is useful if you have many nodes and the whole network is a bit difficult to see clearly. 
Bear in mind that the values can be extracted with the $ symbol, i.e. this is an S3 class output. 
}

\value{
\item{parents}{
The parents of the node of interest.
}
\item{children}{
The children of the node of interest.
}
\item{spouses}{
The spouses of the node of interest. These are the other parents of the children of the node of interest. 
}
\item{relatives}{
Nodes which are connected with the node of interest, but it is not known whether they are parents or children. The edge between them is undirected.
}
\item{markov.blanket}{
The Markov blanket of the node of interest. The collection of all the previous.
}
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
\code{ \link{plota}, \link{nei}, \link{pc.or} 
}
}

\examples{
# simulate a dataset with continuous data
# simulate a dataset with continuous data
y = rdag(1000, 10, 0.3)
tru = y$G 
x = y$x
mod = pc.con(x)
G = pc.or(mod)$G
mb(G, 1, graph = FALSE)
}

\keyword{ Network plot }
\keyword{ Neighbour nodes }