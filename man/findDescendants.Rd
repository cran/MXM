\name{Descendants of a node in a directed graph}
\alias{findDescendants}
\title{
Returns and plots, if asked, the descendants of a node (or variable)
}

\description{
Returns and plots, if asked, the descendants of a node (or variable).
}

\usage{
findDescendants(G, node, graph = FALSE)
}

\arguments{
  \item{G}{
  The graph matrix as produced from \code{\link{pc.or}} or any other algorithm which produces directed graphs. 
}
  \item{node}{
  A numerical value indicating the node (or variable) whose descendants are to be returned. 
}
  \item{graph}{
  A boolean variable. If TRUE the relevant graph will appear (if there are descendants). 
}
}

\details{
The functions searches for the descendants of some node. This is an S3 class output. 
}

\value{
\item{isAnc}{
A matrix of the same dimensions as the original graph matrix with 0s and 1s. isAnc[i, j] = 1 indicates that the i-th node is an ancestor of the j-th node.
}
\item{Ganc}{
A matrix of dimensions equal to the number of descendants of the node with 0s and 1s. 
}
\item{anc}{
The descendants of the node.
}
}

%\references{
% 
%}

\author{
 Anna Roumpelaki

R implementation and documentation: Anna Roumpelaki <anna.roumpelaki@gmail.com>
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
plota(G)
findDescendants(G, 4, graph = FALSE)
}

\keyword{ Network plot }
\keyword{ Descendants }