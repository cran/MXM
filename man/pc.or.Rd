\name{Orientation rules for the PC algorithm}
\alias{pc.or}
\title{
The orientations part of the PC algorithm.
}

\description{
The function takes the outcome of the PC algorithm, as produced by \code{\link{pc.skel}} or \code{\link{pc.con}} and performes the 4 orientation rules. A graph is also possible to visualize.
}

\usage{
pc.or(mod, graph = FALSE) 

}

\arguments{
  \item{mod}{
  An object with the results of the PC algorithm, as produced by \code{\link{pc.skel}} or \code{\link{pc.con}}.   
}
\item{graph}{
Boolean that indicates whether or not to generate a plot with the graph. Package RgraphViz is required.
}
}

\details{
After having calculated the skeleton of the PC algorithm one may wants to perform orientations, leading to causal relationships. The rules as stated in Spirtes, Glymour and Scheines (2001) are
\enumerate{
\item \bold{Rule 1}. For each triple of vertices X, Y, Z such that the pair X, Y and the pair Y, Z are each adjacent in C but the pair X, Z are not adjacent in C, orient X - Y - Z as X -> Y <- Z if and only if Y is not in Sepset(X, Z). 
\item \bold{Rule 2}. If A -> B, B and C are adjacent, A and C are not adjacent, and there is no arrowhead at B, then orient B - C as B -> C.
\item \bold{Rule 3}. If there is a directed path from A to B, and an edge between A and B, then orient A - B as A -> B.
\item \bold{Rule 4}. This was added by Zhang (2008). If A -> B <- C, A - D - C, A and C are not adjacent, and D - B, then orient D - B as D -> B. 
}
The first rule is applied once. Rules 2-4 are applied repeatedly until no more edges can be oriented. 
}

\value{
A list including:
\item{Gini}{
The initial adjacency matrix, no orientations. This is the matrix produced by \code{\link{pc.skel}} or \code{\link{pc.con}}. 
}
\item{G}{
The final adjaceny matrix with the orientations. If G[i, j] = 2 then G[j, i] = 3. This means that there is an arrow from node i to node j. If G[i, j] = G[j, i] = 0; there is no edge between nodes i and j. If G[i, j] = G[j, i] = 1; there is an (undirected) edge between nodes i and j.
}
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
}

\references{
Spirtes P.,  Glymour C. and Scheines R. (2001). Causation, Prediction, and Search. The MIT Press, Cambridge, MA, USA, 3nd edition.

Zhang, Jiji. (2008). On the completeness of orientation rules for causal discovery in the presence of latent confounders and selection bias. Artificial Intelligence 172(16): 1873--1896.
}

\author{
 Michail Tsagris

R implementation and documentation: Giorgos Athineou <athineou@csd.uoc.gr> and Michail Tsagris <mtsagris@yahoo.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{pc.con}, \link{pc.skel}, \link{mmhc.skel}, \link{mb} 
}
}

\examples{
# simulate a dataset with continuous data
y = rdag(1000, 15, 0.3)
tru = y$G 
x = y$x
mod = pc.con(x)
mod$runtime

b = pc.or(mod)
plota(tru) 
dev.new()
plota(b$G)

dev.off()
plota( dag2eg(tru) )  ## essential graph
dev.new()
plota(b$G)

}

\keyword{ pc algorithm }