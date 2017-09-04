\name{Topological sort of a DAG}
\alias{topological_sort}
\title{
Topological sort of a DAG
}

\description{
Topological sort of a DAG.
}

\usage{
topological_sort(dag)
}

\arguments{
\item{dag}{
A square matrix representing a directed graph which contains either 0, 1 or 0, 2, and 3. In the first canse where G[i, j] = 1, means there is an arrow from node i to node j. 
In the second case G[i, j] = 2 and G[j, i] = 3 means that there is an arrow from node i to node j, where the 2 iindcates the arrohead and the 3 inducates the arrowtail. 
}
}

\details{
The function is an R translation from an old matlab code.
}

\value{
A vector with numbers indicating the sorting. 
}

\references{
Chickering, D.M. (1995). A transformational characterization of equivalent Bayesian network structures. Proceedings of the 11th Conference on Uncertainty in Artificial Intelligence, Montreal, Canada, 87-98. 
}

\author{
Ioannis Tsamardinos and Michail Tsagris

R implementation and documentation:  Ioannis Tsamardinos <tsamard@csd.uoc.gr> and Michail Tsagris <mtsagris@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{plotnetwork}, \link{nei}, \link{pc.or} 
}
}

\examples{
# simulate a dataset with continuous data
# simulate a dataset with continuous data
G = rdag(1000, 10, 0.3)$G
G[G == 2] <- 1
G[G == 3] <- 0
topological_sort(G)
}

\keyword{ Topological sort }
\keyword{ DAG }