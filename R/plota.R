#####################
##### Plot a graph using the adjacency matrix
#####
####################

plota <- function(G) {
  if(requireNamespace("Rgraphviz", quietly = TRUE, warn.conflicts = FALSE) == TRUE) {
    am.graph <- new("graphAM", adjMat = G, edgemode = "undirected")
    plot(am.graph, main = "Association network graph")
  } else {
    warning('In order to plot the generated network, package Rgraphviz is required.')
  }
}