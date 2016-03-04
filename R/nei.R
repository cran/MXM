############################
############################
##### Return and plot if requested the neighbours of an undirected graph
#####
############################
############################

nei <- function(G, ver, graph = TRUE) {
  ## G is the adjacency matrix of an UN-DIRECTED graph
  ## ver is a number between 1 and the number of nodes
  ## it is a node whose neighbors you want to find
  ## it can also be avector with more than one nodes
  if ( is.null( colnames(G) ) ) {
    p <- ncol(G)
    colnames(G) <- paste("X", 1:p, sep = "")
  }

    geit <- list()
    for (i in 1:length(ver)) {
      geit[[ i ]] <- which( G[ver[i], ] == 2)
    }
    names(geit) <- colnames(G)[ver]
    ind <- unique( c( ver, as.vector( unlist(geit) ) ) )
    Gnei <- G[ind, ind]

    if ( class(Gnei) == "numeric" ) {
      graph = FALSE
      geit <- paste("There are no neighbours of the chosen node(s)");
    } 

    if (graph == TRUE) {
      if(requireNamespace("Rgraphviz", quietly = TRUE, warn.conflicts = FALSE) == TRUE) {
        am.graph <- new("graphAM", adjMat = Gnei, edgemode = "undirected")
        plot(am.graph, main = "Association network graph")
      } else {
        warning('In order to plot the generated network, package Rgraphviz is required.')
      }
    }
  geit
}