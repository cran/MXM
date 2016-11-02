## finds descenadants of some given node(s)

findDescendants <- function(G, node = NULL, graph = FALSE) {

  n <- dim(G)[1]
  dag <- matrix(0, n, n)
  dag[ G == 3  &  t(G) == 2 ] <- 1
  isDesc <- transitiveClosure(dag)
  
  if ( is.null(node) )  {
  
    res <- isDesc
	
  } else {
    desc <- as.vector( isDesc[, node] )
    desc <- which( desc > 0 )
    Gdesc <- G[c(node, desc), c(node, desc)]
  
    if ( graph == TRUE) {
  
      if ( length(desc) > 0 ) {
	
        if ( requireNamespace("Rgraphviz", quietly = TRUE, warn.conflicts = FALSE) == TRUE ) {
          Gdesc[ Gdesc != 2 ] <- 0
          g <- as( Gdesc, "graphNEL" )
          plot(g, main = paste("Completed partially directed graph with ancestors of ", node ) )
        } else {
          warning('In order to plot the generated network, package Rgraphviz is required.')
        }
      } 
	
    }

    res <- list(isDesc = isDesc, Gdesc = Gdesc, desc = desc)     
	
  }
  
  res
  
}


