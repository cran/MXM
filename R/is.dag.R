is.dag <- function(dag) {
  if ( any( dag == 2 ) ) {
    dag[ dag != 2 ] <- 0
    dag[ dag == 2 ] <- 1
  } 
  a <- topological_sort(dag)
  ( sum(a > 0, na.rm = TRUE) == dim(dag)[2] )
} 
       
     
  
  