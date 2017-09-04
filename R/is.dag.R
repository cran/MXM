is.dag <- function(dag) {
  if ( any( dag == 2 ) ) {
    dag[ dag == 1 ] <- 0
    dag[ dag == 3 ] <- 0
    dag[ dag == 2 ] <- 1
  } 
  !is.na( sum( topological_sort(dag) ) )
} 
       
     
  
  