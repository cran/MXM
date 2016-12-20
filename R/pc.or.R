pc.or <- function(mod, graph = FALSE) {
  ## mod is the outcome of the function pc.con or pc.skel
  ## method is either "pearson", "spearman" or "cat"
  ## alpha is the threshold used in the skeleton of the pc
  ## rob is either TRUE or FALSE
  ## sim is either TRUE or FALSE
  ## graph is either TRUE or FALSE
  G <- mod$G  ## the undirected graph
  n <- ncol(G)  ## how many variables are there
  ina <- 1:n 
  kappa <- mod$kappa
  sep <- mod$sepset  ## the triples with separating variables
  durat <-  proc.time()
  if ( length(sep) == 0 ) {
    G <- G 
    paste("There are no separating sets")
  } else {
    
    len <- length(sep) 
    if ( len == 1 )  {
      if ( !is.matrix( sep[[ 1 ]] ) ) {
        sepa <- c( sep[[ 1 ]][ 1:3 ], sep[[ 1 ]][ 5 ] )
      } else {
        sepa <- cbind( sep[[ 1 ]][, 1:3 ], sep[[ 1 ]][, 5] )
      }
    } else {
      if ( !is.matrix( sep[[ 1 ]] ) ) {
        sepa <- c( sep[[ 1 ]][ 1:3 ], numeric( len - 1 ), sep[[ 1 ]][ 5 ] )
      } else {
        ante <- matrix( numeric( nrow( sep[[ 1 ]] ) * (len - 1) ), ncol = len - 1)
        sepa <- cbind( sep[[ 1 ]][, 1:3 ], ante, sep[[ 1 ]][, 5] )
      }
      for (i in 2:len) {
        if ( !is.matrix( sep[[ i ]] ) ) {
          sep[[ i ]] <- c( sep[[ i ]][ 1:c(2 + i)], numeric( len - i ), sep[[ i ]][ 4 + i ] )
          sepa <- rbind( sepa, sep[[ i ]] ) 
        } else {
          ante <- matrix( numeric( nrow( sep[[ i ]] ) * (len - i) ), ncol = len - i )
          if ( sum( dim(ante) ) != 0 ) {
            la <- cbind( sep[[ i ]][, 1:c(2 + i) ], ante, sep[[ i ]][ , 4 + i ] )
          } else {
            la <- cbind( sep[[ i ]][, 1:c(2 + i) ], sep[[ i ]][ , 4 + i ] )
          }
          sepa <- rbind(sepa, la)
        }
      }
    }
    
    if ( !is.matrix(sepa) )  sepa <- matrix( sepa, nrow = 1 )
    colnames(sepa) <- c("X", "Y", paste("Var", 1:len, sep = ""), "logged p-value")
    m <- nrow(sepa)
    p <- ncol(sepa) - 1
    rownames(sepa) <- colnames(G)[ sepa[, 1] ]
    ## sepa contains all the non significant pairs given at least one variable
    
    vim <- 1  ## 1st step	
    ### Orientation rule 0: find the triplets (v structures) and orient
  	G <- R0(G, ina, sepa) 
	
    ## 1st orientation rule: if there is an arrow and then an edge make it arrow
    G <- R1(G)
	
    ## 2nd orientation rule: if there is a directed path connect the beginning with the end
    G <- R2(G) 
    
    ### 3rd rule, the complicated one with v structures
    G <- R3(G)
    
    G1 <- G
    
    ## 1st orientation rule: if there is an arrow and then an edge make it arrow
    G <- R1(G) 
    
    ## 2nd orientation rule: if there is a directed path connect the beginning with the end
    G <- R2(G)    
    
    ### 3rd rule, the complicated one with v structures
    G <- R3(G)
    
    G2 <- G
    
    vim <- 2 
    while ( sum( abs( G1 - G2 ) ) != 0 ) {
      
	  vim <- vim + 1 
      G1 <- G2   
      ## 1st orientation rule: if there is an arrow and then an edge make it arrow
      G <- R1(G)
      
      ## 2nd orientation rule: if there is a directed path connect the beginning with the end
      G <- R2(G) 
       
      ### 3rd rule, the complicated one with v structures
      G <- R3(G)
	        
      G2 <- G
    }
    
    G <- G2
    
  } ## end of all rules 
  
  
  durat <- proc.time() - durat
  colnames(G) <- rownames(G) <- colnames(mod$G)
  
  Ga <- G 
  Ga[ Ga == 3 ] <- 0
  
  if ( graph ) {
    
    plotnetwork(mod$G, titlos = paste("Skeleton of the PC algorithm for", mod$title ) )
    dev.new()
    plotnetwork(Ga, titlos = paste("CPDAG of the PC algorithm for", mod$title ))
  }
  
  final = list(Gini = mod$G, G = G, runtime = durat) 
  final
}









#############
### Rules  
#############

R0 <- function(G, ina, sepa) {
  
    p <- ncol(sepa) - 1
    l <- Rfast::rowsums( G == 1 )
    id <- which( l >= 0 ) 
    
    if ( length(id) > 0 ) {
      
      for (i in id) {
        adj <- ina[ G[i, ] == 1 ]
        if ( length(adj) > 1 ) {
          sam <-  as.matrix(  t( combn(adj, 2) ) )
          ela <- NULL
          for ( j in 1:nrow(sam) ) {
            if ( G[sam[j, 1], sam[j, 2] ] == 0 ) {
              res <- is.sepset( sam[j, ], i, sepa[, 1:p] )
              if ( !res ) {
                G[ sam[j, 1], i ] = G[ sam[j, 2], i ] = 2 
                G[ i, sam[j, 1] ] = G[ i, sam[j, 2] ] = 3
                
              } else G <- G  
			  
            }
          }
        }
      }
	  
    }
	
	G
	
}
    
	
	
R1 <- function(G) {

   if ( sum( G == 2 ) > 0  &  sum( G == 1 ) > 0 ) {
      tup <- which( G == 2, arr.ind = TRUE )
      nup <- nrow(tup) 
      for (i in 1:nup) {
        can <- tup[i, 2] 
        geit <- which( G[can, ] == 1 )
        # G[can, geit] <- 2
        # G[geit, can] <- 3
        if ( length(geit) > 0 ) {
          G[can, geit] <- 2
          G[geit, can] <- 3
        }  
      }
   } 
	
  G
}   


R2 <- function(G) {

   if ( sum( G == 2 ) > 0  &  sum( G == 1 ) > 0 ) {
      
      poia <- which( G == 1, arr.ind = TRUE ) 
      poia1 <- which( G == 1 )
      poia2 <- which( G == 3 )
      GGG <- G
      GGG[poia1] <- NA
      GGG[poia2] <- NA
      GGG[ GGG == 0 ] <- NA
      g1 <- e1071::allShortestPaths(GGG)
      nu <- nrow(poia)
      
      if ( nu > 0 ) {
        for ( i in 1:nrow(poia) )  {
          aa <- e1071::extractPath( g1, start = poia[i, 1], end = poia[i, 2] )
          if ( length(aa) > 2 ) {
            if ( G[ poia[i, 1], poia[i, 2] ] == 1 ) {
              G[ poia[i, 1], poia[i, 2] ] <- 2
              G[ poia[i, 2], poia[i, 1] ] <- 3
            }
          }
        }   
      }
      
   } 
	
   G
}	


R3 <- function(G) {
    if ( sum( G == 2 ) > 0 & sum( G == 1 ) > 0 ) {
        
       tup <- which( G == 2, arr.ind = TRUE )
       nup <- nrow(tup) 
       poia <- as.factor( tup[, 2] )
       poia <- as.numeric( levels(poia) ) 
       ela <- NULL

       for ( i in poia ) {
         b3 <- which(G[i, ] == 1) 
         if ( length( b3 ) > 0 ) {
           wa <- tup[ tup[, 2] == i, ]
           if ( is.matrix(wa) )  ela <- combn( 1:nrow(wa), 2 )
           if ( !is.null(ela) ) {
             for ( j in 1:ncol(ela) ) {
               wa2 <- as.vector( wa[ ela[, j], 1 ] ) 
               b1 <- which( G[ wa2[1], ] == 2 )
               b2 <- which( G[ wa2[2], ] == 2 ) 
               tomi <- intersect( b1, intersect(b3, b2) )
               if ( length(tomi) > 0 )  {
                 G[tomi, i] <- 2  
                 G[i, tomi] <- 3    
               }
             }
           } 
         }

         ela <- NULL
	   }  
    }
	  
	G 
}


is.sepset <- function(pair, nd, separ) {
  
  a <- which( Rfast::colsums(t(separ[, 1:2]) - pair ) == 0 )
  if ( length(a) > 1 )  { 
    b <- length( intersect(nd, separ[a, -c(1:2)]) ) > 0
  } else b <- FALSE
  
  b
} 