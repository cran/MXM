pc.or <- function(mod, prev.cycles = TRUE) {
  ## mod is the outcome of the function pc.con or pc.skel
  G <- mod$G  ## the undirected graph
  n <- ncol(G)  ## how many variables are there
  ina <- 1:n 
  sep <- mod$sepset  ## the triples with separating variables
  durat <-  proc.time()
  if ( length(sep) == 0 ) {
    G2 <- G 
    paste("There are no separating sets")
  } else {
 
    len <- length(sep) 
    if ( len == 1 )  {
      if ( !is.matrix( sep[[ 1 ]] ) ) {
        sepa <- c( sep[[ 1 ]][ 1:3 ], sep[[ 1 ]][ 5 ] )
      } else if ( dim( sep[[ 1 ]] )[1] == 1 ) {
        sepa <- c( sep[[ 1 ]][ 1:3 ], sep[[ 1 ]][ 5 ] )
      } else {
        sepa <- cbind( sep[[ 1 ]][, 1:3 ], sep[[ 1 ]][, 5] )
      }
    } else {
      if ( !is.matrix( sep[[ 1 ]] ) ) {
        sepa <- c( sep[[ 1 ]][ 1:3 ], numeric( len - 1 ), sep[[ 1 ]][ 5 ] )
      } else {
        ante <- matrix( numeric( nrow( sep[[ 1 ]] ) * (len - 1) ), ncol = len - 1)
        sepa <- cbind( sep[[ 1 ]][, 1:3, drop = FALSE ], ante, sep[[ 1 ]][, 5, drop = FALSE] )
      }
      for (i in 2:len) {
        if ( !is.matrix( sep[[ i ]] ) ) {
          sep[[ i ]] <- c( sep[[ i ]][ 1:c(2 + i)], numeric( len - i ), sep[[ i ]][ 4 + i ] )
          sepa <- rbind( sepa, sep[[ i ]] ) 
        } else {
          ante <- matrix(numeric( NROW( sep[[ i ]] ) * (len - i) ), ncol = len - i )
          if ( sum( dim(ante) ) > 2 ) {
            la <- cbind( sep[[ i ]][, 1:c(2 + i), drop = FALSE ], ante, sep[[ i ]][ , 4 + i, drop = FALSE ] )
          } else {
            la <- cbind( sep[[ i ]][, 1:c(2 + i), drop = FALSE ], numeric( len - i ), sep[[ i ]][ , 4 + i, drop = FALSE ] )
          }
          sepa <- rbind(sepa, la)
        }
      }
    }
    
    if ( !is.matrix(sepa) )  sepa <- matrix( sepa, nrow = 1 )
    colnames(sepa) <- c("X", "Y", paste("Var", 1:len, sep = ""), "logged p-value")
    rownames(sepa) <- colnames(G)[ sepa[, 1] ]
    ## sepa contains all the non significant pairs given at least one variable
    
    ### Orientation rule 0: find the triplets (v structures) and orient
    p <- dim(sepa)[2] - 1
    G1 <- R0(G, ina, sepa[, 1:p, drop = FALSE], prev.cycles = prev.cycles) 
    ## 1st orientation rule: if there is an arrow and then an edge make it arrow
    G1 <- R1(G1, prev.cycles = prev.cycles)
    ## 2nd orientation rule: if there is a directed path connect the beginning with the end
    G1 <- R2(G1, prev.cycles = prev.cycles) 
    ### 3rd rule, the complicated one with v structures
    G1 <- R3(G1, prev.cycles = prev.cycles)
    ## 1st orientation rule: if there is an arrow and then an edge make it arrow
    G2 <- R1(G1, prev.cycles = prev.cycles) 
    ## 2nd orientation rule: if there is a directed path connect the beginning with the end
    G2 <- R2(G2, prev.cycles = prev.cycles)    
    ### 3rd rule, the complicated one with v structures
    G2 <- R3(G2, prev.cycles = prev.cycles)
    
    while ( sum( abs( G1 - G2 ) ) != 0 ) {
      G1 <- G2   
      ## 1st orientation rule: if there is an arrow and then an edge make it arrow
      G2 <- R1(G1, prev.cycles = prev.cycles)
      ## 2nd orientation rule: if there is a directed path connect the beginning with the end
      G2 <- R2(G2, prev.cycles = prev.cycles) 
      ### 3rd rule, the complicated one with v structures
      G2 <- R3(G2, prev.cycles = prev.cycles)
    }
    
  } ## end of all rules 
    
  durat <- proc.time() - durat
  colnames(G2) <- rownames(G2) <- colnames(mod$G)
  
  final = list(Gini = mod$G, G = G2, runtime = durat) 
  final
}


#############
### Rules  
#############

R0 <- function(G, ina, sepa, prev.cycles = TRUE) {
  
    l <- Rfast::rowsums( G == 1 )
    id <- which( l >= 0 ) 

    if ( length(id) > 0 ) {
      
      for (i in id) {
        adj <- ina[ G[i, ] == 1 ]
        if ( length(adj) > 1 ) {
          sam <-  as.matrix(  t( combn(adj, 2) ) )
          for ( j in 1:nrow(sam) ) {
            if ( G[sam[j, 1], sam[j, 2] ] == 0  &  G[sam[j, 1], i ] == 1  &  G[i, sam[j, 2] ] == 1 ) {
              res <- is.sepset( sam[j, ], i, sepa )
              if ( !res ) {
                G[ sam[j, 1], i ] = 2
                G[ i, sam[j, 1] ] = 3
                G[ sam[j, 2], i ] = 2 
                G[ i, sam[j, 2] ] = 3
                if (prev.cycles) {
                  if ( !is.dag(G) )  {
                    G[ sam[j, 1], i ] = 1
                    G[ i, sam[j, 1] ] = 1
                    G[ sam[j, 2], i ] = 1 
                    G[ i, sam[j, 2] ] = 1
                  }
                }  
			        }	 
            } else if ( G[sam[j, 1], sam[j, 2] ] == 0  &  G[sam[j, 1], i ] == 2  &  G[i, sam[j, 2] ] == 1 ) {
              res <- is.sepset( sam[j, ], i, sepa )
              if ( !res ) {
                G[ sam[j, 2], i ] = 2 
                G[ i, sam[j, 2] ] = 3
                if (prev.cycles) {
                  if ( !is.dag(G) )  {
                    G[ sam[j, 2], i ] = 1 
                    G[ i, sam[j, 2] ] = 1
                  }
                }  
			        }
            } else if ( G[sam[j, 1], sam[j, 2] ] == 0  &  G[sam[j, 1], i ] == 1  &  G[i, sam[j, 2] ] == 2 ) {
              res <- is.sepset( sam[j, ], i, sepa )
              if ( !res ) {
                G[ sam[j, 1], i ] = 2 
                G[ i, sam[j, 1] ] = 3
                if (prev.cycles) {
                  if ( !is.dag(G) )  {
                    G[ sam[j, 1], i ] = 1 
                    G[ i, sam[j, 1] ] = 1
                  }
                }  
			        }
            } 
          }
        }
      } 
    }
  G	
}
    
	
	
R1 <- function(G, prev.cycles = TRUE) {

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
        if (prev.cycles) {
          if ( !is.dag(G) )  {
            G[ can, geit] = 1
            G[ geit, can ] = 1
          }
        }  
      }
   } 
	
  G
}   

R2 <- function(G, prev.cycles = TRUE) {

   if ( sum( G == 2 ) > 0  &  sum( G == 1 ) > 0 ) {
      tup <- which( G == 2, arr.ind = TRUE )
      nup <- nrow(tup) 
      for (i in 1:nup) {
        arxi <- tup[i, 1] 
        can <- tup[i, 2] 
        geit <- which( G[can, ] == 2 )
        if ( sum( G[arxi, geit] == 1 ) > 0 ) {
          G[arxi, geit] <- 2
          G[geit, arxi] <- 3
        }
        if (prev.cycles) {
          if ( !is.dag(G) )  {
            G[ arxi, geit ] = 1
            G[ geit, arxi ] = 1
          }
        }  
      }
   } 
	
  G
} 


R3 <- function(G, prev.cycles = TRUE) {
  if ( sum( G == 2 ) > 0 & sum( G == 1 ) > 0 ) {
    a <- which(G == 1, arr.ind = TRUE)
    a <- a[order(a[,1]), ]
    id <- unique(a[, 1])
    b <- as.vector( table(a[, 1]) )
    b <- id[b > 1]
    a2 <- list()
    for (i in b)   a2[[ i ]] <- c( i, a[which(a[, 1] == i), 2] )
    for ( i in b)  {
      arxi <- a2[[ i ]][ 1 ]
      met <- a2[[ i ]][ -1 ]
      ela <- matrix(0, nrow = 1, ncol = 2)
      for ( j in 1:length(met) )  {
        ande <- which( G[ met[j], ] == 2 )
        if (sum(ande) == 0) {
          yp <- cbind(met[j], 0 )
        } else  yp <- cbind(met[j], ande )
        ela <- rbind(ela, yp)
      }
      ela <- cbind(arxi, ela)
      ela <- ela[-1, ]
      dipla <- which( as.vector( table(ela[, 3]) ) > 1 )
      if (length(dipla) > 0) {
        for ( k in 1:length(dipla) ) {
          poia <- Rfast::sort_unique(ela[, 3])[ dipla[k] ]
          poia2 <- ela[ which( ela[, 3] == poia ), ]
          nr <- dim(poia2)[1]
          if (nr > 2) {
            a <- combn(poia2[, 2], 2)
            for ( m in 1:dim(a)[2] ) {
              b <- cbind(i, a[, m], poia2[1, 3])
              if ( G[ b[1, 2], b[2, 2] ] != 1 ) {
                G[ i, b[1, 3] ] <- 2
                G[ b[1, 3], i ] <- 3
                if (prev.cycles) {
                  if ( !is.dag(G) )  {
                    G[ b[1, 3], i ] = 1
                    G[ i, b[1, 3] ] = 1
                  }
                }
              }
            }
            
          } else {
            if ( G[ poia2[1, 2], poia2[2, 2] ] != 1 ) {
              G[ i, poia2[1, 3] ] <- 2
              G[ poia2[1, 3], i ] <- 3
              if (prev.cycles) {
                if ( !is.dag(G) )  {
                  G[ poia2[1, 3], i ] = 1
                  G[ i, poia2[1, 3] ] = 1
                }
              }
            }
          }  
        }
      }
    }
  }  
  G 
}




is.sepset <- function(pair, nd, separ) {
  pair <- sort(pair)
  a <- which( Rfast::colsums( abs( t(separ[, 1:2, drop = FALSE]) - pair ) ) == 0 )
  if ( length(a) >= 1 )  { 
    b <- length( intersect(nd, separ[a, -c(1:2), drop = FALSE]) ) > 0
  } else b <- FALSE
  
  b
} 