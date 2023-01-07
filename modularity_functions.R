modularity_matrix <- function(g, gamma = 1, weighted = TRUE){
  #
  # INPUT:
  #  g,        igraph or matrix object (to work as adjacency matrix)
  #  memb,     community membership
  #  gamma,    (optional) resolution parameter. Default is 1, which gives equivalent to Newman (2004) modularity
  #  weighted, (optional) logical. TRUE for weighted networks
  #     
  # OUTPUT:
  # M, Modularity Matrix
  
  stopifnot(is.igraph(g) | is.matrix(g));
  
  if(is.igraph(g)){                           # get adjacendy matrices
    if(isFALSE(weighted)){
      A <- as.matrix(as_adjacency_matrix(g));
    } else {
      A <- as.matrix(g[]);
    }
  } else {
    A <- g
  }
  
  K <- colSums(A);                            # degrees/strengths
  v <- sum(A);                                # Number of edges (each edge counted twice, 2m)
  M <- A - gamma*outer(K,K,'*')/v;            # modularity matrix
  
  return(M)
}

Q <- function(g, memb, gamma = 1, weighted = TRUE, weights = NULL){
  #
  # INPUT:
  #  g,        igraph or matrix object (to work as adjacency matrix)
  #  memb,     community membership
  #  gamma,    (optional) resolution parameter. Default is 1, which gives equivalent to Newman (2004) modularity
  #  weighted, (optional) logical. TRUE for weighted networks
  #  weights,  user-defined weights. Must be same length as vertices (V(g)). 
  #  if NULL   (default) gets weights from E(g)$weight if present
  #     
  # OUTPUT:
  # Q, Modularity
  
  stopifnot(is.igraph(g) | is.matrix(g));
  
  if(is.igraph(g)){ 
    stopifnot(length(memb) == length(V(g)));
  } else { 
    stopifnot(length(memb) == ncol(g));
  }
  
  if(is.igraph(g)){                           # get adjacendy matrices
    if(isFALSE(weighted)){
      A <- as.matrix(as_adjacency_matrix(g));
    } else {
      if(!is.null(weights)) { E(g)$weight <- weights }
      A <- as.matrix(g[]);
    }
  } else {
    if(!is.null(weights)) { g <- matrix(weights, ncol = ncol(g)) }
    A <- g
  }
  
  K <- rowSums(A);                            # degrees/strengths
  v <- sum(K);                                # Number of edges (each edge counted twice, 2m)
  M <- A - gamma*outer(K,K,'*')/v;            # modularity matrix
  Delta <- outer(memb, memb, '==');           # Delta(M1,M2)
  Q <- (1/v)*sum(M*Delta, na.rm = T);         # modularity
  
  if(is.nan(Q)){ Q <- 0}
  
  return(Q)
}

Q_signed <- function(g, memb, gamma = 1, weighted = TRUE, weights = NULL){
  #
  # INPUT:
  #  g,        igraph or matrix object (to work as adjacency matrix)
  #  memb,     numeric vector with community membership
  #  gamma,    (optional) resolution parameter. Default is 1, which gives equivalent to Newman (2004) modularity
  #  weighted, (optional) logical. TRUE for weighted networks
  #  weights,  user-defined weights. Must be same length as vertices (V(g)). 
  #  if NULL   (default) gets weights from E(g)$weight if present
  #     
  # OUTPUT:
  # Q, Modularity
  #
  # This modularity is that presented on Rubinov & Sporns (2011). Weight-conserving 
  # characterization of complex functional brain networks. Neuroimage 56, 2068.
  
  stopifnot(is.igraph(g) | is.matrix(g))
  
  if(is.igraph(g)){ 
    stopifnot(length(memb) == length(V(g))) 
  } else { 
    stopifnot(length(memb) == ncol(g)) 
  }  
  
  if(is.igraph(g)){                           
    if(isFALSE(weighted)){
      A <- as.matrix(as_adjacency_matrix(g)) 
    } else {
      A <- as.matrix(g[])                     
    }
  } else {
    A <- g
  }
  
  # getting graph with positive weights
  pos <- A*(A > 0);
  vpos <- sum(rowSums(pos));
  
  # getting graph with negative weights
  neg <- -A*(A < 0);
  vneg <- sum(rowSums(neg));
  
  qpos <- Q(pos, memb = memb, weighted = T, weights = weights, gamma = gamma)
  qneg <- Q(neg, memb = memb, weighted = T, weights = weights, gamma = gamma)
  
  qstar <- qpos - (vneg/(vpos+vneg))*qneg
  
  return(qstar)
}
