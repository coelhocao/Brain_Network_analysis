## REQUIRED PACKAGES
if (!('igraph' %in% installed.packages()[,'Package'])){install.packages("igraph")}; require(igraph) #package for network generation and measures
####################################################################################################

participation_coefficient <- function(g, Ci, sym = FALSE, out = 'list', sep_values = T){
  # Participation coefficient is a measure of diversity of intermodular
  # connections of individual nodes. This function was adapted from the 
  # participation_coef_sign function from the Brain connectivity Toolbox 
  # https://sites.google.com/site/bctnet/measures/list 
  # 
  # INPUT:
  #        g          igraph or matrix object picturing a undirected weighted signed graph
  #        Ci         Community assignment vector
  #        sym        Logical. if FALSE(default), applies assymetrical contribution of
  #                   negative and positive coefficients as proposed in the paper
  #                   Rubinov & Sporns (2011). Neuroimage 56, 2068. If TRUE, just 
  #                   their subtraction.
  #
  # OUTPUT:
  #        Ppos      participation coefficients from positive weights
  #        Pneg      participation coefficients from negative weights
  #        P_signed  Signed participation coefficent according to sym arg
  
  stopifnot(out == 'list' | out == 'data.frame')
  stopifnot(is.logical(sym))
  
  if(is.igraph(g)){                           # get adjacendy matrices
    g <- as.matrix(g[]);
  }
  else {
    g <- g
  }
  
  pcoef <- function(g, Ci){
    
    S <- rowSums(g); # strength
    Gc <- (g != 0) %*% (diag(Ci)); # neighbor community affiliation
    Sc2 <- rep(0,ncol(g)); # community-specific neighbors
    
    Sc2 <- rowSums(
      sapply(
        seq(max(Ci)),
        function(i){
          (rowSums(g * (Gc == i))^2);
        }
      )
    )
    
    P = rep(1,ncol(g)) - (Sc2/(S^2));
    P[P %in% NaN] <- 0;
    P[!P] = 0; # p_ind=0 if no (out)neighbors
    
    return(P)
  }
  
  gpos <- g;
  gpos[gpos < 0] <- 0
  
  gneg <- -g
  gneg[gneg < 0] <- 0
  
  Ppos <- pcoef(gpos, Ci);
  Pneg <- pcoef(gneg, Ci);
  
  if(isFALSE(sym)){
    Spos <- rowSums(gpos);
    Sneg <- rowSums(gneg);

    P <- list(Ppos = Ppos, Pneg = Pneg, 
                    P_signed = Ppos - (Sneg/(Spos + Sneg))*Pneg) 
  }
  else{
    P <- list(Ppos = Ppos, Pneg = Pneg, 
                    P_signed = Ppos - Pneg) 
  }
  
  if(isFALSE(sep_values)){ P <- P$P_signed }
  if(out == 'data.frame'){ P <- as.data.frame(P) }
  
  return(P)
}

diversity_coefficient <- function(g, Ci, sym = FALSE, out = 'list', sep_values = T){
  # Shannon-entropy based diversity coefficient measures the diversity of 
  # intermodular connections of individual nodes and ranges from 0 to 1 
  # (the H_signed ranges from -1 to 1). This function was adapted from the 
  # diversity_coef_sign function from the Brain connectivity Toolbox 
  # https://sites.google.com/site/bctnet/measures/list
  # 
  # INPUT:
  #        g          igraph or matrix object picturing a undirected weighted signed graph
  #        Ci         Community assignment vector
  #        sym        Logical. if FALSE(default), applies assymetrical contribution of
  #                   negative and positive-based diversity values as proposed in 
  #                   Rubinov & Sporns (2011). Neuroimage 56, 2068. If TRUE, just 
  #                   their subtraction.
  #
  # OUTPUT:
  #        Hpos      Diversity coefficient based on positive connections
  #        Hneg      Diversity coefficient based on negative connections
  #        H_signed  Signed diversity coefficent according to sym arg
  
  stopifnot(out == 'list' | out == 'data.frame')
  stopifnot(is.logical(sym))
  
  if(is.igraph(g)){                           # get adjacendy matrices
    g <- as.matrix(g[]);
  }
  else {
    g <- g
  }
  
  entropy <- function(g,Ci){
    
    n = ncol(g);
    m <- max(Ci);
    
    # strengths
    S <- rowSums(g);
    # node-to-module degree
    Snm <- matrix(0,n,m);
    # loop over modules
    Snm <- sapply(seq(m), function(i){
      rowSums(as.matrix(g[,Ci==i]));
    })
    Sp <- sapply(seq(m), function(i){ S })
    pnm <- Snm/Sp
    pnm[pnm %in% NaN] <- 0;
    pnm[!pnm] <- 1;
    H = -rowSums(pnm*log(pnm))/log(m);
    
    return(H)
  }
  
  pos <- g*(g>0);
  neg <- -g*(g<0);
  
  Hpos <- entropy(pos, Ci);
  Hneg <- entropy(neg, Ci);
  
  if(isFALSE(sym)){
    Spos <- rowSums(pos);
    Sneg <- rowSums(neg);
    H_signed <- Hpos - (Sneg/(Spos + Sneg))*Hneg;
    H <- list(Hpos = Hpos, Hneg = Hneg, H_signed = H_signed);
    
    # In case isolated regions are included
    if(any(is.nan(H$H_signed))){ H$H_signed[H$H_signed %in% NaN] <- 0 }
  }
  else{
    H <- list(Hpos = Hpos, Hneg = Hneg, H_signed = Hpos - Hneg);
  }
  
  if(isFALSE(sep_values)){ H <- H$H_signed }
  if(out == 'data.frame'){ H <- as.data.frame(H) }
  
  return(H)
}

gateway_coefficient <- function(g, Ci, centtype = c('str', 'bet'), sym = FALSE, out = 'list', sep_values = F){
  #
  # Gateway coefficient (gcoef) is a cariant of participation coefficient. Similar to participation
  # coefficient, gcoef measures the civersity of intermodular connections of indivicual nodes, but
  # this is weighted by how critical these connections are to intermodular connectivity (e.g. if a
  # node is the only connection between it's module and another module, it will have higher gcoef)
  # This function was adapted from the gateway_coef_sign function the Brain connectivity Toolbox 
  # https://sites.google.com/site/bctnet/measures/list
  #
  # INPUT: 
  #        g          igraph or matrix object picturing a undirected weighted signed graph
  #        Ci         Community assignment vector
  #        centtype   Centrality measure to be used. 
  #                   'sth' = strength, 
  #                   'bet' = betweenness
  #        sym        Logical. if FALSE(default), applies assymetrical contribution of
  #                   negative and positive coefficients as proposed in the paper
  #                   Rubinov & Sporns (2011). Neuroimage 56, 2068. If TRUE, just 
  #                   their subtraction.
  # OUTPUT:
  #        GWpos      Gateway coefficient based on positive connections
  #        GWneg      Gateway coefficient based on negative connections
  #        GW_signed  Signed gateway coefficent according to sym arg
  
  stopifnot(out == 'list' | out == 'data.frame')
  stopifnot(is.logical(sym))
  
  if(is.igraph(g)){                           # get adjacendy matrices
    g <- as.matrix(g[]);
  }
  else {
    g <- g
  }
  
  Ci <- match(Ci, unique(Ci)); # Remap module indices to consecutive numbers
  n <- ncol(g); # Number of nodes
  g = g*(diag(n) == 0); #Ensure diagonal is zero
  
  gcoef <- function(g, Ci){
    
    k <- rowSums(g); # strengths
    Gc <- (g != 0) %*% diag(Ci); # Create neighbor community affiliation matrix
    nmod = max(Ci); # Find # of modules
    ks <- matrix(0,n, nmod); # Preallocate space
    kjs <- matrix(0,n,nmod); # Preallocate space
    cs <- matrix(0,n,nmod); # Preallocate space
    
    bet <- function(g) {  # Betweenness function
      inv <- graph.adjacency(g, mode = 'undirected', weighted = T); # makes igraph object
      E(inv)$weight <- 1/E(inv)$weight; # turns weights into distance
      b <- betweenness(inv); #betweenness
      return(b)
    }
    cent <- switch(centtype, # Which centrality measure to use?
                   'str' = rowSums(g), # Node Strength
                   'bet' = bet(g)); # Betweenness centrality
    
    mcn <- 0; # Set max summed centrality per module to 0
    for(i in seq(nmod)){ # for each module
      if(sum(cent[Ci == i]) > mcn){ # If current module has a higher sum
        mcn = sum(cent[Ci == i]); # reassign value
      }
      ks[,i] <- rowSums(g*(Gc == i)); # Compute the total weight of the connections per node to each module
    }
    
    for(i in seq(nmod)){ # For each module 
      if(sum(Ci == i) > 1){ # If there is more than 1 node in a module 
        kjs[Ci == i,] <- as.matrix(rep(1,sum(Ci == i))) %*% 
          as.matrix(t(colSums(ks[Ci == i,]))); # Compute total module-module connections
        kjs[Ci == i,i] <- kjs[Ci == i,i]/2; # Account for redundancy due to double counting within-network work weights
      }
    }
    
    for(i in seq(n)){ # For each node
      if(k[i] > 0){ # If node is connected
        for(j in seq(nmod)){ # For each module
          cs[i,j] <- sum(cent[(Ci*(g[,i] > 0)) == j]); # Sum of centralities of neighbors of a node within a module
        }
      }
    }
    
    ksm <- as.matrix(ks/kjs); # Normalize by total connections
    ksm[kjs == 0] <- 0; # Account for division by 0
    csm <- as.matrix(cs/mcn); # Normalize by max summed centrality
    gs <- as.matrix((1 - (ksm*csm))^2); # Calculated total weighting
    GW <- 1 - rowSums((ks^2)/(k^2)*gs); # Compute gateway coefficient
    GW[GW %in% NaN] <- 0; #account for division by 0
    GW[!GW] <- 0; #Set to 0 if no neighbors
    
    return(GW)
  }
  
  GWpos <- gcoef(g*(g > 0), Ci); # Compute gateway coefficient for positive weights
  GWneg <- gcoef(-g*(g < 0), Ci); # Compute gateway coefficient for negative weights
  
  if(isFALSE(sym)){
    Spos <- rowSums(g*(g > 0));
    Sneg <- rowSums(-g*(g < 0));
    GW_signed <- GWpos - (Sneg/(Spos + Sneg))*GWneg;
    GW <- list(GWpos = GWpos, GWneg = GWneg, GW_signed = GW_signed);
  }
  else{
    GW <- list(GWpos = GWpos, GWneg = GWneg, GW_signed = (GWpos - GWneg));
  }
  
  if(isFALSE(sep_values)){ GW <- GW$GW_signed }
  if(out == 'data.frame'){ GW <- as.data.frame(GW) }
  return(GW)
}

within_community_strength <- function(G, comm, signed = TRUE, normalized = TRUE, sep_values = T){
  #
  # INPUT: igraph object and an object with communities class attr
  #
  # OUTPUT: the within-community strength of each node
  stopifnot(is.logical(signed) | is.logical(normalized) | is.logical(sep_values))
  #stopifnot(length(V(G)) == length(comm))
  if(length(V(G)) != length(comm)){
    print(V(G))
    print(comm)
    print(length(comm))
    stop('length(V(G)) != length(comm)')
  }
  
  subgraphs <- lapply(seq(max(comm)), function(i){
    induced_subgraph(G, V(G)[comm == i])
  })
  
  if(isTRUE(signed)){
    wcs <- lapply(subgraphs, function(i){
      signed_strength(i, out = 'data.frame',
                      normalized = normalized, 
                      sep_values = sep_values)
    })
    wcs <- do.call(rbind, wcs)
  } else {
      wcs <- data.frame(
        wcs = unlist(sapply(subgraphs, function(i){
        strength(i)/length(V(i))
      })))
  }
  
  wcs <-  data.frame(sth = wcs[order(V(G)$name),])
  rownames(wcs) <- V(G)$name
  return(wcs)
}

community_based_metrics <- function(g, comm, signed = T, normalized = T, sep_values = T,
                                    metrics = c('within_community_strength', 'diversity',
                                                'participation', "gateway_bet","gateway_str"), hier_clust = F){
  #
  #
  # consensus, participation, wcs, gateway_bet, gateway_str
  #
  stopifnot(is.igraph(g) | 'communities' %in% attr(comm, 'class'))
  stopifnot(length(E(g)$weight) > 0)
  stopifnot(is.logical(hier_clust))
  
  if(isTRUE(hier_clust)){
    if(is_hierarchical(comm)){
      memb <- comm$merges
      if(isTRUE(sep_values)){
        sep_values = F
        cat('sep_values is set to FALSE in hierarchical clusters for now.\n')
      }
    }
    else{
      memb <- matrix(comm$membership, ncol = 1)
      warning("hier_clust is TRUE but comm$merges doesn't exist. \n
              using comm$membership instead!")
    }
  }
  else{
    memb <- matrix(comm$membership, ncol = 1)
  }
  
  output <- data.frame(name = V(g)$name)

    if('within_community_strength' %in% metrics) {
      output <- data.frame(region = output, wcs = apply(as.matrix(memb), 2, function(i){
        within_community_strength(g, i, signed = signed, normalized = normalized, sep_values = sep_values)
      }))
    }
    if('participation' %in% metrics) { 
      output <- data.frame(output, participation = apply(as.matrix(memb), 2, function(i){
        participation_coefficient(g, i, sep_values = sep_values)
      }))
    }
    if('diversity' %in% metrics) {
      output <- data.frame(output, diversity = apply(as.matrix(memb), 2, function(i){
        diversity_coefficient(g, i, sep_values = sep_values)
      }))
    }
    if('gateway_bet' %in% metrics) { 
      output <- data.frame(output, gate_bet = apply(as.matrix(memb), 2, function(i){
        gateway_coefficient(g, i, centtype = 'bet', sep_values = sep_values)
      }))
    }
    if('gateway_str' %in% metrics) {
      output <- data.frame(output, gate_str = apply(as.matrix(memb), 2, function(i){
        gateway_coefficient(g, i, centtype = 'str', sep_values = sep_values)
      }))
    }

  output <- as.data.frame(output[,2:ncol(output)])
  return(as.data.frame(output))
}
