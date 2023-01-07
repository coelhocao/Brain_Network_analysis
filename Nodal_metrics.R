#############################################################################################
# These codes have been costumly designed for brain network analysis of 1-2 time point data,
# in our case, acquired by immunolabelling of activity-dependent proteins (i.e. c-fos, Arc)
#
# Many functions here have been largely based on some of the codes shared by Justin Kenney,
#which can be found on this link https://github.com/jkenney9a/Networks
#
# NOTE: It is intended mostly for undirected network data
# Most funcitons won't work in structural brain data or networks which support directed graphs 
# or different upper/lower matrix triangles.
#
# NOTE: Designed mostly for comparisons between single-time point graphs or two-time points graphs
# won't work on multiple time points graphs
#
# LAST UPDATED: Jun 11, 2020
#
#
# Cesar A O Coelho
# cebacio@gmail.com
#

if (!('dplyr' %in% installed.packages()[,'Package'])){install.packages("dplyr")}; require(dplyr)
if (!('igraph' %in% installed.packages()[,'Package'])){install.packages("igraph")}; require(igraph) #package for network generation and measures 

#####################################################################################
#######################          Nodal Measures          ############################
#####################################################################################

signed_strength <- function(g, sym = FALSE, out = 'list', normalized = T, sep_values = F){
  #
  # INPUT: 
  #        g          igraph object or symmetric matrix
  #        sym        Logical. if FALSE(default), applies assymetrical contribution of
  #                   negative and positive coefficients as proposed in the paper
  #                   Rubinov & Sporns (2011). Neuroimage 56, 2068. If TRUE, just 
  #                   their subtraction.
  #        out        whether the output should be a 'list' or 'data.frame'
  # OUTPUT: 
  #        sth        Positive strength, negative strength, and normalized strength
  
  stopifnot(out == 'list' | out == 'data.frame')
  stopifnot(is.logical(sym))
  
  stopifnot(is_igraph(g) | is.matrix(g))
  
  if(is_igraph(g)){ 
    A <- as.matrix(g[]) 
  } else {
    A <- g
  }
  if(max(A) > 1 | min(A) < -1){
    warning('Not a correlation matrix. Output will be beyond [-1,1] range.\n
            Consider using sym = T and normalized = T')
  }
  
  n <- ncol(A);
  Spos <- rowSums(A*(A > 0));  # positive strength
  
  Sneg <- rowSums(-A*(A < 0)); # negative strength
  
  nSpos <- Spos*(1/(n-1));
  nSneg <- Sneg*(1/(n-1));
  
  if(isFALSE(sym)){
    signed <- nSpos - (Sneg/(Spos + Sneg))*nSneg;  # signed strength according to Rubinov & Sporns (2011)
    signed[signed %in% NaN] <- 0;
    signed[!signed] <- 0;
    
    sth <- list(Spos = nSpos, Sneg = nSneg, signed_strength = signed)
  }
  else{
    if(isTRUE(normalized)){
      sth <- list(Spos = nSpos, Sneg = nSneg, signed_strength = (nSpos - nSneg)/(nSpos + nSneg))
    }
    else{
      sth <- list(Spos = Spos, Sneg = Sneg, signed_strength = (nSpos - nSneg)) 
    }
  }
  
  if(isFALSE(sep_values)){ 
    sth <- sth$signed_strength
    names(sth) <- colnames(A)
  }
  
  if(out == "data.frame"){ sth <- data.frame(sth) }
  return(sth)
}

leverage <- function(g, weighted = TRUE, signed = TRUE){
  #
  # Calculates the Leverage centrality as defined in  Joyce et al (2010). (PLoS ONE 5(8):e12200)
  # code mostly taken from igraph wiki igraph.wikidot.com/r-recipes#toc10
  #
  #
  # Input: 
  #       g          igraph object
  #       weighted   Logical. whether it should (TRUE, default) use the weight attribute from the 
  #                  igraph object (if any), not (FALSE).
  #       signed     Logical. whether it should (TRUE, default) also use use the negative weights  
  #                  or not (FALSE).
  # 
  # CAUTION: The options of putting weighted and signed in the leverage calculation are not 
  #          predicted in the paper, and were not fully evaluated in regards to its behavior.
  #
  # Output: 
  #        Leverage centrality of each node
  
  if(isFALSE(weighted) & isFALSE(signed)){
    k <- degree(g)
  } else{
    if(isTRUE(weighted) & isFALSE(signed)){
      k <- strength(g)
    } else{
      if(isTRUE(weighted) & isTRUE(signed)){
        k <- signed_strength(g)
      }
    }
  }
  n <- vcount(g)
  
  lev <- sapply(1:n, function(v) { base::mean((k[v]-k[neighbors(g,v)]) / (k[v]+k[neighbors(g,v)])) })
  
  #when k = 0, leverage = NA or NaN. As this can compromise permutations,
  #and both negative and positive values are meaningful, the closest interpretation is -1
  lev[lev %in% NA] <- -1
  lev[lev %in% NaN] <- -1
  names(lev) <- V(g)$name
  
  return(lev) 
}


Nodal_efficiency <- function(G, normalized=FALSE) {
  # Input: 
  #       G             igraph object. Weights will be used if attr weight exists.
  #       normalized    Logical. Whether (TRUE) or not (FALSE/default) to normalize the value
  #                     (where max value = 1)
  #
  # Output: 
  #       Nodal efficiency of each node 
  #       (i.e, average inverse shortest path length from individual node to all other nodes)
  
  if(ecount(G) > 0){
    E(G)$weight <- 1/abs(E(G)$weight) #calculate distance matrix from correlation matrix b/c
  }
  eff <- 1 / distances(G)
  
  eff[!is.finite(eff)] <- 0
  out <- colSums(eff) / (vcount(G) - 1)
  
  if(normalized == TRUE){
    out <- out / max(out)
  }
  return(out)
}

efficiency_loss <- function(g){
  #
  # Calculates the loss of Global Efficiency (GE) and mean Local Efficiency (LE) after removal of each node - 
  # at a time. Recursively removes a node, calculates GE and LE relative to original graph's GE and LE.
  # NOTE: Takes a LONG TIME. it seems it's not optimized. Beware when using with resampling procedures.
  #
  # Input: 
  #       g   igraph object
  #
  # Output: 
  #         dataframe with Loss of GE and LE after node removal
  
  rem <- function(net,r) {
    # Remove vids vertices by random
    vids<-V(net)[r]
    net <- delete.vertices(net,vids)
    return(net)
  }
  
  rem_calc <- function(g,r){
    
    
    gr <- rem(g,r)
    grgef <- gl_eff(gr)
    grlef <- mean(local_eff(gr))
    
    output <- c(grgef, grlef)
    return(output)
  }
  
  m <- as.data.frame(t(sapply (1:vcount(g), FUN= function (x) rem_calc(g,x))))
  geg <- gl_eff(g)
  leg <- mean(local_eff(g))
  if(geg == 0){ m[,1] <- rep(0, length(m[,1]))}
  else{
    m[,1] <- 1 - (m[,1]/geg)
  }
  if(leg == 0){ m[,2] <- rep(0, length(m[,1]))}
  else{
    m[,2] <- 1 - (m[,2]/leg)
  }
  colnames(m) <- c("gleff0/gleffi", "leff0/leffi")
  
  return(m)
}

Nodal_measures <-function(G, normalized=FALSE, efficiency_loss= FALSE, local_eff=FALSE,
                          metrics = c("degree", "strength", "signed_strength", "eigenvector","betweenness",
                                      "burt","closeness", "leverage", "transitivity", "Nodal_eff"), 
                          out = 'list', sep_values = T, sym = F){
  #
  # Generalized function that measures all defined nodal measures with the option to select the desired ones?
  #
  # Input: 
  #       G                 igraph object
  #       normalized        Logical. Whether or not (FALSE, default) to normalize of the values by the largest
  #       efficiency_loss   Logical. Whether or not (FALSE, default) to calculate efficiency loss. 
  #                         It makes use of weights if they exist. For unweighted calculatation, the igraph 
  #                         input must be unweighted. **CAUTION**: efficiency loss currently takes a long time.
  #                         Be careful when using it with resampling functions. 
  #       local_eff         Logical. Whether or not (FALSE, default) to calculate mean local efficiency. 
  #       metrics           Vector indicating which nodal metrics to be measured. The current possible options are:
  #                         c("degree", "strength", "signed_strength", "eigenvector","betweenness",
  #                         "burt","closeness", "leverage", "transitivity", "Nodal_eff")
  #
  # Output: 
  #       A dataframe of centrality measures selected
  #
  # DISCLOSURE: The centralities accepted in this function may be appropriated for different networks.
  #             Ex. metrics based on distance, rather than correlation, like closeness, betweenness,
  #             efficiency, are supposed to be applied in networks where links represent structural connections
  #             like anatomical brain networks, and may be not appropriated to be applied on networks
  #             in which edges represent probabilistics of interaction (i.e. correlation).
  
  if(ecount(G) == 0){
    zeros <- rep(0, vcount(G))
    output <- as.data.frame(cbind("degree"=zeros, "strength"=zeros, "signed_strength"=zeros, "eigenvector"=zeros,"betweenness"=zeros, 
                                  "burt"=zeros,"closeness"=zeros, "leverage" = zeros, "transitivity"=zeros,
                                  "Nodal_eff" = zeros), row.names=V(G)$name)
    if(local_eff == TRUE){output$local_eff = zeros}
    if(efficiency_loss == TRUE){
      output$geff_loss = zeros 
      output$leff_loss = zeros
    }
    return(output)
  }
  
  if(any(metrics %in% 'betweenness' |
     metrics %in% 'burt' |
     metrics %in% 'closeness' |
     metrics %in% 'transitivity' | 
     metrics %in% 'Nodal_eff')){

    G.pos <- G
    E(G.pos)$weight <- abs(E(G.pos)$weight) #Positive numbers are necessary for betweenness and closeness
    E(G.pos)$weight <- 1/abs(E(G.pos)$weight) #bet and clo are based on distance graphs/matrices, not weights. here, dist= 1 - corr
  }

  #Centrality measures
  output <- data.frame(V(G)$name)
  if('degree' %in% metrics){ output <- data.frame(output, degree = igraph::degree(G, normalized=normalized)) }
  if('strength' %in% metrics){  output <- data.frame (output, strength = strength(G))  }
  if('signed_strength' %in% metrics){  output <- data.frame (output, signed_strength(G, out = out, sep_values = sep_values, sym = sym, normalized = normalized))  }
  if('eigenvector' %in% metrics){  output <- data.frame (output, eigenvector = evcent(G)$vector)  }
  if('betweenness' %in% metrics){  output <- data.frame (output, betweenness = betweenness(G.pos, normalized=normalized, weights = E(G.pos)$weight))  }
  if('closeness' %in% metrics){  output <- data.frame (output, closeness = closeness(G.pos, normalized=normalized))  }
  if('burt' %in% metrics){  output <- data.frame (output, burt = constraint(G))  }
  if('leverage' %in% metrics){  output <- data.frame (output, leverage = leverage(G))  }
  if('transitivity' %in% metrics){  output <- data.frame (output, transitivity = transitivity(G,type="local", isolates='zero'))  }
  if(tolower('Nodal_eff') %in% metrics){  output <- data.frame (output, Nodal_eff = Nodal_efficiency(G, normalized = normalized))  }
  
  #when degree = 0, burt = NaN. Set it to a value higher than max since lower values indicate more centrality
  if('burt' %in% metrics) {output$burt <- output$burt[!is.finite(output$burt)] <- max(!is.finite(output$burt)) + 0.1 }
  
  output <- output[,2:ncol(output)]
  # include nodal local eff (different from nodal_eff) in the dataframe
  if(local_eff == TRUE){
    
    leff <- local_eff(G)
    output$local_eff <- leff
  }
  
  #include loss of efficiency (global and local) after that node's removal (different from nodal_eff and local_eff) 
  if(efficiency_loss == TRUE){
    
    loss <- efficiency_loss(G)
    output$geff_loss <- loss[,1]
    output$leff_loss <- loss[,2]
  }
  
  return(output)
}
