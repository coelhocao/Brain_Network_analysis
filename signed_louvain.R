if (!('igraph' %in% installed.packages()[,'Package'])){install.packages("igraph")}; require(igraph)
if (!('psych' %in% installed.packages()[,'Package'])){install.packages("psych")}; require(psych)

signed_louvain <- function(g, gamma = NULL, weighted = T, comm = NULL, B = NULL, seed = NULL,
                        mod = c(NULL, 'modularity', 'potts', 'neg_sym', 'neg_asym'), class = TRUE) {
  
  # Fast and accurate multi-iterative implementation of theLouvain community detection algorithm. 
  # This implementation was adapted from the MatLab function community_louvain.m available on 
  # Brain Connectivity toolbox. https://sites.google.com/site/bctnet/
  # 
  # INPUT: 
  #          G,     a directed/undirected weighted/unweighted igraph object or matrix object
  #          gamma, (optional) resolution parameter. gamma > 1 small modules, 0 < gamma > 1, larger modules 
  #          comm,  (optional) initial community affiliation vector
  #          B,     (optional) custom objective matrix (must have the same dimensions as input graph)
  #          mod,   (optional) Objective function type. it can have the following inputs
  #                 'modularity'  Modularity (default)
  #                 'potts'       Potts-model Hamiltonian (for binary networks)
  #                 'neg_sym'     symmetric treatment of negative weights
  #                 'neg_asym',    asymmetric treatment of negative weights (see Rubinov and Sporns (2011))
  #                 custom objective-function matrix (under substitute(), bquote(), expression())
  #                 bquote(B <- (A - gamma*outer(rowSums(A),colSums(A),'*')/s)/s)
  #          class  Logical. TRUE (default) gives list with class 'communities'. FALSE gives 
  #                 a vector with the community assignments
  #
  # OUTPUT:   list of class 'communities' (if class == TRUE) or just the membership vector
  #           membership  optimal community structure 
  #           modularity  modularity value
  #           names       nodes' names
  #           algorithm   chr 'signed_louvain'
  
  #stopifnot(is.igraph(g) | is.matrix(g));

  if(is.igraph(g)){                           
    if(isFALSE(weighted)){
      A <- as.matrix(as_adjacency_matrix(g));
    } else {
      A <- as.matrix(g[]);
    }
  } else {
    if(is.matrix(g)){
      A <- g 
    } else{
      stop('g must be an igraph or matrix object!')
    }
  }
  
  if(!is.null(seed)){
    set.seed(seed);
  }
  if(is.null(gamma)){ gamma = 1 }
  
  n <- ncol(A) #get number of edges
  s <- sum(A)  #get sum of edges
  
  if(mod == 'neg_sym' | mod == 'neg_asym'){
    # getting graph with positive weights
    pos <- A;
    pos[pos < 0] <- 0;
    vpos <- sum(pos);
    Bpos <- pos - gamma*outer(rowSums(pos),colSums(pos), '*')/vpos; #positive modularity
    
    # getting graph with negative weights
    neg <- -A;
    neg[neg < 0] <- 0;
    vneg <- sum(neg);
    if(vneg != 0){
      Bneg <- neg - gamma*outer(rowSums(neg),colSums(neg), '*')/vneg; #negative modularity
    }
    else{
      Bneg <- 0;
    }
  }
  else{
    if(min(A) < -1e-10){
      stop('Input graph/matrix has negative weights.\n
           Specify "neg_sym" or "neg_asym" in mod argument')
    }
  }
  
  if(mod == 'potts' && any(a != as.logical(a))){
    stop('Potts-model Hamiltonian requires a binary network')
  }
  
  if(is.null(mod)){ 
    mod <- 'modularity' 
  }
  
  if(!is.null(B)){
    if(identical(dim(B), dim(A))){
      B <- as.matrix(B); 
    }
    else{
      stop('Graph arg and B arg must have the same dimensions')
    }
  }
  else{
    if(is.language(mod)){
      eval(mod) # Ex: substitute(B <- (A - gamma*outer(rowSums(A),colSums(A),'*')/s)/s), bquote(), expression()
    }
    else{
      if(mod == 'modularity'){
        B <- as.matrix((A - gamma*outer(rowSums(A),colSums(A),'*')/s)/s);
      }
      else{
        if(mod == tolower('potts')){
          B <- as.matrix(A - gamma*(!as.logical(A)));
        }
        else{
          if(mod == 'neg_sym'){
            B <- as.matrix(Bpos/(vpos + vneg) - Bneg/(vpos+vneg));
          }
          else{
            if(mod == 'neg_asym'){
              B <- as.matrix(Bpos/vpos - Bneg/(vpos+vneg));
            }
            else{
              stop('Must choose "modularity", "potts(binary only)", "neg_sym", "neg_asym"\n 
                   or input a valid objective function for modularity matrix B')
            }
          }
        }
      }
    }
  }
  
  if(is.null(comm)){
    M0 <- as.numeric(c(1:n));
  } 
  else{
    if(length(comm) != n){ 
      stop('Length of comm arg must be the same as number of nodes')
    }
    else{
      M0 <- comm
    }
  }
  
  Mb <- match(M0, unique(M0));
  M <- Mb
  
  B <- as.matrix((B + t(B))/2);                     # symmetrize modularity matrix
  Hnm <- matrix(0,n,n);                             # node-to-module degree
  for(i in 1:max(Mb)){                              # loop over modules
    Hnm[,i] <- rowSums(as.matrix(B[,Mb == i]))
  }
  
  Q0 <- -Inf;
  Q <- sum(B[outer(M0, M0, '==')]);                 # compute modularity
  first_iteration <- TRUE;
  flag1 <- TRUE;                                    # explicit flag for first while loop
  
  while(isTRUE(flag1)){
    flag1 <- FALSE;
    flag2 <- TRUE;                                  # flag for within-hierarchy search
    
    while(isTRUE(flag2)){
      flag2 <- FALSE;
      
      for(i in sample(1:n, replace = F)){           # loop over all nodes randomically
        ma <- Mb[i];                                # current module i
        dQ <- Hnm[i,] - Hnm[i,ma] + B[i,i];
        dQ[ma] <- 0;                                # (line above) algorithm condition
        
        max_dQ <- max(dQ);                          # maximal increase in modularity
        mb <- which(dQ %in% max_dQ)[1];             # and corresponding module (the [1] is bc sometimes Length(mb) > 1)
        
        if(isTRUE(max_dQ > 0)){                     # if maximal increase is positive
          flag2 <- TRUE;
          Mb[i] <- mb;                              # reassign module
          
          Hnm[,mb] <- Hnm[,mb] + B[, i];            # change node-to-module strengths
          Hnm[,ma] <- Hnm[,ma] - B[, i];
        }
      }
    }
    
    Mb <- match(Mb,unique(Mb));                     # new module assignments
    M0 <- M;
    if(isTRUE(first_iteration)){
      M <- Mb;
      first_iteration <- FALSE;
    }
    else{
      for(i in 1:n){                                # Loop through initial module assignments
        M[M0 == i] <- Mb[i];                        # assign new modules
      }
    }
    
    n <- max(Mb);                                   # new number of modules (weird), initially n was supposed to be nodes, not modules
    Bl <- matrix(0,n,n);                            # new weighted matrix
    for(i in 1:n){
      for(z in i:n){
        bm <- sum(B[Mb == i, Mb == z]);             # pull weights of nodes in same module
        Bl[i,z] <- bm;
        Bl[z,i] <- bm;
      }
    }
    B <- Bl;
    
    Mb <- seq_len(n);                               # initial module assignments
    Hnm <- B;                                       # node-to-module strength
    
    Q0 <- as.numeric(Q)
    Q <- psych::tr(B)                               # compute modularity sum(diag()) is the same as trace()
    # De-bugging Q - Q0 >1e-10 non-logical (which means Q or Q0 is na or nan)
    if(is.na(Q) | is.nan(Q0) |
       is.na(Q0) | is.nan(Q0)){
      #print(c(Q, Q0))  #this points can be used to look for bugs across loops
      break
    }
    else{
      if(Q > Q0){
        flag1 <- TRUE
        #print(paste('Q = ', Q))#this points can be used to look for bugs across loops
      }
      #print(paste('Q = ', Q))#this points can be used to look for bugs across loops
    }
  }
  
  if(isTRUE(class)){
    M <- list(membership = M,
              modularity = Q,
              names = colnames(A),
              gamma = gamma,
              algorithm = 'signed_louvain')
    M <- structure(M, class = 'communities')        # To fit in other functions from igraph
  }
  return(M)
}
