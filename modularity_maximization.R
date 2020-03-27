
modularity_maximization <- function(g, gamma = 1, weighted = T, class = TRUE){
  #
  # This function is an R implementation of the modularity maximization function for Matlab
  # in the brain connectivity toolbox (https://sites.google.com/site/bctnet/measures/list as modularity_und.m)
  #
  # INPUT:    g           Adjacendy matrix, can an igraph or a matrix object
  #           gamma       (optional) Numeric, gamma resolution (optional)
  #           weighted    (optional) Logical, whether to use matrix's weights (if present) or not
  #           class       Logical. TRUE (default) gives list list with class 'communities'. FALSE gives 
  #                       avector with the community assignments
  #
  # OUTPUT:   list of class 'communities' (if class == TRUE) or just the membership vector
  #           membership  optimal community structure 
  #           modularity  optimized modularity
  #           names       nodes' names
  #           algorithm   'modularity maximization'
  #
  #
  # Note: (from the document in Brain Connectivity toolbox for Matlab)
  #       This algorithm is essentially deterministic. The only potential
  #       source of stochasticity occurs at the iterative finetuning step, in
  #       the presence of non-unique optimal swaps. However, the present
  #       implementation always makes the first available optimal swap and
  #       is therefore deterministic.
  
  
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
  v <- sum(K);                                # Number of edges (each edge counted twice, 2m)
  M <- A - gamma*outer(K,K,'*')/v;            # modularity matrix
  
  Ci <- rep(1, ncol(M));                      # community indices
  cn <- 1;                                    # number of commuities
  U <- as.numeric(c(1,0));                    # vector of unexamined communities
  Mg <- M;
  Ng <- ncol(M);
  ind <- seq(Ng);
  
  while (as.logical(U[1])) {
    #e <- eigen(Mg);                          # eigenvalues and eigen vectors of modularity matrix
    # eigen() reorders columns but leaves rows unaltered
    e <- fda::Eigen(Mg);                      # eigenvalues and eigen vectors of modularity matrix
    n <- which(e$values == max(Re(e$values)));# eigenvector of the maximal positive eigenvalue of M
    v1 <- Re(e$vectors[,n]);
    
    S <- as.numeric(rep(1, Ng));              # vector of ones of length = ncol M
    S[v1 < 0 %in% S] <- -1;                   # negative eigenvectors of highest eigenvalue = -1
    q <- as.numeric(t(S) %*% (Mg %*% S));     # contribution to modularity
    
    if(q > 1e-16){                            # contribution positive: U(1) is divisible
      qmax <- q;                              # maximal contribution to modularity
      Mg[which(diag(Ng) == 1)] <- 0;          # Zeroing diag to enable fine-tuning
      indg <- rep(1, Ng);                     # array of unmoved indices
      Sit = S;
      
      while (suppressWarnings(
        any(indg, na.rm = T)
      )){                                   # iterative fine-tuning
        Qit <- qmax - (4*Sit)*(Mg %*% Sit); # this line is equivalent to 
        qmax <- max(Qit*indg, na.rm = T);
        imax <- which(Qit*indg == qmax);
        Sit[imax] <- -Sit[imax];
        indg[imax] <- NaN;
        if(qmax > q){
          q <- qmax;
          S <- Sit;
        }  
      }
      
      if(abs(sum(S)) == Ng){                  # unsuccessful splitting of U[1]
        U = U[-(1)];                          # in MatLab, the implementationis U(1) = []...
      } 
      else {
        cn <- cn+1;
        Ci[ind[S == 1]] <- U[1];              # splits old U[1] into new U[1] and updated cn value
        Ci[ind[S == -1]] <- cn;
        U <- c(cn, U);                        # ok<AGROW>
      }
    }
    else {                                    # contribution of nonpositive: U[1] is indivisible
      U = U[-(1)];                            # in MatLab, the implementationis U(1) = []...
    }
    
    ind <- which(Ci == U[1]);                 # indices of unexamined community U[1]
    mg <- M[ind,ind];
    if(length(dim(mg)) > 1){
      Mg <- mg - diag(colSums(mg));           # modularity matrix for U[1] #  
    }
    else{                                     # work around single node communities
      Mg <- mg - diag(mg);
    }
    
    Ng <- length(ind);                        # number of vertices in subgraph U[1]
  }
  
  s <- matrix(rep(Ci, ncol(M)), ncol = ncol(M));
  Delta <- (!(s - t(s)));                     # matrix w TRUE for same community, vice-versa
  Q <- Delta*M/v;                             # new modularity matrix
  Q <- sum(Q);                                # maximazed modularity
  
  
  if(isTRUE(class)){
    Ci <- list(membership = Ci,
               modularity = Q,
               names = colnames(A),
               algorithm = 'modularity maximization')
    Ci <- structure(Ci, class = 'communities')  # To fit in other functions from igraph
  }

  
  return(Ci)
}