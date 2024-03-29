if (!('Matrix' %in% installed.packages()[,'Package'])){install.packages("Matrix")}; require(Matrix)
source('/home/cesar/Dropbox/R/Network sufficiency/Modular_codes/signed_louvain.R')


agreement <- function(comms){
  #
  # This function was adapted from the function agreeement in the Brain connectivity toolbox
  # (https://sites.google.com/site/bctnet/measures/list as agreement.m)
  #
  # INPUT:  
  #         comms    List, matrix or dataframe of community partitions. Matrices 
  #                  and dfs should be [vertex x partitions]
  # OUTPUT: 
  #         d        Matrix of probabilities for community co-assignment
  #
  
  if(is.list(comms)){
    ci <- as.matrix(sapply(comms, function(i){
      as.numeric(membership(i))
    }))
  }
  else{
    if(is.matrix(comms) || is.data.frame(comms)){
      ci <- as.matrix(comms)
    }
    else{
      stop('comms arg must be a list, a matrix or a dataframe variable')
    }
  }
  
  n <- ncol(ci);
  r <- nrow(ci);
  
  ind <- numeric()
  ind <- as.matrix(as.data.frame(sapply(seq(ncol(ci)), function(x){
    sapply(seq(max(ci[,x])), function(y){
      ind <- as.matrix(as.numeric(ci[,x] == y), nrow = r)
    })
  }, simplify = F)))
  
  d <- ind %*% Conj(t(ind))
  d <- d/n
  d <- d*(!diag(ncol(d)))
  
  return(as.matrix(d))
}

permModel <- function(comms){
  #
  # Computes co-classification probabilities using Permutation Model
  # This function was translated from the Matlab Hierarchical Consensus Clustering
  # Jeub L, et all(2018). Scientific Reports 8:3259 | DOI:10.1038/s41598-018-21352-7
  #
  # INPUT: 
  #        S    Community partition matrix[nodes x partitions]
  #
  # OUTPUT:
  #        mu   Vector of probabilities, where 'mui[j]' is the probability for
  #             a pair of nodes to be co-assigned in partition 'j' under the 
  #             permutation model

  if(is.list(comms)){
    S <- as.matrix(sapply(comms, function(i){
      as.numeric(i$membership)
    }))
  }
  else{
    if(is.matrix(comms) || is.data.frame(comms)){
      S <- as.matrix(comms)
    }
    else{
      stop('comms arg must be a list, a matrix or a dataframe variable')
    }
  }
  
  N <- nrow(S);
  L <- ncol(S);
  A <- max(S);

  if(N == 1){
    mu <- matrix(0, N, L);
  }
  else{
    sizes <- matrix(0, A, L);
    for(i in seq(L)){
      G <- Matrix::sparseMatrix(i = seq(N), j = S[,i], x= 1, dims = c(N,A));
      sizes[,i] <- colSums(G);
    }
    
    mu <- colSums((sizes * (sizes -1))/(N*(N - 1)));
  }
  
  return(mu)
}


localPermModel <- function(S){
  #
  # Computes co-classification probabilities using local Permutation Model
  # This function was translated from the Hierarchical Consensus Clustering described in
  # Jeub L, et all(2018). Scientific Reports 8:3259 | DOI:10.1038/s41598-018-21352-7
  #
  # INPUT: 
  #        S    Community partition matrix[nodes x partitions]
  #
  # OUTPUT:
  #        mu   Matrix of probabilities, where 'mui[i,j]' is the probability for
  #             node 'i' to be co-assigned with any other node in partition 'j'
  #             given that the community assignment of 'i' remains fixed.
  
  if(is.list(S)){
    S <- as.matrix(sapply(S, function(i){
      as.numeric(i$membership)
    }))
  }
  else{
    if(is.matrix(S) || is.data.frame(S)){
      S <- as.matrix(S)
    }
    else{
      stop('comms arg must be a list, a matrix or a dataframe variable')
    }
  }
  
  N <- nrow(S);
  L <- ncol(S);
  A <- max(S);
  mu <- matrix(0,N,L);
  
  if(N > 1){
    sizes <- matrix(0, A, L);
    for(i in seq(L)){
      G <- Matrix::sparseMatrix(i = seq(N), j = S[,i], x= 1, dims = c(N,A));
      sizes[,i] <- colSums(G);
    }
    
    for(i in seq(L)){
      mu[,i] <- (sizes[S[,i], i] - 1)/(N - 1);
    }
  }
  
  return(mu)
}

sampleApprox <- function(mu, alpha = 0.05, R = 10000){
  #
  # Computes the sample-approximated null model
  # This model was translated from the Hierarchical consesus Clustering describe in 
  # Jeub L, et all(2018). Scientific Reports 8:3259 | DOI:10.1038/s41598-018-21352-7
  #
  # Computes a sampled approximation to the alpha significance level of the 
  # Poisson-Binomial distribution with parameters mu using R samples for each partition
  #
  # INPUT: 
  #         mu     Parameters for the Poisson-Binomial distribution (PBD). if 'nrow(mu) == 1', 
  #                the PBD is assumed to be the same for each node, otherwise each row of 'mu'
  #                specifies the parameters for the PBD for the corresponding node.
  #         alpha  Confidence level to compute (between o and 1)
  #         R      Number of samples to use to approximate the distribution (default = 10000)
  # OUTPUT:
  #         P      Null model at significance level 'alpha'. this is an
  #                N x N matrix where 'nrow(mu) == N'
  #
  # When the input 'size(mu,1)==1', P will be a scalar. This is the case for the permutation 
  # model. If 'size(mu,1)>1', the function computes the confidence level for each row of mu. 
  # The output 'P(i,j)' is the minimum of the confidence level between rows 'i' and 'j'.
  
  N <- nrow(mu);
  L <- ncol(mu);
  p <- matrix(0,N,1);
  r <- matrix(runif(R*L), R, L);
  
  for(i in seq(N)){
    if(sum(mu[i,] - mu[i,]^2) == 0){
      p[i] <- sum(mu[i,])/L;
    }
    else{
      p[i] <- quantile(rowSums(r < mu[i,]), probs = alpha)/L
    }
  }
  
  P <- matrix(0,N,N);
  for(i in seq(N)){
    for(j in seq(N)){
      P[i,j] <- min(p[i], p[j])
    }
  }
  return(P)
}

normalApprox <- function(mu, alpha = 0.05){
  #
  # Computes the Poisson approximate null model
  # This model was translated from the Hierarchical consesus Clustering describe in 
  # Jeub L, et all(2018). Scientific Reports 8:3259 | DOI:10.1038/s41598-018-21352-7
  #
  # Computes a normal approximation to the alpha-confidence level of the 
  # Poisson-Binomial distribution with parameters mu
  #
  # INPUT: 
  #         mu     Parameters for the Poisson-Binomial distribution (PBD). if 'nrow(mu) == 1', 
  #                the PBD is assumed to be the same for each node, otherwise each row of 'mu'
  #                specifies the parameters for the PBD for the corresponding node.
  #         alpha  Confidence level to compute (between o and 1)
  # OUTPUT:
  #         P      Null model at significance level 'alpha'. this is an
  #                N x N matrix where 'nrow(mu) == N'
  #
  # When the input 'size(mu,1)==1', P will be a scalar. This is the case for the permutation model. 
  # If 'size(mu,1)>1', the function computes the confidence level for each row of mu. 
  # The output 'P(i,j)' is the minimum of the confidence level between rows 'i' and 'j'. The normal
  # approximation takes into account that there is no randomness when 'mu==1' and the output is clipped to [0,1].
  
  N <- nrow(mu);
  L <- ncol(mu);
  p <- matrix(0, N, 1);
  
  for(i in seq(N)){
    a <- mu[i,] == 1;
    k_low <- sum(a);
    
    sigma <- sqrt((sum(mu[i,] - mu[i,]^2)));
    mui <- sum(mu[i,!a])
    if(sigma == 0){
      p[i] <- (k_low + mui)/L;
    }
    else{
      p[i] <- min(max(0,(k_low + qnorm(alpha, mui,sigma))/L),1);
    }
  }
  
  P <- matrix(0, N, N);
  for(i in seq(N)){
    for(j in seq(N)){
      P[i,j] <- min(p[i],p[j]);
    }
  }
  
  return(P)
}

ConsenClust <- function(g = NULL, Ci, perm = localPermModel, distrib = sampleApprox, alpha = 0.05, tau = 0.5, 
                      R = 1000, convergence = TRUE, ALG = signed_louvain, ...){
  #
  # Computes the consensus clustering based on a network's community partition matrix [N x M], where N is
  # the amount of nodes and M different (near degenerate) partitions. An agreement matrix is calculated based on
  # the different partitions. Then a null distribution of partitions is calculated through localPermModel or 
  # PermModel permutations. Then a sample-based or normal-based distribution approximation is used to draw an 
  # agreement matrix under random Null Hypothesis. Then the null Hypothesis agreement matrix (multiplied by alpha) 
  # is subtracted from the empirical one. A community partition algorithm is then used to repartition the resultant
  # agreement matrix R times and the procedure repeats until the matrix converges to a binary matrix where 1s show 
  # co-clustering. Alternatively, this algorithm can be ran only once, instead of recurssively (converge = FALSE), 
  # or be based on a fixed threshold (tau).
  # This algorithm was adapted from 3 papers
  # Jeub L, et al (2018). Scientific Reports 8:3259 | DOI:10.1038/s41598-018-21352-7
  # Betzel RF, et al (2020). NeuroImage 213 (2020) 116687 | DOI: 10.1016/j.neuroimage.2020.116687
  # Lancichinetti A. & Fortunato S., (2012). Scientif Reports. 2, 336 | DOI: 10.1038/srep00336
  #
  #
  # INPUT:
  #        g            (optional) igraph or matrix object. if provided, the output is an
  #                     object of class 'communities' (see below)
  #        Ci           Community partition matrix[nodes x partitions]
  #        tau          Threshold which controls the resolution of the reclustering
  #        R            Numer of repetitions that the clustering algorithm will run
  #        perm         Whether to use localPermModel(default), permModel null model 
  #                     hypothesis (replaces tau), or NULL for a tau-based thresholding
  #        distrib      Approximation for Null Model distribution:
  #        convergence  Logical. repeats the Consensus Clustering until convergence if TRUE (default).
  #                     If FALSE, is performs the B{i,j} = A{i,j} - P{i,j} only once, as described in
  #                     Betzel et al. NeuroImage 213 (2020) 116687
  #        ALG          Function of the community assignemt algorithm used (signed_louvain is default)
  #        ...          Additional arguments are passed to ALG function
  #
  # OUTPUT:   list of class 'communities' (if g is provided) or just the membership vector
  #           membership  optimal community structure 
  #           modularity  optimized modularity
  #           names       nodes' names
  #           algorithm   chr 'Consensus Clustering'
  
  d <- agreement(Ci);
  n <- ncol(d);
  flg = 1;
  
  if(isFALSE(convergence)){
    # As in Betzel 2019 (still on bioRxiv)
    dt <- d - distrib(perm(Ci), alpha = alpha);
  }
  else{
   
    while(flg == 1){
      
      flg <- 0;
      
      if(is.null(perm)){
        dt <- d * (d >= tau) * (!diag(n));
      }
      else{
        # this line does the empirical co-assignment probability - used Null Hypothesis distribution
        dt <- d - distrib(perm(Ci), alpha = alpha);
      }
      
      if(length(dt[dt != 0]) == 0){
        ciu <- seq(n);
      }
      else{
        Ci <- matrix(0, nrow = n, ncol = R);
          Ci <- sapply(seq(R), function(i){
            Ci[,i] <- ALG(dt, class = FALSE, ...)
          })
        Ci <- relabel_partitions(Ci);
        ciu <- unique_partitions(Ci);
        nu <- ncol(ciu);
        if(nu > 1){
          flg <- 1; 
          d <- agreement(as.matrix(Ci));
        }
      }
    }
    ciu <- as.numeric(ciu)   
  }
  
  if(isFALSE(convergence)){
    ciu <- modularity_maximization(dt);
    if(!is.null(g)){
      ciu <- list(
        membership = ciu$membership,
        modularity = Q_signed(g, ciu$membership),
        names = V(g)$name,
        algorithm = 'Consensus Clustering and modularity maximization' # would love to specify ALG here
      )
      ciu <- structure(ciu, class = 'communities')  # To fit in other functions from igraph 
    }
  }
  else{
    if(!is.null(g)){
      ciu <- list(
        membership = ciu,
        modularity = Q_signed(g,ciu),
        names = V(g)$name,
        algorithm = 'Consensus Clustering with signed louvain'
        
      ) 
      ciu <- structure(ciu, class = 'communities')  # To fit in other functions from igraph
    }
    else{
      
    }
  }
  
  return(ciu)
}

relabel_partitions <- function(ci){
  #
  # INPUT :
  #         ci     Community partition matrix[nodes x partitions]
  # OUTPUT: 
  #         cinew  relabelled partitions (to be the same and start at 1)
  
  
  n <- nrow(ci); m <- ncol(ci);
  cinew <- matrix(0, nrow = n, ncol = m);
  for(i in seq(m)){
    a <- ci[,i];
    d <- numeric(length(a));
    count <- 0;
    
    while(sum(d != 0) < n){
      count <- count + 1;
      ind <- which(a != 0)[1];
      tgt <- a[ind];
      rep <- a == tgt;
      d[rep] <- count;
      a[rep] <- 0;
    }
    
    cinew[,i] <- d
  }
  return(cinew)
}

unique_partitions <- function(ci){
  ci <- relabel_partitions(ci);
  ciu <- NULL;
  count <- 0;
  a <- seq(ncol(ci));
  
  while(ncol(ci) > 0){
    count <- count +1;
    tgt <- ci[,1];
    ciu <- cbind(ciu,tgt);
    dff <- colSums(abs(ci - tgt)) == 0; 
    ci <- cbind(ci[,!dff]);
    a <- a[!dff];
  }
  return(ciu)
}