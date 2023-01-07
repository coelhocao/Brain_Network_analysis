source('/home/cesar/Dropbox/R/Network sufficiency/Modular_codes/signed_louvain.R')
source('/home/cesar/Dropbox/R/Network sufficiency/Modular_codes/modularity_functions.R');

if (!('doParallel' %in% installed.packages()[,'Package'])){install.packages("doParallel")}; require(doParallel)
#A <- as.matrix(read.csv('/home/cesar/Dropbox/R/Matlab_toolbox_codes/A_matlab.csv', row.names = 1))

gammaRange <- function(A, mod_fun = modularity_matrix, optimizer = signed_louvain,
                       initial_guess = 1, samples = 20, seed = 1987){
  # gammaRange Compute range of gamma values.
  #
  # Syntax
  #__________________________________________________________________________
  #
  #   [gamma_min,gamma_max]=gammaRange(A)
  #
  #   [gamma_min,gamma_max]=gammaRange(__,Name,Value)
  #
  #
  # Description
  #__________________________________________________________________________
  #
  #   [gamma_min,gamma_max]=gammaRange(A) computes range of 'gamma'
  #       values that result in non-trivial partitions.
  #
  #   [gamma_min,gamma_max]=gammaRange(__,Name,Value) additionally customizes
  #       the behavior of the function by e.g. using a different algorithm to
  #       optimize modularity or using a different modularity-like quality
  #       function.
  #
  #
  # Input Arguments
  #__________________________________________________________________________
  #
  #   A -- Adjacency matrix of the network
  #
  #
  # Name-Value Pair Arguments
  #__________________________________________________________________________
  #
  # Parameter names can be abbreviated and are not case sensitive.
  #
  #   'Optimizer' -- Function for finding "optimal" partitions for a
  #                  modularity-like quality function with modularity matrix
  #                  'B'.
  #                  @(B) iterated_genlouvain(B,[],0,1,'moverandw') (default) |
  #                  function handle
  #
  #   'Modularity' -- Function to compute modularity matrix of the form
  #                   'B=modularity(A,gamma)' where 'A' is an adjacency
  #                   matrix and 'gamma' is a resolution parameter. Note that
  #                   the function assumes that 'B=(A+A')/2-gamma*P' for some
  #                   matrix 'P'.
  #                   @modularity (default) | function handle
  #
  #   'InitialGuess' -- Initial value for 'gamma_min' search.
  #                      1 (default)| scalar
  #
  #   'Samples' -- Number of partitions to sample to test for 'gamma_min' at
  #                proposed value. More samples should result in more
  #                accurate values.
  #
  #
  # Output Arguments
  #__________________________________________________________________________
  #
  #   gamma_min -- Upper bound for smallest value of 'gamma' for which the
  #                network splits into communities.
  #
  #   gamma_max -- Smallest value of 'gamma' for which the network is split
  #                into singleton communities.
  #
  #
  # Implementation
  #__________________________________________________________________________
  #
  # 'gamma_max' is simply the largest value of gamma for which there exist
  # ferromagnetic interactions in the modularity matrix and is easy to
  # compute directly based on the adjacency and modularity matrix.
  #
  # 'gamma_min' does not have a closed-form solution and needs to be
  # approximated numerically. This function uses an iterative algorithm that
  # exploits the linearity of modularity as a function of gamma for a fixed
  # partition. This means that we can directly compute the minimum value of
  # gamma for which a given partition is better than the trivial partition.
  # The iterative algorithm proceeds by first estimating 'gamma_min' using a
  # small sample of partitions. We then sample a new set of patitions using
  # 'gamma=gamma_min-epsilon' (to ensure the previous partitions are strictly
  # non-optimal) and use the new sample to update 'gamma_min' and repeat. The
  # algorithm stops once the new sample consists only of the trivial
  # partition.
  #
  # OBS1 [Cesar]: On R implementation,treating A as sparse isn't very helpful 
  # since multiple functions will not run it as dgCMatrix object from Matrix package.
  # The only solution was to run as.matrix(A) as argument, which makes
  # the output not sparse. Since the sparceness isn't maintained, it helps little
  # for the memory saving.
  #
  #
  # See Also eventSamples, exponentialSamples, hierarchicalConsensus
  
  # Version: 1.1.1
  # Date: Thu  8 Mar 2018 15:34:46 CET
  # Author: Lucas Jeub
  # R implementation by Cesar Coelho
  # Date Jan 2021

  NUM_TOL <- 1e-8;
  
  A <- as(A, 'CsparseMatrix') 
  # find smallest value of gamma necessary to make all interactions negative
  AT <- (as.matrix(A) + Conj(t.default(as.matrix(A))))/2;
  PT <- AT - mod_fun(g = as.matrix(A), gamma = 1);
  gamma_max <- max(div_0(AT, PT));
  # find the largest value of gamma necessary to make all interactions
  # non-negative (necessary to handle the case where A has negative entries)
  g_min <- min(div_0(AT,PT));
  
  # compute partition for gamma=g_min to handle the case where the network is not
  # connected
  B0 <- mod_fun(as.matrix(A), g_min);
  S0 <- optimizer(B0, mod = 'neg_asym', class = FALSE);
  G <- as.matrix(Matrix::sparseMatrix(i = seq(S0), j = S0, x =1));
  a0 <- psych::tr(t(G) %*% AT %*% G);
  p0 <- psych::tr(t(G) %*% PT %*% G);
  
  # make sure initial guess is valid
  gamma_min = Inf;
  set.seed(seed)
  while(gamma_min >= Inf | is.nan(gamma_min)){
    
    S <-matrix(0,nrow(A), samples)
    
    ##############################################################
    #setup parallel backend to use many processors
    cores=detectCores()
    cl <- makeCluster(cores[1]-1) #not to overload your computer
    registerDoParallel(cl)
    
    for(i in seq(samples)){ #parallelize this??
      S[,i] <- optimizer(mod_fun(as.matrix(A), initial_guess - NUM_TOL), mod= 'neg_asym', class = FALSE)
    }
    #stop cluster
    stopCluster(cl)
    ##############################################################

    gamma_min <- gamma_min_bound(AT, PT, S, a0, p0);
    initial_guess <- min(c(2*initial_guess, gamma_max-NUM_TOL));
    if(is.nan(gamma_min)){ #In R, min(Inf, NaN) is NaN, not Inf like in matlab
      warning('initial guess for "gamma_min" did not split the network')
      #Maybe put an user input here for initial_guess not to break the function when error comes
    }
  }
  
  #update using convex hull idea
  gamma_min_new <- gamma_min;
  gamma_min <- Inf;
  set.seed(seed)
  while(gamma_min_new < gamma_min | is.nan(gamma_min_new)){
    if(!is.nan(gamma_min_new)){
      gamma_min <- gamma_min_new; 
    }
    
    ##############################################################
    #setup parallel back end to use many processors
    cores=detectCores()
    cl <- makeCluster(cores[1]-1) #not to overload your computer
    registerDoParallel(cl)
    
    for(i in seq(samples)){
      S[,i] <- optimizer(mod_fun(as.matrix(A), gamma_min-NUM_TOL), mod = 'neg_asym', class = FALSE);
    }
    
    #stop cluster
    stopCluster(cl)
    ##############################################################
    
    gamma_min_new <- gamma_min_bound(AT, PT, S, a0, p0);
  }
  
  return(c(gamma_min = gamma_min, gamma_max = gamma_max))
}


div_0 <- function(A, B){
  # div_0 pairwise division such that 0/0=0
  ind <- which(A != 0)
  A[ind] <- A[ind]/B[ind];
  A[A %in% NaN] <- 0
  return(A)
}

gamma_min_bound <- function(A, PT, S, a0, p0){
  
  gamma_min = Inf;
  for(i in 1:10){#ncol(S)){
    u <- unique(S[,i]);
    e <- match(S[,i], unique(S[,i]));
    if(length(u) > 1){
      G <- Matrix::sparseMatrix(i = seq(nrow(S)), j = e, x=1, dims = c(nrow(S), length(u)));
      a <- psych::tr(as.matrix(t(G) %*% A %*% G));
      p <- psych::tr(as.matrix(t(G) %*% PT %*% G));
      gamma_min <- min(gamma_min, (a0 - a)/(p0 - p))
      #if(is.nan(gamma_min)) print(i)     # This is a bug checker for when min(gamma_min, (a0 - a)/(p0 - p)) = NaN
    }
  }
  return(gamma_min)
}

