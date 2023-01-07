source('/home/cesar/Dropbox/R/Network sufficiency/Modular_codes/signed_louvain.R'); # for optimizer
source('/home/cesar/Dropbox/R/Network sufficiency/Modular_codes/modularity_functions.R'); # for modularity matrix func
source('/home/cesar/Dropbox/R/Network sufficiency/Modular_codes/gamma_range.R'); # to calculate g_min

if (!('doParallel' %in% installed.packages()[,'Package'])){install.packages("doParallel")}; require(doParallel);
if (!('pracma' %in% installed.packages()[,'Package'])){install.packages("pracma")}; require(pracma);
if (!('Hmisc' %in% installed.packages()[,'Package'])){install.packages("Hmisc")}; require(Hmisc);



eventSamples <- function(g, n = 1000, mod_fun = modularity_matrix, optimizer = signed_louvain, 
                         gmin_samples = 100, extrapolFUN = approxExtrap, seed = 1987,
                         optim_Args = list(mod = 'neg_asym', class = F)){
# eventSamples Compute multiresolution ensemble using event sampling.
#
# Syntax
#__________________________________________________________________________
#
#   [S,gammas]=eventSamples(A,n)
#
#   [S,gammas]=eventSamples(__,Name,Value)
#
#
# Description
#__________________________________________________________________________
#
#   [S,gammas]=eventSamples(A,n) computes event sampling ensemble with 'n'
#       partitions.
#
#   [S,gammas]=eventSamples(__,Name,Value) additionally customizes the
#       behavior of the function by e.g. using a different algorithm to
#       optimize modularity or using a different modularity-like quality
#       function.
#
#
# Input Arguments
#__________________________________________________________________________
#
#   A -- Adjacency matrix of the network
#
#   n -- Number of partitions to be generated
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
#                   the function assumes that 'B=(A+A')/2-gamma*P' for some matrix
#                   'P'.
#                   @modularity (default) | function handle
#
#   'GammaMinSamples' -- Number of partitions to sample at each iteration
#                        when estimating 'gamma_min' (see 'gammaRange' for
#                        more details)
#                        10 (default)| scalar
#
#
# Output Arguments
#__________________________________________________________________________
#
#   S -- Ensemble of partitions. This is a matrix where each column
#        corresponds to a partition for a given value of 'gamma'.
#
#   gammas -- Values of 'gamma' corresponding to each partition in 'S'.
#
# See Also hierarchicalConsensus, exponentialSamples

# Version: 1.1.1
# Date: Thu  8 Mar 2018 15:34:46 CET
# Author: Lucas Jeub
# Email: ljeub@iu.edu

  stopifnot(any(is.igraph(g) | 
            is.matrix(g)) |
            is.function(optimizer) |
            is.function(mod_fun) |
            is.function(extrapolFUN)
              )
  
  # gets adjacency matrix from igraph object
  if(is.igraph(g)){
    A <- as.matrix(g[])
  }
  else{
    A <- g
  }
  
  gamma_min <- gammaRange(A, mod_fun, optimizer, samples = gmin_samples, seed = seed);
  gamma_min <- gamma_min[1]
  
  P <- mod_fun(A, 1); 
  A <- (A + t(A))/2;
  P <- A - P;
  
  for(i in seq(max(dim(A)))){
    A[i,i] <- 0;
    P[i,i] <- 0;
  }
  
  # get discrete events where interactions change sign
  gamma_et <- div_0(A,P);
  g_sample <- unique(sort(gamma_et))
  ind <- match(unique(gamma_et), g_sample)
  
  PS <- sum(P);
  AS <- sum(A);
  Pp <- rep(0, length(g_sample));
  Ap <- rep(0, length(g_sample));
  
  for(i in seq(length(g_sample))){
    Pp[i] <- sum(P[ind >= i]);
    Ap[i] <- sum(A[ind >= i]);
  }
  b_sample <- (g_sample * (PS - Pp) - (AS - Ap))/(g_sample * (PS - 2*Pp) + 2*Ap - AS);
  
  #approxExtrap
  b_min <- (gamma_min * approxExtrap(g_sample, PS - Pp, xout = gamma_min)$y - approxExtrap(g_sample,AS - Ap, xout = gamma_min)$y)/
    (gamma_min * approxExtrap(g_sample, PS - 2*Pp, xout = gamma_min)$y + approxExtrap(g_sample, 2*Ap - AS, xout = gamma_min)$y);
  
  
  #predict(loess())
  #b_min <- (gamma_min * predict(loess((PS - Pp) ~ g_sample), data.frame(x = gamma_min)) - 
  #            predict(loess((AS - Ap) ~ g_sample), data.frame(x = gamma_min)))/
  #  (gamma_min * predict(loess((PS - 2*Pp) ~ g_sample), data.frame(x = gamma_min)) + 
  #     predict(loess((2*Ap - AS) ~ g_sample), data.frame(x = gamma_min)));
  #b_min <- min(b_min[b_min >= 0])
  
  #Sometimes, interp1 which doesn't do extrapolation, gives out an error insteat of NA values outside the interpolated interval
  #b_min <- (gamma_min * interp1(g_sample, PS - Pp, gamma_min, 'nearest') - interp1(g_sample,AS - Ap, gamma_min, 'nearest'))/
  #  (gamma_min * interp1(g_sample, PS - 2*Pp, gamma_min, 'nearest') + interp1(g_sample, 2*Ap - AS, gamma_min, 'nearest'));
  
  b_sample1 <- unique(sort(b_sample))
  b_red <- match(b_sample1, unique(b_sample))
  b_sample <- b_sample1
  
  g_sample <- g_sample[b_red];
  Pp <- Pp[b_red];
  Ap <- Ap[b_red];
  
  # Avoid outputting NaN for largest value of gamma due to numerical error
  #and 'next' ('nearest in R implementation') interpolant not supporting extrapolation
  if(b_sample[length(b_sample)] < 1){
    b_sample <- c(b_sample, 1);
    Pp <- c(Pp, Pp[length(Pp)]);
    Ap <- c(Ap, Ap[length(Ap)]);
    g_sample <- c(g_sample, g_sample[length(g_sample)]);
  }
  
  # interpolators/extrapolators to handle events
  # In R, approxExtrap and spline are the closest functions
  # to matlab's gridInterpolant() that perform extrapolation
  b <- pracma::linspace(b_min, 1, n);  # Need this to calculate xout
  
  Pminus <- extrapolFUN (b_sample, PS - Pp, xout = b)$y
  Pplus  <- extrapolFUN (b_sample, Pp     , xout = b)$y
  Aminus <- extrapolFUN (b_sample, AS - Ap, xout = b)$y
  Aplus  <- extrapolFUN (b_sample, Ap     , xout = b)$y
  
  gammas <- (Aminus + b*(Aplus - Aminus)) / ((1-b) * Pminus + b * Pplus)
  
  ##############################################################
  #setup parallel backend to use many processors
  cores=detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  # compute partitions
  Ci <- matrix(0, nrow(as.matrix(g[])), n)
  for(i in seq(n)){
    Ci[,i] <- do.call(optimizer, c(list(A, gammas[i]), optim_Args))
  }
  #stop cluster
  stopCluster(cl)
  ##############################################################
  rownames(Ci) <- rownames(A)
  
  return(list(Ci = Ci, gammas = gammas, b = b, b_events = b_sample, g_events = g_sample)) 
}
