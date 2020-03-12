#############################################################################################
# These codes have been costumly designed for brain network analysis of 1-2 time point data
# generally acquired by immunolabelling of activity-dependent proteins (i.e. c-fos, Arc)
#
# Many functions here have been largely based on some of the codes shared by Justin Kenney,
# which can be found on this link https://github.com/jkenney9a/Networks
#
# NOTE: It is intended mostly for undirected network data
# Most funcitons won't work in structural brain data or networks which support directed graphs 
# or different upper/lower matrix triangles.
#
# NOTE: Designed mostly for comparisons between single-time point graphs or two-time points graphs
# won't work on multiple time points graphs
#
# Cesar A O Coelho
# cebacio@gmail.com
#
#########################################################################################################################################
#TODO
# - calculate consensus or change between two network states (community assignement agreement)
# - Calculate louvain algorithm allowing signed networks and resolution parameter
# - Calculate Multi-resolution consensus function allowing for different algorithms
#    - varia gamma -> pega comunidades -> calcula consensus matrix 
#    -> aleatoriza classificacao de comunidades em cada gamma -> calcula random consensus matrix
#    -> EMPconsensus - RANDOMconsensus -> reclustering 
# - Statistical significance testing. Repeat the above procedure in rAndomized networks N times
#    -> calculate Qz = (Qemp - Qrandom^2)/sqrt(Qrandom^2(1 - Qrandom^2))
#    -> calculate
# - implement gen_louvain from Lucas et al 2018 ?? TOO DIFFICULT FOR NOW

if (!('dplyr' %in% installed.packages()[,'Package'])){install.packages("dplyr")}; require(dplyr)
if (!('Hmisc' %in% installed.packages()[,'Package'])){install.packages("Hmisc")}; require(Hmisc) #better functions for correlation matrices and p-values
if (!('lattice' %in% installed.packages()[,'Package'])){install.packages("lattice")}; require(lattice) #for generating graphs, colors and matrices
if (!('data.table' %in% installed.packages()[,'Package'])){install.packages("data.table")}; require(data.table) #for rbindlist
if (!('ggcorrplot' %in% installed.packages()[,'Package'])){install.packages("ggcorrplot")}; require(ggcorrplot) #for plotting correlation heatmaps
if (!('data.tree' %in% installed.packages()[,'Package'])){install.packages("data.tree")}; require(data.tree) #for hierarchical organization of anatomical regions
if (!('resolution' %in% installed.packages()[,'Package'])){devtools::install_github("analyxcompany/resolution")}; require(resolution) #for community detection varying resolution parameter
if (!('brainGraph' %in% installed.packages()[,'Package'])){install.packages("braingraph")}; require(brainGraph) #participation and Gateway and Diversity coefficients
if (!('igraph' %in% installed.packages()[,'Package'])){install.packages("igraph")}; require(igraph) #package for network generation and measures
if (!('fda' %in% installed.packages()[,'Package'])){install.packages("fda")}; require(fda) #package for network generation and measures
if (!('psych' %in% installed.packages()[,'Package'])){install.packages("psych")}; require(psych) #package for matrix additions


#if (!('boot' %in% installed.packages()[,'Package'])){install.packages("boot")}; require(boot)
#if (!('multcomp' %in% installed.packages()[,'Package'])){install.packages("multcomp")}; require(multcomp)
#if (!('heplots' %in% installed.packages()[,'Package'])){install.packages("heplots")}; require(heplots)
#if (!('effsize' %in% installed.packages()[,'Package'])){install.packages("effsize")}; require(effsize)
#if (!('statGraph' %in% installed.packages()[,'Package'])){install.packages("statGraph")}; require(statGraph) #For Jensen-Shannon divergence between two graphs
#if (!('cocoframer' %in% installed.packages()[, 'Package'])){devtools::install_github("AllenInstitute/cocoframer")}; require (cocoframer) #3D visualization of brain (can it do networks?)

source('C:/Users/CAOC/Dropbox/R/Network sufficiency/Modular_codes/activity_comparisons.R');
source('C:/Users/CAOC/Dropbox/R/Network sufficiency/Modular_codes/signed_louvain.R');
source('C:/Users/CAOC/Dropbox/R/Network sufficiency/Modular_codes/modularity_maximization.R');
source('C:/Users/CAOC/Dropbox/R/Network sufficiency/Modular_codes/community-based_nodal_metrics.R');

#########################################################################################################################################
##########################                  COMMUNITY-BASED MEASURES AND ANALYSIS                 #######################################
#########################################################################################################################################




multi_scale_comm <- function(G, FUN = signed_louvain, gamma = seq(-0.5, 1, 0.05), param_name = 'gamma', sec_param_name = NULL, ...){
  #
  # This function repeats the community detection algorithm given in FUN across a number of resolution values.
  #
  # INPUT: 
  #       G,        igraph object or dataframe structured as detailed in the 
  #       FUN,      Community detection algorithm. It accepts four functions for now
  #                 signed_louvain(default, adapted from Brain connectivity toolbox for matLab), 
  #                 cluster_resolution (resolution pckg), Louvain algorithm, but unoptimized and does not accept signed graphs)
  #                 spinglass.community (igraph pckg)
  #                 modularity maximization (adapted from Brain connectivity toolbox for matLab)
  #       gamma,    parameter resolution (gamma in all functions but cluster_resolution, which is called t)
  #
  # OUTPUT: List(
  #         comm,   community affiliation vactors
  #         gamma,  resolution parameter value for each community affiliation
  #
  
  
  comm <- lapply(gamma, function(i){
    if(is.null(sec_param_name)){
      params <- list(G, i, ...)
      names(params)[2] <- param_name
      FUN(G, i, ...)
    }
    else{
      params <- list(G, i, ...)
      names(params)[2] <- param_name
      names(params)[3] <- sec_param_name
      FUN(G, i, i,...)
    }
  })

  output <- list()
  output$comm <- comm
  output$gamma <- gamma
  
  return(output)
}

agreement <- function(ci){
  
  
  n <- ncol(ci);
  r <- nrow(ci);

  # Can be done with fastDummies package as below
  #ind <- fastDummies::dummy_cols(matrix(as.factor(ci), ncol = n));
  #ind <- ind[,(n+1):ncol(ind)]
  #Or
  ind <- numeric()
  ind <- as.data.frame(sapply(seq(ncol(ci)), function(x){
    sapply(seq(max(ci[,x])), function(y){
      ind <- as.matrix(as.numeric(ci[,x] == y), nrow = r)
    })
  }, simplify = F))
  
  D <- ind %*% Conj(t(ind));

  
}


community_distance <- function(comm, method = 'vi'){
  #
  # INPUT: igraph community object with multiple repetitions of community classification (e.g. louvain, resolution etc)
  # OUTPUT: Mean of all the pairwise information variation (VI). 
  
  n <- nrow(comm$memberships)
  vi <- matrix(0, nrow = n, ncol = (n-1))
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      vi[j,i] <- compare(comm$memberships[i,], comm$memberships[j,], method = method)
      
    }
  }
  return (vi = mean(vi[lower.tri(vi)]))
  
}

multi_scale_comm_stability <- function(comm, stat = NULL, method = 'vi'){
  #
  # INPUT: list of communities from multi_scale_comm() varying resolution parameter
  # OUTPUT: data.frame with mean variation of information and nnumber of communities at each parameter
  #
  m <- sapply(comm$comm, function(i){ community_distance(i, method = method) })
  
  cn <- sapply(comm$comm, function(i){ length(i) })
  
  s <- data.frame(gamma = comm$gamma, community_distance = m, communities = cn)
  
  if(!is.null(stat)) {  s <- data.frame(s,t(stat))  }
  
  return(s)
}


comm_resamp <- function(G, comm, isolated = F){
  #
  # INPUT: an igraph object and a community object
  # 
  # OUTPUT: maximazed modularity of the randomized community assignments
  #
  
  # Delete isolated nodes
  if(isolated == F){ G <- del_isolated(G) }
  
  rc <- comm
  
  # shuffles the community labels in all repetitions of the community assignments
  rc$memberships <- sapply(seq(nrow(rc$memberships)), function(x){
    rc$memberships[x,] <- sample(rc$memberships[x,], replace = T)
  })
  # transpose the df
  rc$memberships <- sapply(seq(nrow(rc$memberships)), function(x){
    rc$memberships[x,] <- t(rc$memberships[x,])
  })
  # recalculates the modularity in each repetition
  rc$modularity <- sapply(seq(nrow(rc$memberships)), function(x){
    rc$modularity[x] <- modularity(G, rc$memberships[x,])
  })
  
  #elects the max modularity
  rc$modularity <- modularity(rc)
  
  return(rc)  
}

comm_permutation_test <- function(G, comm, isolated = F, seed = 1987, R = 1000){
  #
  # INPUT: an igraph object and a community object
  #
  # OUTPUT: max modularity, std error and p_value against distribution of community assignments
  
  set.seed(seed)
  
  random_comms <- lapply(comm$comm, function(i){
    replicate(R, comm_resamp(G, i, isolated = isolated))
  })
  
  #random_comms is supposed to be a nested t lists of R lists. but
  #Sometimes, for some reason, one of the t lists will be a list of 6*R lists of single elements and not R lists with community attribute
  # this piece solves it but the R* lits lose the communities attribute.
  random_comms <- lapply(random_comms, function(i){
    if(length(i)  == R){
      i = i
    }
    else { 
      lapply((seq(R)*6-5), function(j){
        list(membership = i[[j]], memberships = i[[j+1]], modularity = i[[j+2]], 
             names = i[[j+3]], vcount = i[[j+4]], algorithm = i[[j+5]])
      })
    }
  })
  
  
  modularity_distrib <- sapply(random_comms, function(i){
    sapply(i, function(j){
      j$modularity
    })
  })
  
  modularity_distrib <- list(t = modularity_distrib,
                             t0 = sapply(seq_along(comm$comm), function(i){ modularity(comm$comm[[i]]) }))
  
  results <- rbind(p.value(modularity_distrib),
                   lower_ci = sapply(seq_along(comm$comm), function(i){ quantile(modularity_distrib$t[,i], probs = 0.025) }),
                   upper_ci = sapply(seq_along(comm$comm), function(i){ quantile(modularity_distrib$t[,i], probs = 0.975) }))
  colnames(results) <- comm$resolution
  
  return(results)
}

from_graph_to_comm_analysis <- function(G, comm = NULL, isolated = F, resolution = c(seq(0.5, 10,0.5)), rep = 5, RandomOrder = T, method = 'vi', 
                                        bootstrap = T, seed = 1987, R = 10){
  #
  # INPUT: igraph object or list of igraph objects. Optionally, the community structure
  #
  # OUTPUT: communities structures, the Louvain resolution used in each, and a df with their stability data
  
  if(!is.null(comm)){
    
    if(bootstrap == T){
      comm$comm <- comm$communities
      comm_stat <- comm_permutation_test(G, comm, isolated = isolated, R = R, seed = seed)  
      comm_stability <- multi_scale_comm_stability(comm, comm_stat, method = method)
    }
    else {
      comm_stability <- multi_scale_comm_stability(comm)
    }
    comms_data <- comm
    comms_data$stability <- comm_stability
    
  }
  else{
    comm <- multi_scale_comm(G, isolated = isolated, resolution = resolution, rep = rep, RandomOrder = RandomOrder)
    if(bootstrap == T){
      comm_stat <- comm_permutation_test(G, comm, isolated = isolated, R = R, seed = seed)
      comm_stability <- multi_scale_comm_stability(comm, comm_stat, method = method)
    }
    else {
      comm_stability <- multi_scale_comm_stability(comm)
    }
    
    comms_data <- list()
    comms_data$communities <- comm$comm
    comms_data$resolution <- comm$resolution
    comms_data$stability <- comm_stability
  }
  return(comms_data)
}


get_stability <- function(comms){
  #
  # INPUT: datalist from  from_graph_to_comm_analysis() function
  #
  # OUTPUT: dataframe to be used in plot_community_stability() function
  
  output <- list()
  output <- lapply(comms, function(i){
    lapply(i, function(j){
      ouput <- j$stability
    })
  })
  
  data <- melt(output, id.vars = names(output[[1]][[1]]))
  names(data) <- c(names(output[[1]][[1]]), 'threshold', 'signal')
  
  return(data)
}

consensus_coefficient <- function(comm){
  #
  # INPUT: community object from igraph which provides multiple repetitions of the community assignment
  # ex. cluster_louvain, walktrap, cluster_resolution, etc...
  #
  # OUPUT: D matrix with the proportion each node falls into the same 
  
  n <- ncol(comm$memberships)
  D <- matrix(0,ncol = ncol(comm$memberships), nrow = ncol(comm$memberships))
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      D[j,i] <- mean(comm$memberships[,j] == comm$memberships[,i])
    }
  }
  D <- as.matrix(D)
  D[upper.tri(D)] <- t(D)[upper.tri(D)]
  colnames(D) <- comm$names
  rownames(D) <- comm$names
  return(D)
}

#########################################################################################################################################
#########################################################################################################################################
###########           FUNCTIONS THAT NEED SOME WORK


consensus_to_convergence <- function(G, FUN = modularity_maximization, tau = 0.5, max = 10, ...){
  # NOT WORKING FOR SOME REASON
  # INPUT: an igraph object, a clustering function, a consensus threshold
  #
  # OUTPUT: a community object outputed from a repeated consensus trimming until convergence (total consensus)
  # NOTE: This function is based on the paper: 
  # "Lancichinetti & Fortunato, 2011 - Consensus clustering in complex networks"
  
  #community assignment
  
  comm <- FUN(G, ...)
    
  
  # consensus matrix of community co-assignment
  D <- consensus_coefficient(comm)
  
  # threshold to consensus
  D[D < tau] <- 0
  
  # connect unconnected nodes
  if(0 %in% colSums(D)){
    D[colSums(D) == 0, max(colSums(D))] <- 0.01
  }
  
  Dash <- D
  n = 0
  s = 0
  while(n == 0 & s < max){
    
    Dc <- FUN(graph.adjacency(as.matrix(Dash), mode = 'undirected', weighted = T), ...)
    #Dc <- cluster_resolution(Dash, RandomOrder = T, rep = 100, t=3)
    Dash <- consensus_coefficient(Dc)
    
    Dash[Dash < tau] <- 0    
    
    # connect unconnected nodes
    if(0 %in% colSums(Dash)){
      Dash[colSums(Dash) == 0, max(colSums(Dash))] <- 0.01
    }
    
    if(length(Dash[Dash > 0 & Dash < 1]) > 0){ 
      n = 1  
    }
    s <- s+1
  }
  print(s)
  return(Dc)
}

consensus_prob <- function(G, FUN= cluster_resolution, isolated = F, seed = 1987, R = 10, mode = 'undirected', weighted = T, ...){
  #
  # This is an alternative to find your communities. it is still based on consensus coefficient
  # INPUT: consensus matrix and community structure object
  # 
  # OUTPUT: consensus of community co-assignament matrix at random
  #
  # NOTE: the rationale here is that presented in Betzel et al, 2019 - Temporal fluctuations in the brain's modular architecture during movie-watching
  # But for some reason subtracting this matrix from D gives out negative matrices.
  #
  
  if(isolated == T){ G <- del_isolated(G)}
  emp_comm <- FUN(G, ...)
  #emp_comm <- cluster_resolution(G, t = 3, rep = 100, RandomOrder = T)
  D <- consensus_coefficient(emp_comm)
  set.seed(seed)
  rand_comm_distrib <- replicate(R, comm_resamp(G, emp_comm, isolated = isolated))
  
  D_rand_distrib <- lapply(rand_comm_distrib, function(x){
    consensus_coefficient(x)
  })
  
  # turns into array
  D_rand_distrib <- array(unlist(D_rand_distrib), dim = c(nrow(D_rand_distrib[[1]]), ncol(D_rand_distrib[[1]]), length(D_rand_distrib)))
  
  # probability D = random
  p_D <- matrix(0, nrow = nrow(D), ncol = ncol(D))
  n <- ncol(D)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      p_D[j,i] <- mean(D_rand_distrib[j,i,])##this line is the one I have to figure out
    }
  }
  
  #makes p_D simmetric matrix
  p_D <- as.matrix(p_D)
  p_D[upper.tri(p_D)] <- t(p_D)[upper.tri(p_D)]
  colnames(p_D) <- colnames(D)
  rownames(p_D) <- rownames(D)
  
  # diff between empirical and random
  D <- D - p_D
  D[D < 0] <- 0
  
  #community structure
  D <- graph.adjacency(as.matrix(D), mode = mode, weighted = weighted)
  D <- FUN(D, ...)
  return(D)
}
