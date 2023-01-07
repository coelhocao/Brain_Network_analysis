#############################################################################################
# This is an implementation of the anogva procedure described by 
# Fujita et al, 2017 DOI: 10.3389/fnins.2017.00066
# 
# Important noteÇ In Fujita 2017, this procedure is applied on a population of graphs representing
# subjectss (fMRI data) that are further grouped. And the bootstrap procedure resamples the graph sample 
# across groups. Here, the procedure is calculated given one single network per group, and the bootstrap 
# procedure resamples the subjects before network construction, then network is constructed, 
# spectral entropy calculated for each network and similarity is calculated by Kullback-Leibler divergence,
# which is a generalization of the Jensen-Shannon divergence when there are more than 2 networks..
#
# Cesar A O Coelho
# cebacio@gmail.com
#
# TODO
# - s20, g2 v
# - s50, g2 (N) v
# - s40, g4 (group) v
# - s200, g4 (group * N)
# - repeat with diff_samples
#   
#

if (!('igraph' %in% installed.packages()[,'Package'])){install.packages("igraph")}; require(igraph)
if (!('statGraph' %in% installed.packages()[,'Package'])){install.packages("statGraph")}; require(statGraph)

entropy_similarity <- function(graph_list, labels = NULL, bandwidth='Silverman', all_divergences = FALSE){
  # Calculates the spectral entropy similarity between all matrices graph_list.
  # Uses Jensen-Shannon divergence (JS) for 2 matrices and Kullback-Leibler divergence (KL) for more.
  # Function adapted from anogva() function in statGraph package.
  # Based on Takahashi et al, 2012. DOI: 10.1371/journal.pone.0049949
  # and on FUjita et al, 2017. DOI: 10.3389/fnins.2017.00066
  # 
  # Input: 
  #        graph_list       List of matrices, each being a network of a given group. Unweighted graphs
  #                         should be binary matrices, weighted graphs, continuous weights. 
  #        labels           Optional. NULL(default). vector with factors. Must have same length as graph_list
  #        bandwidth        String indicating which criterion will be used to choose the bandwidth during
  #                         the spectral density estimation. The available criteria are "Silverman" (default)
  #                         and "Sturges".
  #        all_divergences  Logical. Indicates whether (TRUE,default) to provide the individual
  #                         contribution of each network to the Kullback-Leibler Divergence, or (FALSE) just
  #                         the Kullback-Leibler Divergence(mean of all contributions).
  #
  # Output: 
  #        Kullback-Leibler divergence statistic
  stopifnot(is.logical(all_divergences))
  
  if(is.igraph(graph_list[[1]])){
    graph_list <- matrix_list(graph_list) # Here I can save a function by applying Multi_net_measures(graph_list, FUN = function(i){as.matrix(i[])})
  }
  if(!is.null(labels)){
    stopifnot(length(labels) == length(graph_list))
  }
  else{
    if(is.null(names(graph_list))){
      labels <- seq(graph_list)
    }
    else{
      labels <- c(names(graph_list))
    }
  }
  
  
  f <- statGraph:::nSpectralDensities(graph_list, bandwidth = bandwidth)
  band <- length(f$x)
  media <- f$densities
  colnames(media) <- labels
  mediaAll <- list(y = array(0, band), 
                   x = f$x)
  for(i in 1:band){
    mediaAll$y[i] <- mean(media[i,])
  }
  
  distOrig <- 0
  group_media <- list(x = f$x)
  distOrig <- sapply(labels, function(i){
    group_media$y <- media[,i]
    distOrig <- distOrig + statGraph:::KL(group_media,mediaAll)
    return(distOrig)
  })
  
  names(distOrig) <- labels
  if(isTRUE(all_divergences)){
    distOrig <- c(distOrig, KL = sum(distOrig)/length(labels))
  }
  else{
    distOrig <- sum(distOrig)/length(labels)
  }
  
  return(distOrig)
}

adapted_anogva <- function(df, to_sample = 'group', resample = 10000, bandwidth = 'Silverman', all_divergences = FALSE,
                          pairwise = TRUE, seed=1987, null_distrib = FALSE, ...){
  #
  # ANOGVA DESCRIPTION
  #
  # Input:  
  #     df                 List of dataframes output of mfile_to_dataframe() or similar. One elements of the 
  #                        list is $counts, a NxM dataframe, where N = subjects and M = brain_regions, and 
  #                        the data is a proxy of brain activity. Another element, $labels, has all the grouping variables
  #     to_sample          Vector with names of grouping variables
  #     resample           number of resamplings in the null hypothesis distribution
  #     bandwidth          String indicating which criterion will be used to choose the bandwidth during
  #                        the spectral density estimation. The available criteria are "Silverman" (default)
  #                        and "Sturges".
  #     all_divergences    Logical. Indicates whether (TRUE,default) to provide the individual
  #                        contribution of each network to the Kullback-Leibler Divergence, or (FALSE) just
  #                        the Kullback-Leibler Divergence(mean of all contributions).
  #     pairwise           Logical. Whether (TRUE default) or not calculate KL of all combinations of 
  #                        groups 2x2.
  #     seed               seed for simulations
  #     null_distrib       Logical. Whether or not (FALSE, default) to provide the resampled distribution
  #                        used in the permutation test
  #
  # Output: 
  #

  stopifnot(is.logical(all_divergences) & is.logical(pairwise) & is.logical(null_distrib))

  if(length(to_sample) > 1){
    df$labels$interac <- factor(with(df$labels, interaction(df$labels[,match(to_sample, names(df$labels))])))
    to_sample <- 'interac'
  }
  
  if(length(levels(df$labels[,to_sample])) < 2){
    stop("ERROR! It's a between network design. Requires at least TWO labels necessary to compare")
  }
  
  set.seed(seed)
  resamp_distrib <- sapply(seq(resample), function(i){
    resamp_KL(df, to_sample, bandwidth = bandwidth, pairwise = pairwise, 
            all_divergences = all_divergences, ...)
  })
  
  #empirical networks generation and differences
  emp <- net_list(df, factors = to_sample, ...)
  
  #Extracting measures
  emp_measures <- entropy_similarity(emp, bandwidth = bandwidth, all_divergences = all_divergences)
  
  if(isTRUE(pairwise)){
    comb_label <- combn(names(emp),m=2, FUN = paste0)
    emp_sim <- list()
    emp_comb <- combn(emp,2, simplify = FALSE)
    emp_sim <- lapply(seq(emp_comb), function(i) {
      entropy_similarity(emp_comb[[i]], bandwidth = bandwidth, 
                         all_divergences = F)
    })
    names(emp_sim) <- sapply(1:ncol(comb_label), function(x){ paste(comb_label[,x], collapse = " ~ ")})
    emp_measures <- c(KL = emp_measures, unlist(emp_sim))
  }

  #stat <- p.value(list(t0 = emp_measures, t = data.frame(t(resamp_distrib))))
  stat <- list(statistic = emp_measures, 
               p = sapply(seq_along(emp_measures), function(i)
                 #length(which(resamp_distrib[i,] >= emp_measures[i]))/(nrow(resamp_distrib) +  1)
                 mean(abs(resamp_distrib[i,] > emp_measures[i]))
                 )
               )
  names(stat$p) <- names(stat$statistic)

  if(isTRUE(null_distrib)){
    stat$random_distrib = resamp_distrib
  }
  return(stat)
}

resamp_KL <- function(df, to_sample = 'group', bandwidth = 'Silverman', 
                    all_divergences = FALSE, pairwise = TRUE,...){
  
  #shuffle the grouping labels
  if(length(to_sample) > 1){
    df$labels$interac <- factor(with(df$labels, interaction(df$labels[,match(to_sample, names(df$labels))])))
    to_sample <- 'interac'
  }
  df$labels[[to_sample]] <- sample(df$labels[[to_sample]], replace = F)
  
  #randomized networks
  nh_nets <- net_list(df, factors = to_sample,...)
  
  #Extracting measures
  nh_KL <- entropy_similarity(nh_nets, bandwidth = bandwidth, 
                                     all_divergences = all_divergences)
  
  if(isTRUE(pairwise)){
    comb_label <- combn(names(nh_nets),m=2, FUN = paste0)
    nh_sim <- list()
    nh_comb <- combn(nh_nets,2, simplify = FALSE)
    nh_sim <- lapply(seq(nh_comb), function(i) {
      entropy_similarity(nh_comb[[i]], bandwidth = bandwidth, 
                         all_divergences = F)
    })
    names(nh_sim) <- sapply(1:ncol(comb_label), function(x){ paste(comb_label[,x], collapse = " ~ ")})
    nh_KL <- c(nh_KL, unlist(nh_sim))
  }
  
  return(nh_KL)
}

adapted_anogva_test <- function(subjects = 50, groups = 2, size = 100, runs = 1, resample = 1000, seed = 1987, p_distrib = FALSE, 
                                diff_samples = TRUE, means = c(400,400), sds = c(10,50)){
  #
  # this function tests statistical parameters of adapted_anogva 
  #
  #
  stopifnot(is.numeric(subjects) | is.numeric(groups) | is.numeric(size) |
      is.numeric(runs) | is.numeric(resample) | is.logical(p_distrib) |
      is.logical(diff_samples))
  set.seed(seed)
  seed = as.integer(rnorm(runs, mean = 1000, sd = 100))
  
  simulations <- sapply(seq(runs), function(i){
    df = list()
    df$labels <- data.frame(group = sample(as.factor(rep(seq(groups), subjects/groups)), replace = FALSE))
    #df$counts <- as.data.frame(matrix(rbinom(subjects*size, size = 400, prob = probs[1]), subjects, size))
    df$counts <- as.data.frame(matrix(rnorm(subjects*size, means[1], sds[1]), subjects, size))
    if(isTRUE(diff_samples)){
      df$counts[df$labels$group %in% c(1:round(as.numeric(groups/2))),] <- 
        #as.data.frame(matrix(rbinom(length(df$counts[df$labels$group %in% c(1:round(as.numeric(groups/2))),]), size = 400, prob = probs[2]),
        as.data.frame(matrix(rnorm(length(df$counts[df$labels$group %in% c(1:round(as.numeric(groups/2))),]), means[2], sds[2]),
                             nrow(df$counts[df$labels$group %in% c(1:round(as.numeric(groups/2))),]),
                             ncol(df$counts[df$labels$group %in% c(1:round(as.numeric(groups/2))),])))
    }
    
    sim <- adapted_anogva(df = df, to_sample = 'group', resample = resample, seed = seed[i],
                          annotated = F, pairwise = F, all_divergences = T, null_distrib = FALSE)
    return(as.matrix(sim$p[length(sim$p)]))
  })
  rj <- seq(0,1,0.02)
  rejection_rate <- sapply(rj, function(i){
      mean(abs(simulations < i))
  })
  if(isTRUE(p_distrib)){
    output = list(x = rejection_rate, y = rj, perm = simulations)
  }
  else{
    output = list(x = rejection_rate, y = rj)
  }
  output$params <- c(subjects = subjects, groups = groups, size = size, runs = runs, resample = resample, means, sds)
  return(output)
}

#random_anogva_equal <- adapted_anogva_test(subjects = 40, groups = 4, size = 100, runs = 1000, resample = 1000, seed = 1987, 
#                                            diff_samples = F, means = c(400,400), sds = c(10,10))
#random_anogva_N1_equal <- adapted_anogva_test(subjects = 100, groups = 4, size = 100, runs = 1000, resample = 1000, seed = 1987, 
#                                           diff_samples = F, means = c(400,400), sds = c(10,10))
#random_anogva_N2_equal <- adapted_anogva_test(subjects = 200, groups = 4, size = 100, runs = 1000, resample = 1000, seed = 1987, 
#                                              diff_samples = F, means = c(400,400), sds = c(10,10))


#random_anogva_diff <- adapted_anogva_test(subjects = 40, groups = 4, size = 100, runs = 1000, resample = 1000, seed = 1000, 
#                                          diff_samples = T, means = c(400,400), sds = c(10,20))
#random_anogva_diff <- adapted_anogva_test(subjects = 50, groups = 4, size = 100, runs = 1000, resample = 1000, seed = 1000, 
#                                               diff_samples = T, means = c(400,400), sds = c(10,20))
#random_anogva_N0.5_diff <- adapted_anogva_test(subjects = 60, groups = 4, size = 100, runs = 1000, resample = 1000, seed = 1000, 
#                                           diff_samples = T, means = c(400,400), sds = c(10,20))
#random_anogva_N1_diff <- adapted_anogva_test(subjects = 100, groups = 4, size = 100, runs = 1000, resample = 1000, seed = 1987, 
#                                              diff_samples = T, means = c(400,400), sds = c(10,20))
#random_anogva_N2_diff <- adapted_anogva_test(subjects = 200, groups = 4, size = 100, runs = 1000, resample = 1000, seed = 1987, 
#                                              diff_samples = T, means = c(400,400), sds = c(10,20))
#random_anogva_size_N_diff <- adapted_anogva_test(subjects = 100, groups = 4, size = 300, runs = 1000, resample = 1000, seed = 1987, 
#                                                  diff_samples = T, means = c(400,400), sds = c(10,20))


        