
resamp_nodals <- function(df, to_sample = 'group', efficiency_loss = FALSE, local_eff = FALSE, 
                          metrics = c('signed_strength'), 
                          net_list_args = c(...), multi_net_measures_args = c(...)){
  
  #shuffle the grouping labels
  if(length(to_sample) > 1){
    df$labels$interac <- factor(with(df$labels, interaction(df$labels[,match(to_sample, names(df$labels))])))
    to_sample <- 'interac'
  }
  df$labels[[to_sample]] <- sample(df$labels[[to_sample]], replace = F)
  
  #randomized networks
  rnets <- do.call(net_list, list(df, factors = to_sample, unlist(net_list_args)))
  
  # Nodal measures for the networks
  nets_measures <- do.call(multi_net_measures,
                           list(rnets, FUN = Nodal_measures, list_levels = 1, 
                                      efficiency_loss = efficiency_loss, local_eff=local_eff,
                                      metrics = metrics, different_sizes = TRUE, out = 'data.frame', 
                                      unlist(multi_net_measures_args)))
  if(is.list(nets_measures[[1]])){
    nets_measures <- multi_net_measures(nets_measures, FUN = as.data.frame, list_levels = 1)
  }
  
  diff_label <- combn(names(nets_measures), m=2, FUN = paste0)
  diffs <- list()
  comb <- combn(nets_measures,2, simplify = FALSE)
  diffs <- sapply(comb, function(x){ 
    as.data.frame(x[[2]] - x[[1]])
  }, simplify = FALSE)
  names(diffs) <- sapply(seq(ncol(diff_label)), function(x){ 
    paste(diff_label[2:1,x], collapse = " ~ ")
  })
  diffs <- cbind.data.frame(diffs)
  return(as.matrix(diffs))
}

nodal_permutation_test <- function(df, to_sample = 'group', seed = 1987, resample = 10000, efficiency_loss= FALSE, 
                                   local_eff=FALSE, metrics = c('signed_strength'), 
                                   net_list_args = list(),multi_net_measures_args = list()){
  #
  # Input: dataframe, grouping labels and resampling number
  #
  # resample labels, regenerate graphs, calculate network measures, group diffs among them
  # generate empirical networks, calculate measures and group differences
  # calculate the proportion of resampled diffs > empirical differences (p-value)
  #
  # NOTE: The parameters to be put are both for the empirical and resampled graphs
  # NOTE2: The function accepts multiple groups, but Pvalue is not corrected for multiple comparisons
  #
  # Output: List with empirial differences, 2-by-2 and permutation p-values
  #
  
  if(length(to_sample) > 1){
    df$labels$interac <- factor(with(df$labels, interaction(df$labels[,match(to_sample, names(df$labels))])))
    to_sample <- 'interac'
  }
  if(length(levels(df$labels[,to_sample])) < 2){
    stop("ERROR! It's a between network design. Requires at least TWO labels necessary to compare")
  }
  
  set.seed(seed)
  resamp_distrib <- sapply(seq(resample), function(i){
    do.call(resamp_nodals, list(df, to_sample, efficiency_loss= efficiency_loss, local_eff=local_eff,
                  metrics = metrics, multi_net_measures_args = multi_net_measures_args,
                  net_list_args = net_list_args))
  }, simplify = 'array')
  
  #empirical networks generation and differences
  nets <- do.call(net_list, list(df, factors = to_sample, unlist(net_list_args)))
  
  # Nodal measures for the networks
  emp_measures <- do.call(multi_net_measures,
                          list(nets, FUN = Nodal_measures, efficiency_loss = efficiency_loss, 
                                     local_eff=local_eff, metrics = metrics, list_levels = 1, 
                                     different_sizes = TRUE, unlist(multi_net_measures_args)))
  if(is.list(emp_measures[[1]])){
    emp_measures <- multi_net_measures(emp_measures, FUN = as.data.frame, list_levels = 1)
  }
  
  diff_label <- combn(names(emp_measures), m=2, FUN = paste0)
  diffs <- list()
  comb <- combn(emp_measures,2, simplify = FALSE)
  diffs <- sapply(seq(comb), function(x){ 
    as.data.frame(comb[[x]][[2]] - comb[[x]][[1]])
  }, simplify = FALSE)
  names(diffs) <- sapply(seq(ncol(diff_label)), function(x){ 
    paste(diff_label[2:1,x], collapse = " ~ ")
  })
  diffs <- cbind.data.frame(diffs)
  
  #Calculating the p-values 
  ps <- matrix(0, nrow(diffs), ncol(diffs))
  for(j in seq(ncol(diffs))){
    for(i in seq(nrow(diffs))){
      ps[i,j] <- mean(abs(resamp_distrib[i,j,]) > abs(diffs[i,j]))
    }
  }
  rownames(ps) <- rownames(diffs)
  colnames(ps) <- colnames(diffs)
  
  return(list(statistic = diffs, pvalue = ps))
}
