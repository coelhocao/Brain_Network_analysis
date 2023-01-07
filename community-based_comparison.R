source('/home/cesar/Dropbox/R/Network sufficiency/Modular_codes/community-based_nodal_metrics.R');
source('/home/cesar/Dropbox/R/Network sufficiency/Modular_codes/multi_scale_comm.R');

cbm_permutation <- function(df, memb, to_sample = c('signal', 'group') , seed=1987, resample=10000,
                                metrics = c('within_community_strength', 'diversity','participation', "gateway_bet","gateway_str"), 
                                normalized = T, null_distrib = FALSE, sep_values = T, ...){
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
    stop("ERROR! Requires at least TWO levels in to_sample factor(s) to compare")
  }
  
  set.seed(seed)
  resamp_distrib <- sapply(seq(resample), function(i){
    cbm_resamp(df, memb = memb, to_sample = to_sample,
              metrics = metrics, normalized = normalized, sep_values = sep_values, ...)
  }, simplify = 'array')
  
  #empirical networks generation and differences
  emp <- net_list(df, factors = to_sample,...)
  
  #Extracting measures
  emp_measures <- multi_net_measures(emp, sec_param = memb, FUN = community_based_metrics, list_levels = 1, 
                                     add_param = 'comm', metrics = metrics, normalized = normalized, 
                                     different_sizes = T, sep_values = sep_values)
  emp_measures <- lapply(emp_measures, as.data.frame)
  
  diff_label <- combn(names(emp_measures), m=2, FUN = paste0)
  emp_diffs <- list()
  comb <- combn(emp_measures,2, simplify = FALSE)
  emp_diffs <- sapply(comb, function(x){ 
    as.data.frame(x[[1]] - x[[2]])
  }, simplify = FALSE)
  names(emp_diffs) <- sapply(seq(ncol(diff_label)), function(x){ 
    paste(diff_label[,x], collapse = " ~ ")
  })
  
  emp_diffs <- cbind.data.frame(emp_diffs)
  
  #p_values
  p <- matrix(0, nrow = nrow(emp_diffs), ncol = ncol(emp_diffs))
  
  for(i in seq(nrow(emp_diffs))){
    for(j in seq(ncol(emp_diffs))){
      p[i,j] <- mean(abs(resamp_distrib[i,j,]) > abs(emp_diffs[i,j]))
    }
  }
  colnames(p) <- colnames(emp_diffs)
  rownames(p) <- rownames(emp_diffs)
  p <- as.data.frame(p)
  
  output <- list(statistic = emp_diffs, pvalue = p)
  if(isTRUE(null_distrib)){
    output$null_distrib = resamp_distrib
  }
  return(output)
}

cbm_resamp <- function(df, memb, to_sample = 'group', sep_values = sep_values,
                       metrics = c('within_community_strength', 'diversity','participation', "gateway_bet","gateway_str"), 
                       normalized = T, ...){
  
  #shuffle the grouping labels
  if(length(to_sample) > 1){
    df$labels$interac <- factor(with(df$labels, interaction(df$labels[,match(to_sample, names(df$labels))])))
    to_sample <- 'interac'
  }
  df$labels[[to_sample]] <- sample(df$labels[[to_sample]], replace = F)
  
  #randomized networks
  nets <- net_list(df, factors = to_sample, ...)
  
  #
  same <- unlist(sapply(seq(nets), function(i){
    c(length(memb[[i]]$names[!(memb[[i]]$names %in% V(nets[[i]])$name)]) > 0,
      length(V(nets[[i]])$name[!(V(nets[[i]])$name %in% memb[[i]]$names)]) > 0)
  }))
  if(any(same != 0)){
    nets <- match_net_to_comm(nets, memb)
  }
    
  # Nodal measures for the networks
  nets_measures <- multi_net_measures(nets,  sec_param = memb, FUN = community_based_metrics, list_levels = 1, 
                                      add_param = 'comm', metrics = metrics, normalized = normalized, 
                                      different_sizes = T, sep_values = sep_values)
  nets_measures <- lapply(nets_measures, function(i){
    as.data.frame(i)
  })
  
  diff_label <- combn(names(nets_measures), m=2, FUN = paste0)
  diffs <- list()
  comb <- combn(nets_measures,2, simplify = FALSE)
  diffs <- sapply(comb, function(x){ 
    as.data.frame(x[[1]] - x[[2]])
  }, simplify = FALSE)
  names(diffs) <- sapply(seq(ncol(diff_label)), function(x){ 
    paste(diff_label[,x], collapse = " ~ ")
  })
  
  diffs <- cbind.data.frame(diffs)
  
  return(as.matrix(diffs))
}

match_net_to_comm <- function(G, memb){
  
  stopifnot(is.list(G) & is.list(memb))
  stopifnot(unlist(lapply(G, is.igraph)) | 'communities' %in% attr(memb, 'class'))
  
  mnet <- lapply(seq(G), function(i){
    if(length(memb[[i]]$names[!(memb[[i]]$names %in% V(G[[i]])$name)]) > 0){
      nv = length(memb[[i]]$names[!(memb[[i]]$names %in% V(G[[i]])$name)])
      G[[i]] <- add_vertices(G[[i]],nv, attr =  list(name = c(memb[[i]]$names[!(memb[[i]]$names %in% V(G[[i]])$name)])))
    }
    if(length(V(G[[i]])$name[!(V(G[[i]])$name %in% memb[[i]]$names)]) > 0){
      G[[i]] <- delete.vertices(G[[i]], V(G[[i]])$name[!(V(G[[i]])$name %in% memb[[i]]$names)])
    }
    return(G[[i]])
  })
  names(mnet) <- names(G)
  return(mnet)
}
