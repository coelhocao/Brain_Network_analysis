source('/home/cesar/Dropbox/R/Network sufficiency/Modular_codes/community-based_nodal_metrics.R');
source('/home/cesar/Dropbox/R/Network sufficiency/Modular_codes/multi_scale_comm.R');

cbm_permutation2 <- function(df, memb, to_sample = c('signal', 'group') , seed=1987, resample=10000,
                                metrics = c('within_community_strength', 'diversity','participation', "gateway_bet","gateway_str"), 
                                normalized = T, null_distrib = FALSE, sep_values = T, base_comm = 1, ...){
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
  
  # Get Null Model distribution
  set.seed(seed)
  resamp_distrib <- sapply(seq(resample), function(i){
    cbm_comp(df, memb = memb, to_sample = to_sample, resamp = T, base_comm = base_comm,
              metrics = metrics, normalized = normalized, sep_values = sep_values, ...)
  }, simplify = 'array')
  
  # Get empirical Statistics
  emp_metrics <- cbm_comp(df, memb = memb, to_sample = to_sample, resamp = F, base_comm = base_comm,
                          metrics = metrics, normalized = normalized, sep_values = sep_values, ...)
  
  #p_values
  p <- matrix(0, nrow = nrow(emp_metrics), ncol = ncol(emp_metrics))
  
  for(i in seq(nrow(emp_metrics))){
    for(j in seq(ncol(emp_metrics))){
      p[i,j] <- mean(abs(resamp_distrib[i,j,]) > abs(emp_metrics[i,j]))
    }
  }
  colnames(p) <- colnames(emp_metrics)
  rownames(p) <- rownames(emp_metrics)
  p <- as.data.frame(p)
  emp_metrics <- as.data.frame(emp_metrics)
  
  output <- list(statistic = emp_metrics, pvalue = p)
  if(isTRUE(null_distrib)){
    output$null_distrib = resamp_distrib
  }
  return(output)
}

cbm_comp <- function(df, memb, to_sample = 'group', sep_values = sep_values,
                       metrics = c('within_community_strength', 'diversity','participation', "gateway_bet","gateway_str"), 
                       normalized = T, resamp = F, base_comm = 1, ...){
  stopifnot(is.logical(resamp))
  stopifnot(base_comm == 1 | base_comm == 2)
  
  #shuffle the grouping labels
  if(length(to_sample) > 1){
    df$labels$interac <- factor(with(df$labels, interaction(df$labels[,match(to_sample, names(df$labels))])))
    to_sample <- 'interac'
  }
  if(isTRUE(resamp)){
    df$labels[[to_sample]] <- sample(df$labels[[to_sample]], replace = F)
  }
  
  #randomized networks
  nets <- net_list(df, factors = to_sample, ...)
  
  #
  same <- unlist(sapply(seq(nets), function(i){
    c(length(memb[[i]]$names[!(memb[[i]]$names %in% V(nets[[i]])$name)]) > 0,
      length(V(nets[[i]])$name[!(V(nets[[i]])$name %in% memb[[i]]$names)]) > 0)
  }))
  if(any(isTRUE(same))){
    nets <- match_net_to_comm(nets, memb)
  }
  
  # GET label combinations two by two
  diff_label <- combn(names(nets), m=2, FUN = paste0)
  
  #get metrics for each two by two combination and
  # according to base_comm (which one of the)
  comb <- lapply(seq(ncol(diff_label)), function(i){
    #print(diff_label[,i])  #sanity checks that labels are correct
    #print(diff_label[,i][base_comm]) #sanity checks that labels are correct
    lapply(diff_label[,i], function(j){
      community_based_metrics(nets[diff_label[,i]][[j]], memb[diff_label[,i]][[base_comm]],
                              metrics = metrics, normalized = normalized, sep_values = F)
    })
  })
  
  diffs <- list()
  diffs <- sapply(comb, function(x){ 
    data.frame(x[[2]]) - data.frame(x[[1]])
  }, simplify = FALSE)
  names(diffs) <- sapply(seq(ncol(diff_label)), function(x){ 
    paste(diff_label[2:1,x], collapse = " ~ ")
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
