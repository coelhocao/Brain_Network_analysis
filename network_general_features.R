#############################################################################################
# These codes have been costum designed for brain network analysis of 1-2 time point data
# generally acquired by immunolabelling of activity-dependent proteins (i.e. c-fos, Arc)
#
# Many functions here have been largely based on some of the codes shared by Justin Kenney,
#which can be found on this link https://github.com/jkenney9a/Networks
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
#########################################################################################################################

general_measures <- function(G, memb = NULL, mod_fun, metrics = c('mean_FC', 'gc_size', 'pct_gc_size', 'density', 'modularity'),
                             match_to_memb = TRUE, ...){
  #
  # Input: an igraph object or list of igraph objects
  #
  # Output: a dataframe with overall functional connectivity, size of giant component and link density
  #
  stopifnot(is_igraph(G))
  output <- data.frame(0)
  if('mean_FC' %in% metrics){ output <- data.frame(output, mean_FC = as.numeric(mean(E(G)$weight))) }
  if('gc_size' %in% metrics){ output <- data.frame(output, gc_size = as.numeric(max(components(G)[[2]]))) }
  if('pct_gc_size' %in% metrics){ output <- data.frame(output, pct_gc_size = as.numeric((max(components(G)[[2]])/vcount(G)) * 100)) }
  if('density' %in% metrics){ output <- data.frame(output, density = as.numeric(edge_density(G))) }
  

  if('modularity' %in% metrics){
    if(is.null(memb) & !is.null(V(G)$anatomical_groups)){
      memb <- list()
      memb$membership <- as.numeric(as.factor(V(G)$anatomical_groups))
      memb$names <- V(G)$name
    }
    if(isTRUE(match_to_memb)){
      if(length(memb$names[!(memb$names %in% V(G)$name)]) > 0){
        nv = length(memb$names[!(memb$names %in% V(G)$name)])
        G <- add_vertices(G,nv, attr =  list(name = c(memb$names[!(memb$names %in% V(G)$name)])))
      }
      if(length(V(G)$name[!(V(G)$name %in% memb$names)]) > 0){
        G <- delete.vertices(G, V(G)$name[!(V(G)$name %in% memb$names)])
      }
      output <- data.frame(output, modularity = mod_fun(G,memb$membership, ...))
    }
  }
  
  return(output[-1])
}

general_permutation <- function(df, memb = NULL, to_sample = 'signal', seed=1987, resample=10000, mod_fun = Q_signed,  
                                metrics = c('mean_FC', 'gc_size', 'pct_gc_size', 'density', 'modularity'),
                                match_to_memb = TRUE, distrib = F, ...){
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
  stopifnot(is.logical(match_to_memb) |
              is.logical(distrib))
  
  if(length(to_sample) > 1){
    df$labels$interac <- factor(with(df$labels, interaction(df$labels[,match(to_sample, names(df$labels))])))
    to_sample <- 'interac'
  }
  if(length(levels(df$labels[,to_sample])) < 2){
    stop("ERROR! Requires at least TWO levels in to_sample factor(s) to compare")
  }
  
  set.seed(seed)
  resamp_distrib <- sapply(seq(resample), function(i){
    gf_resamp(df, memb = memb, to_sample = to_sample, mod_fun = mod_fun, 
              metrics = metrics, match_to_memb = match_to_memb, ...)
  }, simplify = 'array')

  #empirical networks generation and differences
  emp <- net_list(df, factors = to_sample, ...)

  #Extracting measures
  emp_measures <- multi_net_measures(emp, sec_param = memb, FUN = general_measures, list_levels = 1, 
                                     add_param = 'memb', mod_fun = mod_fun, metrics = metrics,
                                     match_to_memb = match_to_memb)
  
  emp_diff_label <- combn(names(emp_measures),m=2, FUN = paste0)
  emp_diffs <- list()
  emp_comb <- combn(emp_measures,2, simplify = FALSE)
  emp_diffs <- lapply(emp_comb, function(x) abs(x[[1]] - x[[2]]))
  emp_diffs <- data.frame(rbindlist(emp_diffs))
  rownames(emp_diffs) <- sapply(1:ncol(emp_diff_label), function(x){ paste(emp_diff_label[,x], collapse = " ~ ")})
  
  #Calculating the p-values
  if(!is.null(ncol(resamp_distrib))){
    ps <- matrix(0, nrow(emp_diffs), ncol(emp_diffs))
    colnames(ps) <- colnames(emp_diffs)
    rownames(ps) <- rownames(emp_diffs)
    
    for(i in seq(nrow(emp_diffs))){
      for(j in seq(ncol(emp_diffs))){
        ps[i,j] <- mean(abs(resamp_distrib[i,j,] >=  emp_diffs[i,j]), na.rm = T)
      }
    } 
  }
  else{ 
    ps <- mean(abs(resamp_distrib >  emp_diffs[,1]), na.rm = T) 
  }
  
  if(isTRUE(distrib)){
    return(list(statistic = emp_diffs, pvalue = ps, random_distrib = resamp_distrib)) 
  }
  else{
    return(list(statistic = emp_diffs, pvalue = ps))
  }
}

gf_resamp <- function(df, memb = NULL, to_sample = 'group',mod_fun = Q_signed, 
                      metrics = c('mean_FC', 'gc_size', 'pct_gc_size', 'density', 'modularity'),
                      match_to_memb = TRUE, ...){
  
  #shuffle the grouping labels
  if(length(to_sample) > 1){
    df$labels$interac <- factor(with(df$labels, interaction(df$labels[,match(to_sample, names(df$labels))])))
    to_sample <- 'interac'
  }
  
  #some re-samplings end up generating NaNs in either the measurements or in the difference calculation
  # the while loop below ensures these resamplings are remade and doesn't go forward
  nonan <- 1
  while(nonan == 1){
    
    df$labels[[to_sample]] <- sample(df$labels[[to_sample]], replace = F)
    
    #randomized networks
    nets <- net_list(df, factors = to_sample, ...)
    
    # Nodal measures for the networks
    nets_measures <- multi_net_measures(nets, sec_param = memb, FUN = general_measures, list_levels = 1,
                                        add_param = 'memb', mod_fun = mod_fun, metrics = metrics,
                                        match_to_memb = match_to_memb)
    
    diff_label <- combn(names(nets_measures), m=2, FUN = paste0)
    diffs <- list()
    comb <- combn(nets_measures,2, simplify = FALSE)
    diffs <- lapply(comb, function(x){ abs(x[[1]] - x[[2]]) })
    diffs <- as.matrix(rbindlist(diffs))
    rownames(diffs) <- sapply(1:ncol(diff_label), function(x){ paste(diff_label[,x], collapse = " ~ ")})
    if(!any(is.nan(diffs))){
      nonan = 0
    }
  }
  
  return(diffs)
}
