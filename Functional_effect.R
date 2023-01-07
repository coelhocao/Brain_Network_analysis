
if (!('pracma' %in% installed.packages()[,'Package'])){install.packages("pracma")}; require(pracma)
if (!('psych' %in% installed.packages()[,'Package'])){install.packages("psych")}; require(psych)

#trocar esse nome por matrix_difference
functional_effect <- function(df, to_sample = c('signal', 'group'), resamp = FALSE, seed = 1987,
                              normalization = 'max_norm', nans = FALSE, ...){
  #
  # INPUT:
  #   df           List containing $counts as brain data, $labels containing of the factors (see
  #                mfiles_to_dataframe function)
  #   to_sample    Char defining the names of variables to use as factors to build the matrices and
  #                perform matrix difference(s). If more than one is provided, an additional var will be
  #                created as the interaction of all provided.
  #   resamp       Logical (FALSE, default). Whether to perform resampling of the levels of the factors
  #                provided in to_sample.
  #   fisherz      Logical (TRUE, default). Whether to re-scale the correlation matrices to z-score before
  #                calculating matrix difference.
  #   ...          additional args to be passed to corr_matrix_threshold function (ex. zeroed_regions, etc)
  #
  # OUTPUT:
  #                List of all combinations 2 x 2 of levels from factors in to_sample. Each element is 
  #                a matrix of the absolute difference between the two input matrices.
  
  stopifnot(is.logical(resamp))
  
  if(length(to_sample) > 1){
    df$labels$interac <- factor(with(df$labels, interaction(df$labels[,match(to_sample, names(df$labels))])))
    to_sample <- 'interac'
  }
  if(length(levels(df$labels[,to_sample])) < 2){
    stop("ERROR! Requires at least TWO levels in to_sample factor(s) to compare")
  }
  if(isTRUE(resamp)){
    df$labels[,to_sample] <- sample(df$labels[,to_sample], replace = FALSE)
  }
  
  # Leave only groups of the combination above in each list element
  mats <- lapply(levels(df$labels[,to_sample]), function(i){
    m <- df$counts[df$labels[,to_sample] == i,]
    m <- corr_matrix_threshold(m, normalization = normalization, nans = nans, ...)
    
    return(m)
  })
  names(mats) <- levels(df$labels[,to_sample])
  
  #All pairwise comparisons
  lab <- combn(names(mats), m = 2, FUN = paste0)
  # Compute matrix differences
  comb <- combn(mats, m=2, simplify = F)
  comps <- list()
  comps <- lapply(comb, function(i){ 
    i[[1]] - i[[2]]
  })
  names(comps) <- sapply(seq(ncol(lab)), function(x){ 
    paste(lab[,x], collapse = " ~ ")
  })
  return(comps)
}

functional_effects_permutation <- function(df, to_sample = c('signal', 'group'), seed=1987, 
                                           resample=10000, normalization = 'max_norm', nans = F, ...){
  #
  # INPUT:
  #   df           List containing $counts as brain data, $labels containing of the factors (see
  #                mfiles_to_dataframe function)
  #   to_sample    Char defining the names of variables to use as factors to build the matrices and
  #                perform matrix difference(s). If more than one is provided, an additional var will be
  #                created as the interaction of all provided.
  #   seed         Seed number for the resampling. important for reproducibility
  #   resample     Number of resamplings to perform
  #   fisherz      Logical (TRUE, default). Whether to re-scale the correlation matrices to z-score before
  #                calculating matrix difference.
  #   ...          additional args to be passed to corr_matrix_threshold function (ex. zeroed_regions, etc)
  
  
  set.seed(seed)
  resamp_effect <- lapply(1:resample, function(i){
    functional_effect(df, to_sample = to_sample, resamp = TRUE, normalization = normalization, nans = nans, ...)
  })
  
  n <- names(resamp_effect[[1]])
  resamp_effect <- simplify2array(resamp_effect)
  if(!is.null(nrow(resamp_effect))){
    resamp_effect <- split(resamp_effect, 1:nrow(resamp_effect))
    names(resamp_effect) <- n
    
    resamp_effect <- lapply(resamp_effect, function(i){
      array(unlist(i), dim=c(dim(i[[1]]), length(i)),
            dimnames = dimnames(i[[1]]))
    })
  }
  else{
    resamp_effect <- list(array(unlist(resamp_effect), dim=c(dim(resamp_effect[[1]]), length(resamp_effect)),
                           dimnames = dimnames(resamp_effect[[1]])))
    names(resamp_effect) <- n
  }
  
  emp_effect <- functional_effect(df, to_sample = to_sample, resamp = FALSE, normalization = normalization, nans = nans, ...)
  
  p_list <- lapply(emp_effect, function(i){ matrix(0, nrow(i), ncol(i)) })
  
  for(i in seq(p_list)){
    for(j in seq(ncol(p_list[[i]]))){
      for(z in seq(nrow(p_list[[i]]))){
        p_list[[i]][z,j] <- mean(abs(resamp_effect[[i]][z,j,]) > abs(emp_effect[[i]][z,j]))
      }
    }
  }
  
  output <- lapply(seq(length(emp_effect)), function(i){
    list(Functional_Effect = emp_effect[[i]], 
         p_value = p_list[[i]])
  })
  names(output) <- names(emp_effect)
  return(output)
}
