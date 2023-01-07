#############################################################################################
# These codes have been costumly designed for brain network analysis of 1-2 time point data
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

######## TODO
# calculated participation coefficient and shannon-entropy for signed network as defined by Rubinov & Sporns (2011)
# CONSTRUCTING A NETWORK BASED ON MUTUAL INFORMATION
# CONSTRUCTING A NETWORK BASED ON PARTIAL CORRELATION WITH BEHAVIOR ###  MAYBE NOT
# ON corr_matrix_threshold and all derivates
## TO BE DONE: i WANT TO SET THRESH.PARAM = 'NONE' GIVING OUT JUST THE MATRIX.
# THERE IS A PROBLEM WITH fISHER_Z AND NEGS. AS THEY CAN FUNCTION AS A THRESHOLDING OF EVERYTHING BELOW THE MEAN OF THE CORRELATIONS (BELOW ZERO).
# Create a function to delete vertices not present in another graph to be compraed to

########################################################################################################################
if (!('dplyr' %in% installed.packages()[,'Package'])){install.packages("dplyr")}; require(dplyr)
if (!('Hmisc' %in% installed.packages()[,'Package'])){install.packages("Hmisc")}; require(Hmisc) #better functions for correlation matrices and p-values
if (!('lattice' %in% installed.packages()[,'Package'])){install.packages("lattice")}; require(lattice) #for generating graphs, colors and matrices
if (!('data.table' %in% installed.packages()[,'Package'])){install.packages("data.table")}; require(data.table) #for rbindlist
if (!('igraph' %in% installed.packages()[,'Package'])){install.packages("igraph")}; require(igraph) #package for network generation and measures
if (!('boot' %in% installed.packages()[,'Package'])){install.packages("boot")}; require(boot)
if (!('psych' %in% installed.packages()[,'Package'])){install.packages("psych")}; require(psych) #fisherz

########################################################################################################################
###########  Correlation matrices and defining the graph  ##############################################################
########################################################################################################################

corr_matrix <- function(df, p.adjust.method = 'none', type='pearson', zeroed_regions = TRUE){
  #
  # Input: Dataframe
  #
  # Whether or not to apply a p-value adjustment method drawn from the 'p.adjust' function
  # (e.g, 'fdr', 'bonferroni' etc.). 
  # whether to use 'pearson' or 'spearman' correlation
  # NOTE: assume undirected graph!
  #   
  # Output: List of two dataframes corresponding to 1) all pairwise Pearson 
  # correlations and 2) all associated p-values, adjusted or not(default)
  
  # remove regions with no data across all subjects
  stopifnot(is.logical(zeroed_regions))
  
  if(isFALSE(zeroed_regions)){
    df <- df[,which(colSums(df) >0)]
  }
  
  corr <- rcorr(as.matrix(df), type=type)$r
  
  #adjust p-values if necessary
  p_values <- as.matrix(rcorr(as.matrix(df), type=type)$P)
  p_values[lower.tri(p_values)] <- p.adjust(p_values[lower.tri(p_values)], method=p.adjust.method) 
  names(p_values[])=colnames(corr)
  
  return(list("corr" = corr,"pvalue" = p_values))
}

boot_ratio_normalization <- function (dfs, brn_R = 1000, seed = 87) {
  #
  # INPUT: dataframe to construct correlation matrices, number of resampling (brn_R) and seed
  #
  # OUTPUT: bootstrapped ratio (empirical/Stand dev booted distrib) of the correlations
  # OBS: this matrix is to be used as the pvalue matrix in the thresholding process in for ex. the df_to_graph function
  # OBS2: occasionally, the resample may result in correlations of 1, which give FisherZ = Inf. Here, 
  # the highest score is give (4) for 1000 replicates.
  #
  ## Does it have any meaning at all?? Evaluate that and other options maybe the range of the sd (but what is the parameter)??
  
  # resampling function
  resamp <- function(dfs){
    d <- dfs[base::sample(nrow(dfs),nrow(dfs), replace = TRUE), ]
    r <- corr_matrix(d)[['corr']]
    r <- fisherz(r)
    return(as.matrix(r))
  }
  
  set.seed(seed)
  # Resampled and Z transformed correlations distribution
  rdf <- replicate(brn_R, resamp(dfs))
  
  # If resample corr = 1, fisher_z will be Inf, which returns an error
  # the value 4 to Inf values gives the highest score with 1000 replicates
  rdf[which(is.infinite(rdf))] <- 4 
  
  #Empirical Z transformed correlation matrix
  fdf <- fisher_z(corr_matrix(dfs)[['corr']])
  
  #bootstrap ratio empirical corrs/sd_booted_corrs
  sd_boot <- apply(rdf, c(1,2), sd)
  ratio <- as.matrix(fdf/sd_boot)
  
  ratio[which(is.nan(ratio))] <- 0
  ratio[which(is.infinite(ratio))] <- 0
  ratio[which(is.na(ratio))] <- 0
  
  return(as.data.frame(ratio))
}

corr_matrix_threshold <- function(dfs, negs = NULL, thresh = NULL, thresh.param = NULL, p.adjust.method = 'none', 
                                  type='pearson', brn_R = 1000, seed = 1987, nans = c(F,T,numeric), 
                                  zeroed_regions = FALSE, normalization = NULL){
  #
  # Input: 
  #  df              Dataframe with all the columns(vars) to be correlated
  #  negs            Defines handling of negative correlations. NULL (defatul) for nothing, 'zero' to set them to zero,
  #                  or 'abs' for absolute value.
  #  thresh          numeric vector of at least length 1 indicating a cut off point for correlations, based on a parameter
  #  thresh.param    Either 'none' (default), 'p', 'r' or 'cost' indicating parameter for the threshold. Notice that 
  #                  'p' will cut off values below thresh value, whereas 'r' and 'cost' will cut off values above thresh value
  #
  # zero (set to zero, default), abs (absolute value), sqrt (squareroot), resc (rescale to 0-2), or none
  # set threshold value AND parameter (r, p, cost, boot_ratio), the p-value adjustment method to use (if any),
  # the type of correlation (pearson or spearman).
  # UPDATE: threshold can now be a range of thresholds, a vector.
  # nans can be either logical or 'wea_links', which is intended to complete the graphs into a completely conneted graph
  # by adding weak_links to it.
  #
  # IMPORTANT: p.adjust.method WILL modify your p_value matrix used for thresh/thresh.param. BE AWARE!!
  #   
  # Output: [List of] Dataframe(s) of correlations thresholded at p < threshold(s).
  
  ## TO BE DONE: i WANT TO SET THRESH.PARAM = 'NONE' GIVING OUT JUST THE MATRIX.
  # THERE IS A PROBLEM WITH fISHER_Z AND NEGS. AS THEY CAN FUNCTION AS A THRESHOLDING OF EVERYTHING BELOW THE MEAN CORRELATIONS, which is ZERO after fisherZ transform.
  
#  if(thresh.param == 'boot_ratio' & p.adjust.method != 'none'){
#    p.adjust.method = 'none'
#    warning('thresh.param = boot_ratio, forcing p.adjust.method to none!')
#  }
  
  stopifnot(is.logical(zeroed_regions))
  
  stopifnot(is.null(normalization) |
              normalization == 'fisherz' |
              normalization == 'max_norm'|
              normalization == 'max_min_norm')
  stopifnot(is.null(negs) |
              negs == 'zero' |
              negs == 'abs')
  
  # remove regions with no data across all subjects
  if(isFALSE(zeroed_regions)){
    dfs <- dfs[,which(colSums(dfs) > 0)]
  }
  
  df <- corr_matrix(dfs, p.adjust.method=p.adjust.method, type=type)
  df_corr <- df$corr
  df_pvalue <- df$pvalue
  
  #remove any NaNs, infs or NAs (sometimes happens with bootstrapping; not sure why)
  #NaNs also happen when all data is zero, then corrs are NaNs and p_value is inf
  #The option of weak_link can force a fully connected graph
  #When applying thresholds, use nans = FALSE
    if(is.logical(nans)){
      if(isFALSE(nans)){
        df_pvalue[mapply(is.infinite, df_pvalue)] <- 1
        df_pvalue[mapply(is.nan, df_pvalue)] <- 1  
        df_pvalue[mapply(is.na, df_pvalue)] <- 1
        df_corr[mapply(is.infinite, df_corr)] <- 1
        df_corr[mapply(is.nan, df_corr)] <- 0  
        df_corr[mapply(is.na, df_pvalue)] <- 0
        df_corr[mapply(is.na, df_corr)] <- 0
      }  
    }
    else{
      if(!is.null(thresh.param) | !is.null(thresh)){
        stop("weak_link is intended to force a fully connected network\n
             which is not compatible with thresh and thresh.param args \n
             being different than NULL!")
      }
      if(is.numeric(nans)){
        df_corr[mapply(is.infinite, df_corr)] <- 1
        df_corr[mapply(is.nan, df_corr)] <- 1/(10^nans)
        df_corr[mapply(is.na, df_corr)] <- 1/(10^nans)
      }
    }

  #What to do with negative correlations
  #BE CAREFUL when using Fisher Z transofm here. The results are very different than when using correlations.
  if(!is.null(negs)){
    if(negs == 'zero'){
      df_corr <- df_corr * (df_corr > 0)
    }
    else if(negs == 'abs'){
      df_corr <- abs(df_corr)
    }
  }
  
  if(!is.null(thresh) & !is.null(thresh.param)){
    #Bootstrap ratio matrix
    # NOT WORKING YET FOR SOME REASON
    #if(tolower(thresh.param) == 'boot_ratio'){
    #  ratio <- boot_ratio_normalization(dfs, brn_R = brn_R, seed = seed)
    #}
    
    if(tolower(thresh.param)=='p'){#apply p-value threshold to correlation matrix
      df_corr <- lapply(thresh, function(x) {df_corr[mapply(">=", df_pvalue, x)] <- 0; df_corr})
    }
    else{
      if(tolower(thresh.param)=='r'){ # to the coefficient itself
        df_corr <- lapply(thresh, function(x) {df_corr[mapply("<=", df_pvalue, x)] <- 0; df_corr})
      }
      else{
        if(tolower(thresh.param)=='cost'){ #to a proportional thresholding
          r.threshold <- quantile(abs(df_corr), probs=1-thresh, na.rm=TRUE)
          df_corr <- lapply(r.threshold, function(x) {df_corr[mapply("<=", df_pvalue, x)] <- 0; df_corr})
        }
        #else{
        #  if(tolower(thresh.param)=='boot_ratio'){ #based on resampling
        #    df_corr <- lapply(thresh, function(x) {df_corr[mapply("<=", ratio, x)] <- 0; df_corr})
        #  }
        #  else{
        #    stop("Invalid thresholding parameter") 
        #  }
        #}
      }
    }
  }
  
  if(length(thresh) > 1){
    names(df_corr) <- thresh
  }
  
  #Fisher Z transform for correlations
  if(!is.null(normalization)){
    if(normalization == 'fisherz'){
      df_corr <- psych::fisherz(df_corr)
    }
    else{
      if(normalization == 'max_norm'){
        df_corr <- df_corr/max(df_corr)
      }
      else{
        if(normalization == 'max_min_norm'){
          df_corr <- (df_corr - min(df_corr)) / (max(df_corr) - min(df_corr))
        }
      }
    } 
  }
  
  return(df_corr)
} 

threshold <- function(mat, p.mat=NULL, thresh = 0.05, thresh.param = 'p'){
  #Input a correlation matrix, a p.value matrix (thresh param=p), a threshold, and a threshold parameter
  # e.g, (r, p, cost or none)
  #
  #Output: the correlation matrix thresholded as appropriate
  
  #Bootstrap ratio matrix
  # NOT WORKING YET FOR SOME REASON
  
  if(tolower(thresh.param)=='p'){#apply p-value threshold to correlation matrix
    df_corr <- lapply(thresh, function(x) {df_corr[mapply(">=", df_pvalue, x)] <- 0; df_corr})
  }
  else{
    if(tolower(thresh.param)=='r'){ # to the coefficient itself
      df_corr <- lapply(thresh, function(x) {df_corr[mapply("<=", df_pvalue, x)] <- 0; df_corr})
    }
    else{
      if(tolower(thresh.param)=='cost'){ #to a proportional thresholding
        r.threshold <- quantile(abs(df_corr), probs=1-thresh, na.rm=TRUE)
        df_corr <- lapply(r.threshold, function(x) {df_corr[mapply("<=", df_pvalue, x)] <- 0; df_corr})
      }
      else{
        stop("Invalid thresholding parameter") 
      }
    }
  }
  
  return(mat)
}

file_to_graph <- function(file, clean_data = F, negs = NULL, thresh = NULL, thresh.param = NULL, p.adjust.method = 'none', type='pearson', 
                          brn_R = 1000, seed = 87, weighted = TRUE, nans = FALSE, zeroed_regions = FALSE,  normalization = NULL){
  #
  #Input: csv, txt file
  #
  # cleans the data, removes or normalize negative correlations:
  # zero (set to zero, default), abs (absolute value), sqrt (squareroot), resc (rescale to 0-2), or none
  # set threshold value AND parameter (r, p, or cost), the p-value adjustment method to use (if any),
  # the type of correlation (pearson or spearman) and if graph should be waited or not.
  #
  # Output: igraph graph
  
  df <- load_data(file)
  if(isTRUE(clean_data)){
    df <- clean_data(df)
  }
  
  if(any(c(!is.null(negs), !is.null(thresh), !is.null(thresh.param), p.adjust.method != 'none',
           isTRUE(fisherz), isFALSE(nans), nans == 'weak_link',  isTRUE(max_norm) ))){
    df_G <- corr_matrix_threshold(df, negs = negs, thresh=thresh, thresh.param=thresh.param, 
                                  p.adjust.method=p.adjust.method, type=type,
                                  brn_R = brn_R, seed = seed, nans = nans, zeroed_regions = zeroed_regions,
                                  normalization =  normalization)
  }
  else{
    df_G <- corr_matrix(df, type = type, zeroed_regions = zeroed_regions)$corr
  }
  
  
  
  G <- graph.adjacency(as.matrix(df_G), mode="undirected", weighted=weighted, diag = F)
  return(G)
}


df_to_graph <- function(df, clean_data = FALSE, negs = NULL, thresh = NULL, thresh.param = NULL, p.adjust.method = 'none', type='pearson', 
                        brn_R = 1000, seed = 87, weighted = TRUE, nans = FALSE, zeroed_regions = FALSE,  normalization = NULL,
                        annotated = T, to_annotate = cbind(anatomical_groups, colors)){
  #
  # Input: dataframe var
  #
  # cleans the data, removes or normalize negative correlations:
  # zero (set to zero, default), abs (absolute value), sqrt (squareroot), resc (rescale to 0-2), or none
  # set threshold value AND parameter (r, p, or cost), the p-value adjustment method to use (if any),
  # the type of correlation (pearson or spearman) and if graph should be waited or not.
  # NOTE: RN THIS FUNCTION IS SET TO BE USED ONLY WITH A SINGLE THRESHOLD AND TO BE CALLED BY NETS_LIST AS MULTIPLE THRESH LIST OF GRAPHS
  # Output: igraph graph
  
  df <-as.data.frame(df)
  
  if(isTRUE(clean_data)){
    df <- clean_data(df)
  }
  
  if(any(c(!is.null(negs), !is.null(thresh), !is.null(thresh.param), p.adjust.method != 'none',
     isFALSE(nans), nans == 'weak_link', !is.null(normalization)))){
    df_G <- corr_matrix_threshold(df, negs = negs, thresh=thresh, thresh.param=thresh.param, 
                                  p.adjust.method=p.adjust.method, type=type,
                                  brn_R=brn_R, seed=seed, nans = nans, zeroed_regions = zeroed_regions,
                                  normalization = normalization)
  }
  else{
    df_G <- corr_matrix(df, type = type, zeroed_regions = zeroed_regions)$corr
  }
  
  G <- graph.adjacency(as.matrix(df_G), mode = "undirected", weighted = weighted, diag=FALSE)
  
  if(isTRUE(annotated)){
    G <- set.vertex.attribute(G, name = 'counts', 
                              value = sapply(seq_along(df), function(z){ mean(df[,z]) }))
    
    # Am getting an error here whenever I exclude regions using zeroed_regions = FALSE or nans = F
    # Some anatomical or color designations may be wrong.
    G <- set.vertex.attribute(G, name = 'anatomical_groups', value = as.character(to_annotate[,1]))
    G <- set.vertex.attribute(G, name = 'colors', value = as.character(to_annotate[,2]))
    
  }
  return(G)
}

net_list <- function(datalist, factors = c('group', 'signal'), thresh = NULL, annotated = T, quant = 'counts',
                     clean_data=F, negs = NULL, thresh.param = NULL, p.adjust.method = NULL, type='pearson', 
                     brn_R = 1000, seed = 87, weighted = TRUE, nans = FALSE,
                     zeroed_regions = TRUE,  normalization = NULL, nested = FALSE){
  # INPUT: a datalist with $labels and $counts dataframes, and factors (within $labels df) to put separated in different lists
  # OUTPUT: a list of lists of networks. May come with mean_signal annotated as an vertex attribute
  # NOTES: rn to_annotate is set automatically to anatomical_groups and colors from the datalist$annotation list
  
  # generating the networks
  datalist$labels$interac <- factor(with(datalist$labels, interaction(datalist$labels[,match(factors, names(datalist$labels))])))
  nets <- list()
  nets <- lapply(levels(datalist$labels$interac), function(i){
    if(is.null(thresh)){
      df_to_graph(datalist[[quant]][datalist$labels$interac == i,], to_annotate = data.frame(datalist$annotation$anatomical_groups, datalist$annotation$colors), 
                  negs = negs, thresh = thresh, thresh.param = thresh.param, p.adjust.method = p.adjust.method, type = type,
                  brn_R = brn_R, seed = seed, weighted = weighted, annotated = annotated, nans = nans,
                  zeroed_regions = zeroed_regions,  normalization = normalization)
    }
    else{
      lapply(thresh, function(j){
        df_to_graph(datalist[[quant]][datalist$labels$interac == i,], to_annotate = data.frame(datalist$annotation$anatomical_groups, datalist$annotation$colors), 
                    negs = negs, thresh = j, thresh.param = thresh.param, p.adjust.method = p.adjust.method, type = type,
                    brn_R = brn_R, seed = seed, weighted = weighted, annotated = annotated, nans = nans,
                    zeroed_regions = zeroed_regions, normalization = normalization)
      }) 
    }
  })
  names(nets) <- levels(datalist$labels$interac)
  
  if(!is.null(thresh)){
    nets <- lapply(nets, function(i){
      names(i) <- thresh
      i
    })
    if(isFALSE(nested)){
      nets <- unlist(nets, recursive = FALSE)
    }
  }
  
  return(nets)
}


del_isolated <- function(g){
  #
  # Deletes non-connected nodes of the network
  
  isolated <- which(degree(g)==0)
  g2 = delete.vertices(g, isolated)
  return(g2)
}

