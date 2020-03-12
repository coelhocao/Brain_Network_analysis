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


########  THINGS TO BE ADDED
# calculated participation coefficient and shannon-entropy for signed network as defined by Rubinov & Sporns (2011)
# CONSTRUCTING A NETWORK BASED ON MUTUAL INFORMATION
# CONSTRUCTING A NETWORK BASED ON PARTIAL CORRELATION WITH BEHAVIOR ###  MAYBE NOT
# make networks annotated with an anatomical groups of regions and color codes.

if (!('dplyr' %in% installed.packages()[,'Package'])){install.packages("dplyr")}; require(dplyr)
if (!('Hmisc' %in% installed.packages()[,'Package'])){install.packages("Hmisc")}; require(Hmisc) #better functions for correlation matrices and p-values
if (!('lattice' %in% installed.packages()[,'Package'])){install.packages("lattice")}; require(lattice) #for generating graphs, colors and matrices
if (!('data.table' %in% installed.packages()[,'Package'])){install.packages("data.table")}; require(data.table) #for rbindlist
if (!('igraph' %in% installed.packages()[,'Package'])){install.packages("igraph")}; require(igraph) #package for network generation and measures
if (!('boot' %in% installed.packages()[,'Package'])){install.packages("boot")}; require(boot)
#if (!('multcomp' %in% installed.packages()[,'Package'])){install.packages("multcomp")}; require(multcomp)
#if (!('heplots' %in% installed.packages()[,'Package'])){install.packages("heplots")}; require(heplots)
#if (!('effsize' %in% installed.packages()[,'Package'])){install.packages("effsize")}; require(effsize)
#if (!('statGraph' %in% installed.packages()[,'Package'])){install.packages("statGraph")}; require(statGraph) #For Jensen-Shannon divergence between two graphs
#if (!('brainGraph' %in% installed.packages()[,'Package'])){install.packages("braingraph")}; require(brainGraph) #participation and Gateway and Diversity coefficients

#####################################################################################
#######################          Nodal Measures          ############################
#####################################################################################

signed_strength <- function(g){
  #
  # INPUT: an igraph object
  # OUTPUT: normalized strength considering negative weights
  # NOTE: this function employes the rationale presented in Rubinov & Sporns (2011).
  
  stopifnot(is_igraph(g) | is.matrix(g))
  
  if(is_igraph(g)){ 
    A <- as.matrix(g[]) 
  } else {
    A <- g
  }
  
  n <- ncol(A);
  Spos <- A;
  Spos[Spos < 0] <- 0;
  Spos <- rowSums(Spos);  # positive strength
  
  Sneg <- A;
  Sneg[Sneg > 0] <- 0;
  Sneg <- rowSums(abs(Sneg));  # negative strength
  
  nSpos <- Spos*(1/(n-1));
  nSneg <- Sneg*(1/(n-1));
  
  signed <- nSpos - (Sneg/(Spos + Sneg))*nSneg  # signed strength according to Rubinov & Sporns (2011)
  signed <- sapply(signed, function(x){ ifelse(is.nan(x), 0, x) })
  
  return(list(Spos = Spos, Sneg = Sneg, signed_strength = signed))
}

leverage <- function(g, weighted = FALSE, signed = FALSE){
  #
  # Input: igraph graph
  #
  # Calculates centrality as defined in  
  # Joyce et al (2010). A New Measure of Centrality for Brain Networks (PLoS ONE 5(8):e12200)
  # code mostly taken from igraph wiki igraph.wikidot.com/r-recipes#toc10
  # CAUTION: The options of putting weighted and signed in the leverage calculation are not predicted in the paper, and
  # were not fully evaluated in regards to its behavior.
  #
  # Output: Leverage centrality
  #
  
  if(isFALSE(weighted) & isFALSE(signed)){
    k <- degree(g)
  } else{
    if(isTRUE(weighted) & isFALSE(signed)){
      k <- strength(g)
    } else{
      if(isTRUE(weighted) & isTRUE(signed)){
        k <- signed_strength(g)
      }
    }
  }
  n <- vcount(g)
  
  lev <- sapply(1:n, function(v) { mean((k[v]-k[neighbors(g,v)]) / (k[v]+k[neighbors(g,v)])) })
  
  #when k = 0, leverage = NA or NaN. As this can compromise permutations,
  #and both negative and positive values are meaningful, the closest interpretation is -1
  lev[lev %in% NA] <- -1
  lev[lev %in% NaN] <- -1
  names(lev) <- V(g)$name
  
  return(lev) 
}


Nodal_efficiency <- function(G, normalized=FALSE) {
  # Input: igraph graph, weights will be used if existent
  # whether or not to normalize (where max value = 1)
  
  # Output: data frame of nodal efficiency (i.e, average inverse shortest path length 
  # from individual node to all other nodes)
  
  if(ecount(G) > 0){
    E(G)$weight <- 1 - abs(E(G)$weight) #calculate distance matrix from correlation matrix b/c
  }
  eff <- 1 / distances(G)
  
  eff[!is.finite(eff)] <- 0
  out <- colSums(eff) / (vcount(G) - 1)
  
  if(normalized == TRUE){
    out <- out / max(out)
  }
  return(out)
}

efficiency_loss <- function(g){
  # Input: igraph graph
  #
  # removes a node, calculate the Global and mean local efficiency relative to original graph's values
  # NOTE: takes a LONG TIME. it seems it's not optimized. Better not use with resampling procedures.
  #
  # Output: dataframe with the above calculation for each region
  
  rem <- function(net,r) {
    # Remove vids vertices by random
    vids<-V(net)[r]
    net <- delete.vertices(net,vids)
    return(net)
  }
  
  rem_calc <- function(g,r){
    
    
    gr <- rem(g,r)
    grgef <- gl_eff(gr)
    grlef <- mean(local_eff(gr))
    
    output <- c(grgef, grlef)
    return(output)
  }
  
  m <- as.data.frame(t(sapply (1:vcount(g), FUN= function (x) rem_calc(g,x))))
  geg <- gl_eff(g)
  leg <- mean(local_eff(g))
  if(geg == 0){ m[,1] <- rep(0, length(m[,1]))}
  else{
    m[,1] <- 1 - (m[,1]/geg)
  }
  if(leg == 0){ m[,2] <- rep(0, length(m[,1]))}
  else{
    m[,2] <- 1 - (m[,2]/leg)
  }
  colnames(m) <- c("gleff0/gleffi", "leff0/leffi")
  
  return(m)
}

Nodal_measures <-function(G, normalized=FALSE, efficiency_loss= FALSE, local_eff=FALSE, min_max_normalization=FALSE,
                          metrics = c("degree", "strength", "signed_strength", "eigenvector","betweenness","burt","closeness", "leverage", 
                                      "transitivity", "Nodal_eff")){
  # Input: An igraph graph
  #
  # Output: A dataframe of centrality measures:
  # degree, strength, betweenness, eigenvector, closeness, Burt's constraint, leverage,
  # transitivity, nodal efficiency,
  #
  # Optional: local efficiency, loss of efficiency (global and mean local) after region removal
  # provides how much the network loses of Gl and Lc efficiency by losing a region (different than nodal efficiency)
  # HOWEVER, this calculation takes a really long time. Need to be optimized!
  #
  # NOTE: nodal efficiency makes use of weights if they exist. For unweighted measures, create unweighted graphs!!
  # Optional: normalize all values according to min and max values
  #
  #
  
  if(ecount(G) == 0){
    zeros <- rep(0, vcount(G))
    output <- as.data.frame(cbind("degree"=zeros, "strength"=zeros, "signed_strength"=zeros, "eigenvector"=zeros,"betweenness"=zeros, 
                                  "burt"=zeros,"closeness"=zeros, "leverage" = zeros, "transitivity"=zeros,
                                  "Nodal_eff" = zeros), row.names=V(G)$name)
    if(local_eff == TRUE){output$local_eff = zeros}
    if(efficiency_loss == TRUE){
      output$geff_loss = zeros 
      output$leff_loss = zeros
    }
    return(output)
  }
  
  
  G.pos <- G
  E(G.pos)$weight <- abs(E(G.pos)$weight) #Positive numbers are necessary for betweenness and closeness
  E(G.pos)$weight <- 1 - E(G.pos)$weight #bet and clo are based on distance graphs/matrices, not weights. here, dist= 1 - corr
  
  # I had some issues with betweenness() saying weight vector should be positive,
  # During permutation, some iterations may result in r = or ~ 1, which
  # The last one is that, as 0 is not accepted either, a value irrelevant enough (a millionth) is provided as approximation
  
  E(G.pos)$weight[which(is.na(E(G.pos)$weight))] <- 1
  E(G.pos)$weight[which(is.nan(E(G.pos)$weight))] <- 1
  E(G.pos)$weight[which(is.infinite(E(G.pos)$weight))] <- 1
  E(G.pos)$weight[E(G.pos)$weight %in% 0] <- 0.0000000001
  
  
  #Centrality measures
  output <- data.frame(V(G)$name)
  if('degree' %in% metrics){ output <- data.frame(output, degree = igraph::degree(G, normalized=normalized)) }
  if('strength' %in% metrics){  output <- data.frame (output, strength = strength(G))  }
  if('signed_strength' %in% metrics){  output <- data.frame (output, signed_strength = signed_strength(G))  }
  if('eigenvector' %in% metrics){  output <- data.frame (output, eigenvector = evcent(G)$vector)  }
  if('betweenness' %in% metrics){  output <- data.frame (output, betweenness = betweenness(G.pos, normalized=normalized, weights = E(G.pos)$weight))  }
  if('closeness' %in% metrics){  output <- data.frame (output, closeness = closeness(G.pos, normalized=normalized))  }
  if('burt' %in% metrics){  output <- data.frame (output, burt = constraint(G))  }
  if('leverage' %in% metrics){  output <- data.frame (output, leverage = leverage(G))  }
  if('transitivity' %in% metrics){  output <- data.frame (output, transitivity = transitivity(G,type="local", isolates='zero'))  }
  if(tolower('Nodal_eff') %in% metrics){  output <- data.frame (output, Nodal_eff = Nodal_efficiency(G, normalized = normalized))  }
  
  #when degree = 0, burt = NaN. Set it to a value higher than max since lower values indicate more centrality
  if('burt' %in% metrics) {output$burt <- output$burt[!is.finite(output$burt)] <- max(!is.finite(output$burt)) + 0.1 }
  
  output <- output[,2:ncol(output)]
  # include nodal local eff (different from nodal_eff) in the dataframe
  if(local_eff == TRUE){
    
    leff <- local_eff(G)
    output$local_eff <- leff
  }
  
  #include loss of efficiency (global and local) after that node's removal (different from nodal_eff and local_eff) 
  if(efficiency_loss == TRUE){
    
    loss <- efficiency_loss(G)
    output$geff_loss <- loss[,1]
    output$leff_loss <- loss[,2]
  }
  
  if(min_max_normalization==TRUE){
    output <- apply(output, MARGIN=2, function(x) {(x-min(x)) / (max(x) - min(x))})
    output <- as.data.frame(output)
    if(!is.null(output$transitivity)){
      output$transitivity <- transitivity #Do not normalize transitivity b/c it is already normalized
    }
    if(!is.null(output$eigenvector)){
      output$eigenvector <- eigenvector #eigenvector is also a normalized index ranging from 0-1
    }
    if(!is.null(output$leverage)){
      output$leverage <- leverage #leverage is also normalized from -1 to 1.
    }
  }
  
  return(output)
}

#######################################################################################################################
#######################   End of Nodal Measures          ##############################################################
#######################################################################################################################


multi_net_measures <- function(G_list, sec_param, FUN = cluster_louvain, nested_levels = 2, add_param = NULL, ...){
  #
  # INPUT: A graph, a list of graphs or a list of lists of graphs, 
  # A function to apply over an igraph object and whether to delete isolated nodes
  #
  # OUTPUT: Function applied over graphs. Used here for centrality measures and community classifications
  
  if (nested_levels == 0){
    
    if(is.null(add_param)){
      output <- FUN(G_list, ...)
    }
    else{
      params <- list(G_list, sec_param, ...)
      names(params)[2] <- add_param
      output <- do.call(FUN, params)
      
    }
    
  }
  else{
    if (nested_levels == 1){
      
      output <- list()
      if(is.null(add_param)){
        output <- lapply(G_list, function(j){ 
          FUN(j, ...) 
        })
      }
      else{
        output <- list()
        output <- lapply(seq(G_list[[i]]), function(j){
          params <- list(G_list[[i]], sec_param[[i]], ...)
          names(params)[2] <- add_param
          output <- do.call(FUN, params)
          
        })
      }
    }
    else{
      if(nested_levels == 2){
        output <- list()
        if(is.null(add_param)){
          output <- lapply(G_list, function(i){
            lapply(i, function(j){
              FUN(j, ...) 
            })
          })
        }
        
        else{
          output <- list()
          output <- lapply(seq(G_list), function(i){
            lapply(seq(G_list[[i]]), function(j){
              params <- list(G_list[[i]][[j]], sec_param[[i]][[j]], ...)
              names(params)[2] <- add_param
              output <- do.call(FUN, params) 
            })
          })
        }
      }
    }
  }
  
  return(output)
}

measures_scores <- function(a, thresh = expr(mean(as.numeric(a)) + sd(as.numeric(a)))){
  #
  # Input: vector with nodal measures
  #
  # Output: nodes above the given threshold are given the value 1 in a new vector
  # obs: threshold is a formula based on the dataframe or a value
  #alternatively can do 20% highest a[order(a, decreasing = T)][1:(length(a)*0.2)]
  s <- eval(thresh)
  d <- a[a >= s]
  s <- rep(0, length(a))
  s[a %in% d] <- 1
  return(s)
}

hubness_ <- function(a, thresh =  expr(mean(as.numeric(a)) + sd(as.numeric(a))), score = T){
  #
  # Input: data frame with nodal measures in each column
  #
  # Output: Sums the number of measures in which a given node is above the threshold
  #
  
  d <- sapply(a, function(i){  measures_scores(i, thresh = thresh)  })

  
  rownames(d) <- rownames(a)
  
  if(score == T){
    d <- cbind(d, hub_score = rowSums(d)) 
  }
  
  return(as.data.frame(d))
}

resamp <- function(df, labels,negs = 'zero', thresh = 0.05, thresh.param = 'p', 
                   p.adjust.method = 'none', type = 'pearson', weighted = T, normalized = FALSE,
                   efficiency_loss = FALSE, local_eff = FALSE, min_max_normalization = FALSE, nans = 'zero',
                   metrics = c('degree', 'strength', 'signed_strength', 'eigenvector', 
                               'betweenness', 'leverage', 'Nodal_eff')){
  
  #shuffle the grouping labels
  newlabels <- factor(sample(labels, replace=F))
  
  #randomized networks
  nets <- list() 
  nets <- lapply(seq_along(levels(newlabels)), function(x){
    nets$x <- df_to_graph(df[which(newlabels == levels(newlabels)[x]),], negs=negs, thresh=thresh,
                          thresh.param=thresh.param, p.adjust.method=p.adjust.method, type=type, 
                          weighted=weighted, annotated = F, to_annotate = NULL, nans = nans)
  })
  
  # Nodal measures for the networks
  # Uses multi_net_measures for list with groups instead of list with different thresholds
  nets_measures <- multi_net_measures(nets, FUN = Nodal_measures, normalized = normalized, efficiency_loss = efficiency_loss, 
                                  metrics = metrics,local_eff = local_eff, min_max_normalization = min_max_normalization, nested_levels = 1)
  
  names(nets_measures) <- levels(labels)
  
  diff_label <- combn(names(nets_measures), m=2, FUN = paste0)
  diffs <- list()
  comb <- combn(nets_measures,2, simplify = FALSE)
  diffs <- lapply(comb, function(x) abs(x[[1]] - x[[2]]))
  
  dim_names <- list()
  dim_names$dim1 <- rownames(nets_measures[[1]])
  dim_names$dim2 <- names(diffs[[1]])
  dim_names$dim3 <- sapply(1:ncol(diff_label), function(x){ paste(diff_label[,x], collapse = "-")})
  diffs <- array(as.numeric(unlist(diffs)), dim=c(length(diffs[[1]][,1]), length(diffs[[1]]), 
                                                  length(diffs)), dimnames = dim_names)
  
  return(diffs)
}

nodal_permutation_test <- function(df, labels, seed = 1235, resample = 1000, negs = 'none', thresh = 1, 
                                   thresh.param='p', p.adjust.method='none', type='pearson', weighted=TRUE,
                                   normalized=FALSE, efficiency_loss= FALSE, local_eff=FALSE, min_max_normalization=FALSE,
                                   metrics = c('strength', 'eigenvector', 'betweenness', 'Nodal_eff'), nans = 'zero'){
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
  
  
  labels <- as.factor(labels)

  
  
  if(length(levels(labels)) < 2){
    stop("ERROR! It's a between network design. Requires at least TWO labels necessary to compare")
  }
  
  set.seed(seed)
  resamp_distrib <- replicate(resample, resamp(df$counts,labels,negs=negs, thresh=thresh, thresh.param=thresh.param,
                                               p.adjust.method=p.adjust.method, type=type, weighted=weighted, nans = nans,
                                               normalized=normalized,efficiency_loss= efficiency_loss, metrics = metrics,
                                               local_eff=local_eff, min_max_normalization=min_max_normalization))
  
  #empirical networks generation and differences
  emp <- list()
  emp <- lapply(seq_along(levels(labels)), function(x){
    emp$x <- df_to_graph(df$counts[which(labels == levels(labels)[x]),], negs=negs, thresh=thresh,
                         thresh.param=thresh.param, p.adjust.method=p.adjust.method, type=type, 
                         weighted=weighted, annotated = F, nans = nans)
  })
  
  
  
  emp_measures <- multi_net_measures(emp, FUN = Nodal_measures, normalized = normalized, efficiency_loss = efficiency_loss, 
                                 metrics = metrics,local_eff = local_eff, min_max_normalization = min_max_normalization, nested_levels = 1)
  
  names(emp_measures) <- levels(labels)
  
  emp_diff_label <- combn(names(emp_measures),m=2, FUN = paste0)
  emp_diffs <- list()
  emp_comb <- combn(emp_measures,2, simplify = FALSE)
  emp_diffs <- lapply(emp_comb, function(x) abs(x[[1]] - x[[2]]))
  
  dim_names <- list()
  dim_names$dim1 <- rownames(emp_measures[[1]])
  dim_names$dim2 <- names(emp_measures[[1]])
  dim_names$dim3 <- sapply(1:ncol(emp_diff_label), function(x){ paste(emp_diff_label[,x], collapse = "-")})
  
  emp_diffs <- array(as.numeric(unlist(emp_diffs)), dim=c(length(emp_diffs[[1]][,1]), length(emp_diffs[[1]]),
                                                          length(emp_diffs)), dimnames = dim_names)
  
  #Calculating the p-values
  ps <- array(0, dim = dim(emp_diffs))
  
  for(x in 1:dim(emp_diffs)[3]){
    for(y in 1:dim(emp_diffs)[2]){
      for(z in 1:dim(emp_diffs)[1]){
        ps[z,y,x] <- mean(abs(resamp_distrib[z,y,x,]) > abs(emp_diffs[z,y,x]))
      }
    }
  }
  
  dimnames(ps) <- dimnames(emp_diffs)
  
  return(list(diffs = emp_diffs, pvalue = ps))
}

################################################################################################################################
###############################################################################################################################
##################         BEWARE,,,,,  FUNCTIONS BEYOND HERE DO NOT WORK WELL








node_diff_count <- function(nets_nodal_measures, node_comparisons, multi_thresholds = F){
  #
  #   THIS FUNCTIONS IS NOT WELL DEFINED. 
  # INPUT: output of nodal_measures function and output of nodal_permutation_test
  #
  # OUTPUT: matrices to be used for a heatmap of significant differences in centrality comparisons
  if(multi_thresholds == T){
    p_array <- lapply(node_comparisons, function(x){ x$pvalue })
    names(p_array) <- names(node_comparisons)
    
    p_array <- array(unlist(p_array), dim = c(nrow(p_array[[1]]), ncol(p_array[[1]]), length(p_array)), 
                     dimnames = list(rownames(p_array[[1]]), colnames(p_array[[1]]), names(p_array)))
    
    highers <- array(0, dim = dim(p_array), dimnames = dimnames(p_array))
    highers <- sapply(seq(length(nets_nodal_measures[[1]])), function(x){
      sapply(seq(ncol(nets_nodal_measures[[1]][[1]])), function(y){
        sapply(seq(nrow(nets_nodal_measures[[1]][[1]])), function(z){
          highers[z,y,x] <- ifelse(nets_nodal_measures[[1]][[x]][z,y] > nets_nodal_measures[[2]][[x]][z,y], 1, 
                                   ifelse(nets_nodal_measures[[1]][[x]][z,y] < nets_nodal_measures[[2]][[x]][z,y], -1, 0))
        }, simplify = "array")  
      }, simplify = "array")
    }, simplify = "array")
    dimnames(highers) <- dimnames(p_array)
    
    
    cfos_highers <- sapply(seq(dim(highers)[1]), function(i){
      sapply(seq(dim(highers)[2]), function(j){
        cfos_highers <- abs(sum(highers[i,j,p_array[i,j,] <= 0.05] ==1))
      })
    })
    cfos_highers <- t(cfos_highers)
    colnames(cfos_highers) <- colnames(highers)
    rownames(cfos_highers) <- rownames(highers)
    
    cfos_highers <- cbind(cfos_highers, sum_score = sapply(seq(nrow(cfos_highers)), function(i){sum(abs(cfos_highers[i,])>1)}))
    
    rfp_highers <- sapply(seq(dim(highers)[1]), function(i){
      sapply(seq(dim(highers)[2]), function(j){
        rfp_highers <- abs(sum(highers[i,j,p_array[i,j,] <= 0.05] == -1))
      })
    })
    rfp_highers <- t(rfp_highers)
    colnames(rfp_highers) <- colnames(highers)
    rownames(rfp_highers) <- rownames(highers)
    rfp_highers <- cbind(rfp_highers, sum_score = sapply(seq(nrow(rfp_highers)), function(i){sum(abs(rfp_highers[i,])>1)}))
  }
  else{
    p_array <- node_comparisons$pvalue 
    
    highers <- sapply(seq(ncol(nets_nodal_measures[[1]])), function(y){
      sapply(seq(nrow(nets_nodal_measures[[1]])), function(z){
       ifelse(nets_nodal_measures[[1]][z,y] > nets_nodal_measures[[2]][z,y], 1, 
                                ifelse(nets_nodal_measures[[1]][z,y] < nets_nodal_measures[[2]][z,y], -1, 0))
      })  
    })
    rownames(highers) <- rownames(p_array)
    colnames(highers) <- colnames(p_array)
    
    cfos_highers <- sapply(seq(nrow(highers)), function(i){ 
      sapply(seq(ncol(highers)), function(j){
        ifelse(highers[i,j] == 1 & p_array[i,j,1] < 0.05,1,0)
      })
    })
    cfos_highers <- t(cfos_highers)
    colnames(cfos_highers) <- colnames(highers)
    rownames(cfos_highers) <- rownames(highers)
    cfos_highers <- cbind(cfos_highers, sum_score = sapply(seq(nrow(cfos_highers)), function(i){sum(cfos_highers[i,])}))
    
    rfp_highers <- sapply(seq(nrow(highers)), function(i){
      sapply(seq(ncol(highers)), function(j){
        ifelse(highers[i,j] == -1 & p_array[i,j,1] < 0.05,1,0)
      })
    })
    rfp_highers <- t(rfp_highers)
    colnames(rfp_highers) <- colnames(highers)
    rownames(rfp_highers) <- rownames(highers)
    rfp_highers <- cbind(rfp_highers, sum_score = sapply(seq(nrow(rfp_highers)), function(i){sum(rfp_highers[i,])}))
    
  }
  
  
  return(list(cfos = cfos_highers,
              rfp = rfp_highers))
}




