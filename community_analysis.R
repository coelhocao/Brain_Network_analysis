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
