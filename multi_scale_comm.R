if (!('aricode' %in% installed.packages()[,'Package'])){install.packages("aricode")}; require(aricode);
source('/home/cesar/Dropbox/R/Network sufficiency/Modular_codes/consensus_community.R');

multi_scale_comm <- function(g, ALG = signed_louvain, gamma = seq(-1, 2, 0.05), R = 1, seed = 87, param_name = 'gamma', sec_param_name = NULL, ...){
  #
  # This function applis the community detection algorithm provided in FUN , R times for each gamma
  # value provided.
  #
  # INPUT: 
  #       g,              igraph object or dataframe structured as detailed in the 
  #       FUN,            Community detection algorithm. It accepts four functions for now
  #                       signed_louvain(default, adapted from Brain connectivity toolbox for matLab), 
  #                       cluster_resolution (resolution pckg), Louvain algorithm, but unoptimized and does not accept signed graphs)
  #                       spinglass.community (igraph pckg)
  #                       modularity maximization (adapted from Brain connectivity toolbox for matLab)
  #       gamma,          parameter resolution (gamma in all functions but cluster_resolution, which is called t)
  #       sec_param_name  
  #
  # OUTPUT: List(
  #         comm,   community affiliation vectors
  #         gamma,  resolution parameter value for each community affiliation
  #
  
  
  gamma <- rep(gamma, each = R)
  set.seed(seed)
  comm <- lapply(gamma, function(i){
    if(is.null(sec_param_name)){
      params <- list(g, i, ...)
      names(params)[2] <- param_name
      ALG(g, i, ...)
    }
    else{
      params <- list(g, i, ...)
      names(params)[2] <- param_name
      names(params)[3] <- sec_param_name
      ALG(g, i, i,...)
    }
  })
  
  #gamma <- sapply(seq_along(gamma), function(i){
  #  paste(i, gamma[i], sep = '_')
  #})
  
  names(comm) <- gamma
  
  return(comm)
}

ami_comps <- function(partitions, clusters, nets){
  #
  #
  output <- lapply(seq(partitions), function(i){ 
      clust_comp(partitions[[i]], clusters[[i]]$membership, V(nets[[i]])$anatomical_groups) 
    })
  names(output) <- names(partitions)
  return(output)
}

clust_comp <- function(Ci, memb, anato){
  
  # INPUT:
  #
  #
  #
  # OUTPUT:
  #    gamma      gamma value of the corresponding cluster classification
  #    ami        vector with the Adjusted mutual information of each sample in Ci
  
  memb <- cbind(as.numeric(as.factor(memb)), as.numeric(as.factor(anato)))
  gamma = Ci$gammas; 
  b = Ci$b;
  b_events = Ci$b_events;
  g_events = Ci$g_events
  Ci = Ci$Ci
  
  d <- sapply(seq(ncol(memb)), function(j){
    sapply(seq(ncol(Ci)), function(i){
      ami = aricode::AMI(memb[,j], Ci[,i])
    })
  })
  d <- data.frame(gamma = gamma, b = b, HCC = d[,1], Anatomical = d[,2])
  return(d)
}

group_ami <- function(clusters, hierarchical = F, co_classification = F){
  
  stopifnot(is.logical(hierarchical))
  
  #All pair wise comparisons
  lab <- combn(names(clusters), m=2, simplify = FALSE)
  
  if(isFALSE(co_classification)){
    if(isTRUE(hierarchical) & 
       all(sapply(clusters, is.hierarchical))){
      #Get levels of clustering
      Ci <- lapply(clusters, function(i){ i$merges })
      s <- max(sapply(Ci, ncol))
      #make the Cis the same number of columns (levels)
      Ci <- lapply(Ci, function(i){
        if(ncol(i) < s){
          i <- cbind(i, rep(i[,ncol(i)], (s - ncol(i))))
        }
        else{
          i
        }
      })
      
      #calculate ami
      gami <- data.frame(lapply(lab, function(i){
        sapply(2:s, function(j){
          aricode::AMI(as.numeric(Ci[[i[[1]]]][,j]), 
                       as.numeric(Ci[[i[[2]]]][,j]))
        })
      }))
    }
    else{ #if not hierarchical
      #Get cluster allegiances
      Ci <- lapply(clusters, function(i){ i$membership })
      
      # calculate ami
      gami <- data.frame(lapply(lab, function(i){ 
        aricode::AMI(as.numeric(Ci[[i[[1]]]]), 
                     as.numeric(Ci[[i[[2]]]]))
      }))
    } 
  }
  else{
    if(isTRUE(hierarchical) & 
       all(sapply(clusters, is.hierarchical))){
      #Get levels of clustering
      Ci <- lapply(clusters, function(i){ i$merges })
      s <- max(sapply(Ci, ncol))
      #make the Cis the same number of columns (levels)
      Ci <- lapply(Ci, function(i){
        if(ncol(i) < s){
          i <- cbind(i, rep(i[,ncol(i)], (s - ncol(i))))
        }
        else{
          i
        }
      })
      
      #obtain co-classification matrix for each level of each cluster
      Ci <- lapply(Ci, function(i){
        lapply(2:ncol(i), function(j){
          outer(i[,j], i[,j], '==')
        })
      })
      
      #calculate ami
      gami <- data.frame(lapply(lab, function(i){
        sapply(2:s, function(j){
          aricode::AMI(as.numeric(Ci[[i[[1]]]][[j]]), 
                       as.numeric(Ci[[i[[2]]]][[j]]))
        })
      }))
    }
    else{
      #Get cluster allegiances
      Ci <- lapply(clusters, function(i){ i$membership })
      
      #calculate co_classification
      co_classification <- lapply(clusters, function(i){
        outer(i$membership, i$membership, '==')
      })
      
      # calculate ami
      gami <- lapply(lab, function(i){ 
        aricode::AMI(as.numeric(Ci[[i[[1]]]]), 
                     as.numeric(Ci[[i[[2]]]]))
      })
    }
  }
 
  
  names(gami) <- sapply(seq(length(lab)), function(x){
    paste(lab[[x]], collapse = " ~ ")
  })
  
  return(gami)
}

afami <- function(nets, clusters){
  af <- lapply(names(nets), function(i){ 
    aricode::AMI(clusters[[i]]$membership, 
                 V(nets[[i]])$anatomical_groups) })
    names(af) <- names(nets)
    return(af)
}
