source('/home/cesar/Dropbox/R/Network sufficiency/Modular_codes/nested_list_functions.R')
source('/home/cesar/Dropbox/R/Network sufficiency/Modular_codes/consensus_community.R')
source('/home/cesar/Dropbox/R/Network sufficiency/Modular_codes/signed_louvain.R')
source('/home/cesar/Dropbox/R/Network sufficiency/Modular_codes/modularity_functions.R');
source('/home/cesar/Dropbox/R/Network sufficiency/Modular_codes/multi_scale_comm.R');

#source('/home/cesar/Dropbox/R/Network sufficiency/Modular_codes/Figure_functions.R');

require(igraph)
require(doParallel)
require(Matrix)

HierConsenClust <- function(g = NULL, Ci, perm = localPermModel, distrib = sampleApprox, alpha = 0.05,
                            optimizer = signed_louvain,...){ 
  #
  #
  #
  #TODO
  # DESCRIPTION
  if(!is.null(g)){
    stopifnot(is.igraph(g) | is.matrix(g))
  }
  #The args for each function need to be lists
  #perm = localPermModel
  #distrib = normalApprox
  #alpha = 0.05
  #optimizer = signed_louvain
  if(is.list(Ci)){
    Ci = Ci$Ci
  }
  
  A <- agreement(Ci);                                                        # Co-classification Matrix
  N <- nrow(Ci);
  L <- ncol(Ci);
  Sc <- matrix(1, N, 1);
  coms_old <- 0;
  coms_new <- 0;
  coms <- 1;
  level <- 1;
  G <- matrix();
  Cit <- Sc;
  # Loop until no new comms are found
  while(coms_new != coms){
    
    G <- Matrix::sparseMatrix(1:N, Sc)
    #G <- 1*as.matrix(Matrix::sparseMatrix(1:N, Sc))
    coms_new <- coms;
    cat('Level:', level,'\n');
    cat('Spliting', coms_new - coms_old, 'new communities','\n');
    for(i in (coms_old+1):ncol(G)){
      ind <- G[,i]
      ind <- which(ind != 0)
      if(length(ind) > 1){
        A_it <- A[ind,ind];
        p <- distrib(perm(Ci[ind,]), alpha);
        if(any(rowSums(p) == 0)){
          cat('Null model has 0-valued entries.\n
          That means some regions did not vary their cluster allegiance. \n
          Condider increasing partition sample or reducing significance level! \n')
        }
        B <- A_it - p;
        
        while(TRUE){
          if(any(B <= 0)){
            Sc_it <- matrix(0, nrow(B), L)
            
            # paralelize this for
            Sc_it <- sapply(seq(L), function(i){
              Sc_it[,i] <- optimizer(B, mod = 'neg_asym', class=F)
            })
            
            A_it <- agreement(Sc_it)
            if(identical(A_it, 1*(A_it ==1))){
              Sc_it <- Sc_it[,1]
              break
            }
            else{
              B <- A_it - distrib(perm(Sc_it), alpha = alpha)
            }
          }
          else{
            Sc_it <- rep(1, nrow(B))
            break
          }
        }
        
        if(max(Sc_it) > 1){
          Sc[ind] <- Sc_it+coms
          c_it <- max(Sc_it)
          #Simtype code taken off
          #val <- mean(A[ind,ind])
          coms <- coms + c_it
        }
      }
    }
    Cit <- cbind(Cit, Sc)
    level <- level +1
    coms_old <- coms_new
  }
  
  rownames(Cit) <- rownames(Ci)
  # This block identifies the unique communities across
  # all levels of the hierarchy
  # it reorders the matrix so that we can see all degrees of agreement
  # each region is on
  # will need to put this transformation in the supplemental..
  e <- sapply(seq(ncol(Cit)), function(i){
    as.numeric(factor(do.call(paste, data.frame(Cit[,1:i], sep = ""))))
  })
  rownames(e) <- rownames(Cit)
  
  # This block identifies the final communities as
  # the last level before it became a singleton (isolated)
  # Maybe I just leave this as the Cit. Under evaluation.
  s <- numeric()
  for(i in 1:ncol(Cit)){
    s <- cbind(s, as.numeric(factor(do.call(paste, data.frame(Cit[,1:i], sep = "")))))
    if(any(table(s[,ncol(s)]) ==1)){
      t <- table(s[,ncol(s)]);
      t <- as.numeric(names(t)[t ==1]);
      s[s[,ncol(s)] %in% t, ncol(s)] <- s[s[,ncol(s)] %in% t, ncol(s)-1]
    }
  }
  s <- sapply(seq(ncol(s)), function(i){
    as.numeric(as.factor(s[,i]))
  })
  rownames(s) <- rownames(Cit)
  ### End of part to maybe remove
  
  if(is.null(g)){
    mt <- list(
      merges = s,
      membership = s[,ncol(s)],
      names = rownames(Cit),
      merges_Cit = Cit,
      merges_plot = e,
      algorithm = paste('Hierarchical Consensus Clustering')
    )
    mt <-  structure(mt, class = 'communities')
  }
  else{
    mt <- list(
      merges = s,
      membership = s[,ncol(s)],
      #modularities = sapply(seq(ncol(s)), function(i){ 
      #  Q_signed(g, s[,i]) 
      #}),
      modularity = Q_signed(g, s[,ncol(s)]),
      names = rownames(Cit),
      merges_Cit = Cit,
      merges_plot = e,
      algorithm = paste('Hierarchical Consensus Clustering')
    )
    mt <-  structure(mt, class = 'communities')
    
  }
  return(mt)
}

########################################################################################################################
# Testing the function
#g <- as.matrix(read.csv('/home/cesar/Documents/Matlab_toolbox_codes/A_matlab.csv', row.names = 1))
#g <- graph.adjacency(g, mode = 'undirected', weighted = T, diag = F)
#Ci <- eventSamples(g, 1000)[['Ci']]
#cluster <- HierConsenClust(g, Ci)

#set.seed(1987)
#g <- matrix(sample(c(rnorm(150, mean = 0.1, sd = 0.2),rnorm(250, mean=0.5, sd=0.2)), replace = F),20)
#rownames(g) <- 1:nrow(g)
#colnames(g) <- 1:nrow(g)
#g <- g + t(g);
#diag(g) <- 0
#g <- graph.adjacency(g, mode = 'undirected', weighted = T, diag = F)
#Ci <- eventSamples(g, 1000)[['Ci']]
#cluster <- HierConsenClust(g, Ci)

#is_hierarchical(cluster)
#plot_matrices(g, cluster, is.corr = F, anato_first = F, hier_clust = T)
#plot_matrices(comm = cluster, is.corr = F, anato_first = F)
