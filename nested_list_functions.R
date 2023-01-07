
multi_net_measures <- function(G_list, sec_param, FUN = cluster_louvain, list_levels = 1, add_param = NULL, 
                               different_sizes = F, ...){
  #
  # INPUT:
  #     G_list       A list, or a list of lists, of igraph objects.
  #     sec_param    Optional (NULL, default). Some functions will require a minimum of two arguments 
  #                  (ex. Consensus_cluster). Since this is a generic function and accepts many functions,
  #                  when this argument has an input, the name of the parameter will have to be added
  #                  to the add_param argument.
  #     FUN          a function to recursively be applied to all the graphs in the list.
  #     list_levels  Either 1, for a list of igraph objects, or 2, for a list of lists of igraph objects.
  #     add_param    Optional (NULL as default). if there is a second parameter argument, the name of this
  #                  argument should be put in add_para under quotes (ex. sec_param = clusters, add_param = 'comm'). 
  #
  #
  # OUTPUT:
  #    List, or list of lists, of the output of the function in FUN.
  
  if (list_levels == 1){
    
    output <- list()
    if(is.null(add_param)){
      output <- lapply(G_list, function(j){ 
        FUN(j, ...) 
      })
    }
    else{
      output <- list()
      output <- lapply(seq(G_list), function(i){
        params <- list(G_list[[i]], sec_param[[i]], ...)
        names(params)[2] <- add_param
        output <- do.call(FUN, params)
      })
      names(output) <- names(G_list)
    }
  }
  else{
    if(list_levels == 2){
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
        names(output) <- names(G_list)
        for(i in seq(G_list)){
          names(output[[i]]) <- names(G_list[[i]])
        }
      }
    }
  }
  
  if(isTRUE(different_sizes)){
    output <- nodes_match(output, list_levels = list_levels)
  }
  
  return(output)
}

nodes_match <- function(G_list, list_levels = 2){
  #
  #
  #
  #
  # OBS: when using multi_net_measures() on lists of networks, networks of different
  # sizes may need to be matched in their number and ID of nodes. This function gets the 
  # largest network and adds nodes to the smaller ones with the same name as in the largest
  # and gives value of zero to it.
  stopifnot(list_levels <= 2 | list_levels >= 1)
  
  if(list_levels == 1){
    
    #gets number of nodes
    msizes <- unlist(lapply(G_list, function(i){
      nrow(as.data.frame(i))
    }))
    
    # identifies largest network
    nodes <- rownames(as.data.frame(G_list[which(msizes %in% max(msizes))]))
    # identifies missing nodes in each network
    lacking_nodes <- lapply(G_list, function(i){
      nodes[!nodes %in% rownames(as.data.frame(i))]
    })
    
    #create subs with value 0 for each column present
    replacing_nodes <- lapply(seq(lacking_nodes), function(i){
      matrix(0, nrow = length(lacking_nodes[[i]]), 
             ncol = ncol(as.data.frame(G_list[[i]])))
    })
    # Give names to replacing nodes
    replacing_nodes <- lapply(seq(replacing_nodes), function(i){ 
      rownames(replacing_nodes[[i]]) <- c(lacking_nodes[[i]])
      colnames(replacing_nodes[[i]]) <- names(G_list[[i]])
      return(replacing_nodes[[i]])
    })
    
    #merging all
    output <- lapply(seq(G_list), function(i){
      rbind(as.data.frame(G_list[[i]]),
            as.data.frame(replacing_nodes[[i]]))
    })
    
    # Putting names
    names(output) <- names(G_list)
    
    #alligning node names across lists(graphs)
    output <- lapply(output, function(i){
      i[order(match(rownames(i),nodes)),]
    })
    
    if(is.list(G_list[[1]])){
      output <- lapply(output, function(i){
        lapply(i, function(x){
          x[x != ""]
          names(x) <- nodes
          return(x)
        })
      })
    }
  }
  else{
    if(list_levels == 2){
      
      #gets number of nodes
      msizes <- unlist(lapply(G_list, function(i){
        lapply(i, function(j){
          nrow(as.data.frame(j))
        })
      }))
      
      # identifies largest network
      nodes <- rownames(as.data.frame(G_list[which(msizes %in% max(msizes))]))
      # identifies missing nodes in each network
      lacking_nodes <- lapply(G_list, function(i){
        lapply(i, function(j){
          nodes[!nodes %in% rownames(as.data.frame(j))]
        })
      })
      
      #create subs with value 0 for each column present
      replacing_nodes <- lapply(seq(lacking_nodes), function(i){
        lapply(seq(lacking_nodes[[i]]), function(j){ 
          matrix(0, nrow = length(lacking_nodes[[i]][[j]]), 
                 ncol = ncol(as.data.frame(G_list[[i]][[j]])))
        })
      })
      # Give names to replacing nodes
      replacing_nodes <- lapply(seq(replacing_nodes), function(i){ 
        lapply(seq(replacing_nodes[[i]]), function(j){ 
          rownames(replacing_nodes[[i]][[j]]) <- c(lacking_nodes[[i]][[j]])
          colnames(replacing_nodes[[i]][[j]]) <- names(G_list[[i]][[j]])
          return(replacing_nodes[[i]][[j]])
        })
      })
      
      #merging all
      output <- lapply(seq(G_list), function(i){
        lapply(seq(G_list[[i]]), function(j){
          rbind(as.data.frame(G_list[[i]][[j]]),
                as.data.frame(replacing_nodes[[i]][[j]]))
        })
      })
      
      # Putting names
      output <- lapply(seq(G_list), function(i){
        names(output[[i]]) <- names(G_list[[i]]);
        return(output[[i]])
      })
      names(output) <- names(G_list)
      
      #alligning node names across lists(graphs)
      output <- lapply(output, function(i){
        lapply(i, function(j){
          j[order(match(rownames(j),nodes)),]
        })
      })
      
      if(is.list(G_list[[1]][[1]])){
        output <- lapply(output, function(i){
          lapply(i, function(j){
            lapply(j, function(x){
              x[x != ""]
              names(x) <- nodes
              return(x)
            })
          })
        })
      }
    } 
  }
  return(output)
}

matrix_list <- function(nets, list_levels = 1){
  if(list_levels == 1){
    nets <- lapply(nets, function(j){
      as.matrix(j[])
    }) 
  }
  else{
    if(list_levels == 2){
      nets <- lapply(nets, function(j){
        lapply(j, function(i){
          as.matrix(i[])
        })
      }) 
    }
  }
  return(nets)
}

split.along.dim <- function(a, n = length(dim(a)), output = c('list', 'dataframe')){
  #
  # this function splits every dimension of an array into elements of a list, 
  # or list of lists.
  #
  
  if(output == 'list'){
    s <- setNames(lapply(split(a, arrayInd(seq_along(a), dim(a))[, n]),
                         array, dim = dim(a)[-n], dimnames(a)[-n]),
                  dimnames(a)[[n]])
    return(s)
  }
  else{
    if(output == 'dataframe'){
      s <- sapply(split(a, arrayInd(seq_along(a), dim(a))[, n]),
                  array, dim = dim(a)[-n], dimnames(a)[-n])
      s <- data.frame(t(s))
      #colnames(s) <- dimnames(nets_generals_comp$diffs)$dim2
      #s <- data.frame(comparions = dimnames(nets_generals_comp$diffs)$dim3, s)
      return(s)
    }
  }
}
