
significant_comps <- function(metrics_list = NULL, comp_list, alpha = 0.05, trim = FALSE,
                              coi = NULL){
  
  # metrics_list    List of network metrics
  # comp_list       output of the cbn_permutation() or nodal_permutation_test() functions.
  #                 this output has to be from the same data input in metrics_list().
  # alpha           confidence level at which to discriminate significance.
  # trim            Logical (default, TRUE), whether to return only the significant values
  #                 or to return the significant comparisons as 1, in a separate binary vector
  # coi             Comparisons Of Interest (coi) is a Char vector. As the comparison functions
  #                 will make all possible pairwise comparisons, sometimes only some are of interest.
  #                 coi should be a char var as follows coi = c('group1 ~ group2', 'group2 ~ group3'),
  #                 or whatever are the names in ns var below or names(sig_h).

  stopifnot(is.logical(trim))

  if(!is.null(metrics_list)){
    h <- higher_metric(metrics_list)  
  }
  else{
    h <- higher_metric2(comp_list)
  }
  
  sig_h <- lapply(seq(ncol(h)), function(x){
    data.frame(comparison = rep(colnames(h)[x], nrow(h)),
      region = rownames(h),
      diff = as.numeric(comp_list$statistic[,x]),
      p.value = as.numeric(comp_list$pvalue[,x]), 
      higher_group = h[,x])
  })
  
  if(isTRUE(trim)){
    sig_h <- lapply(sig_h, function(x){
      x[x$p.value <= alpha,]
    }) 
  }
  else{
    sig_h <- lapply(sig_h, function(i){
      sig <- rep(0, nrow(i))
      sig[i$p.value <= alpha] <- 1
      i <- cbind(i, significant = sig)
    })
  }
  ns <- sapply(sig_h, function(i){ i$comparison[1] })
  names(sig_h) <- ns
  
  if(!is.null(coi)){
    sig_h <- sig_h[grepl(paste(coi, collapse = '|'), names(sig_h))] 
  }
  
  sig_h <- do.call(rbind.data.frame, sig_h)
  s <- t(as.data.frame(strsplit(sig_h$comparison, '\\.')))
  if(ncol(s)>2){
    s <- data.frame(comparison = factor(sapply(seq(nrow(s)), function(i){
      paste(s[i,1:(ncol(s)-1)], collapse = '.')
    })),
    metric = factor(s[,ncol(s)], levels = unique(s[,ncol(s)])))
  }
  else{
    s <- data.frame(comparison = s[,1], metric = s[,2], stringsAsFactors = T)
  }
  
  s <- data.frame(s, sig_h[,-1], stringsAsFactors = F)
  rownames(s) <- seq(nrow(s))
  
  return(s)
}

higher_metric <- function(metrics_list){
  #aux function for significant_comps()
  if(all(sapply(metrics_list, is.list))){
    metrics_list <- sapply(metrics_list, as.data.frame, simplify = F)
  }
  
  diff_label <- combn(names(metrics_list), m=2, FUN = paste0)
  emp_diffs <- list()
  comb <- combn(metrics_list,2, simplify = FALSE)
  emp_diffs <- sapply(comb, function(x){ 
    as.data.frame(x[[1]] - x[[2]])
  }, simplify = FALSE)
  
  h <- sapply(seq(emp_diffs), function(i){
    a <- emp_diffs[[i]] > 0
    a[a == 'TRUE'] <- diff_label[1,i]
    a[a == 'FALSE'] <- diff_label[2,i]
    return(a)
  }, simplify = F)
  
  names(h) <- sapply(seq(ncol(diff_label)), function(x){ 
    paste(diff_label[,x], collapse = " ~ ")
  })
  
  h <- cbind.data.frame(h)
  
  return(h)
}

higher_metric2 <- function(comp_list){
  #aux function for significant_comps()
  
  s <- t(as.data.frame(strsplit(colnames(comp_list$statistic), '\\.')))
  if(ncol(s)>2){
    s <- sapply(seq(nrow(s)), function(i){
      paste(s[i,1:(ncol(s)-1)], collapse = '.')
    })
  }
  else{
    s <- s[,1:(ncol(s)-1)]
  }
  s <- t(as.data.frame(strsplit(s, ' ~ ')))
  
  h <- sapply(seq(ncol(comp_list$statistic)), function(i){
    a <- comp_list$statistic[,i] > 0
    a[a == TRUE] <- s[i,1]
    a[a == FALSE] <- s[i,2]
    return(a)
  })
  colnames(h) <- colnames(comp_list$statistic)
  rownames(h) <- rownames(comp_list$statistic)
  return(as.data.frame(h))
}


comparison_hubs <- function(sig_cbn){
  
  sig <- subset(sig_cbn, select = -c(diff, p.value))
  
  # get all groups involved
  grupos <- unique(unlist(strsplit(as.character(sig$comparison), ' ~ ')))
  
  #for each group
  hubs <- lapply(seq(grupos), function(i){
    # get all their comparisons
    s <- droplevels(sig[grepl(grupos[i], sig$comparison),])
    
    # b and g will have all the significant and the higher_group infos respectively
    b <- g <- data.frame(region = sig[sig$comparison == levels(sig$comparison)[1], 'region'],
                         metric = sig[sig$comparison == levels(sig$comparison)[1], 'metric'],
                         matrix(0, nrow = length(unique(s$region))*length(levels(s$metric)), ncol = length(levels(s$comparison))))

    for(j in seq(levels(s$comparison))){ 
      b[,j+2] <- s[s$comparison == levels(s$comparison)[j], 'significant']
      a <- unlist(strsplit(levels(s$comparison)[j], ' ~ '))
      a <- a[a!=grupos[i]]
      colnames(b)[j+2] <- a
      g[,j+2] <- s[s$comparison == levels(s$comparison)[j], 'higher_group']
      colnames(g)[j+2] <- a
        }
    # values of b for grupo[i] significantly smaller than the compared group is set to -1
    for(j in 3:ncol(b)){
      b[ b[,j]==1 & g[,j] == colnames(b)[j],j] <- -1
    }
    b$hubs <- sapply(seq(nrow(b)), function(j){ sum(b[j,3:ncol(b)]) })
    return(b)
  })
  
  names(hubs) <- grupos
  return(hubs)
}


sig_func_effect <- function(nets, func_effect, p_value = 0.05){
  
  FF <- func_effect
  #FF <- results$func_effect
  
  m <- lapply(nets, function(i){ as.matrix(i[])})
  
  labels <- combn(names(m), 2, FUN = paste0)
  comb <- combn(m,2,simplify = F)
  diffs <- lapply(comb, function(x){ 
    a <- x[[1]]
    a[] = 1                 #makes a matrix of 1s
    diag(a) <- 0
    a[x[[1]] < x[[2]]] <- -1 # and turn to 2s all correlations where second name is higher
    return(a)
  })
  names(diffs) <- sapply(seq(ncol(labels)), function(x){ 
    paste(labels[,x], collapse = " ~ ")
  })
  
  sig_diffs <- lapply(names(diffs), function(i){
    x <- diffs[[i]]
    x[FF[[i]]$p_value > p_value] <- 0 # here you can put a different p_value
    return(x)
  })
  #insert an fdr_correction here
  names(sig_diffs) <- names(diffs)
  return(sig_diffs)
}



