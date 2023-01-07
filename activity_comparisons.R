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
  
if (!('dplyr' %in% installed.packages()[,'Package'])){install.packages("dplyr")}; require(dplyr)
if (!('boot' %in% installed.packages()[,'Package'])){install.packages("boot")}; require(boot)
if (!('multcomp' %in% installed.packages()[,'Package'])){install.packages("multcomp")}; require(multcomp)
if (!('heplots' %in% installed.packages()[,'Package'])){install.packages("heplots")}; require(heplots)
if (!('effsize' %in% installed.packages()[,'Package'])){install.packages("effsize")}; require(effsize)
if (!('Hmisc' %in% installed.packages()[,'Package'])){install.packages("Hmisc")}; require(Hmisc) #better functions for correlation matrices and p-values
if (!('lattice' %in% installed.packages()[,'Package'])){install.packages("lattice")}; require(lattice) #for generating graphs, colors and matrices
if (!('rcompanion' %in% installed.packages()[,'Package'])){install.packages("rcompanion")}; require(rcompanion) #for pseudo r-squared in glms
if (!('lmtest' %in% installed.packages()[,'Package'])){install.packages("lmtest")}; require(lmtest) # For Likelyhood ratio in glms
  
########################################################################################################################
###########  Activity Analysis of all the regions and factors  #########################################################
########################################################################################################################

scale_factor <- function(datalist, factor = 'Cohort', ave_across_brain = T){
  
  if(ave_across_brain){ #scales the factor across the whole brain
    datalist$counts <- ave(datalist$counts, datalist$labels$Cohort, FUN = scale)
  }else{
    
    for(i in 1:ncol(datalist$counts)){ # scales the factor within each brain region
      datalist$counts[,i] = ave(datalist$counts[,i], datalist$labels[factor], FUN = scale)
    }
  }
  
  return(datalist)
}


factor_boot <- function(df, factor = factor, i){
  # INPUT: df with a factor column 
  # OUTPUT: the t value of the comparison with resampled data
  # NOTE: the variable names used in the formula come from another function that calls this one.
  
  #resampling $counts list
  df$counts <- df$counts[i]
  #performing T test on resampled data
  t <- t.test(counts ~ hems, data = df, var.equal = F)$statistic
  # In case there are more than 2 groups here is a GLHT to test of the groups (they have to be in a single var though)
  return(t)
}

## p-value calculation
p.value <- function(booted, fixed = FALSE){
  # INPUT: output of a boot function
  # extracts the proportion of values from resampled distribution that's bigger than empirical value OR
  # bigger than a fixed value (default = FALSE)
  # OUTPUT: a data.frame with the empirical estimate, std error and p-value of the booted variable.
  emp <- booted$t0
  dist <- booted$t
  p <- numeric(length(emp))
  std.err <- numeric(length(emp))
  if(fixed == FALSE){
    for (i in 1: length(p)){
      p[i] <- mean(abs(dist[,i]) > abs(emp[i]))
      std.err[i] <- sd(dist[,i])
    }
  } else{ 
    if(is.numeric(fixed)){
      for (i in 1: length(p)){
        p[i] <- mean(abs(dist[,i]) > abs(fixed))
        #std.err[i] <- sd(dist[,i])
        std.err[i] <- summary(booted$bootSE)
      }
    } else{
      stop('fixed should be FALSE or NUMERIC!!!')
    }
  }
  results <- rbind(emp, std.err, p)
  return(results)
}

CI <- function(booted){
  
  cis <- sapply(1:length(booted$t0), function(x){
    boot.ci(booted, conf = 0.95, type = 'basic', index = c(x,x))$basic[4:5] ##'bca' will only work if p >0
  })
  cis <- rbind(cis)
  rownames(cis) <- c('lower', 'upper')
  colnames(cis) <- names(booted$t0)
  return(cis)
}

##Effect sizes for multiple factors
eff_size <- function(dd, factors = c('group', 'day'), quant = 'correct', type = c('one-sample', 'two-sample'), H0 = NULL){ 
  # INPUT: data and a list of comparisons 1x1 to be made 
  #OUTPUT: cohen's d effect size of every comparison made
  #NOTE: the VARs have hierarchy, such that that f1 is inside f2 which is inside f3
  #EX: Comparison of behav response btwn 2 contexts (behav ~ context) within each day (class2) and each group (class1)
  
  # New Var interacting Vars to included in the multiple comparison
  dd <- droplevels(dd)
  if(length(factors) > 1) {
    dd$interac <- interaction(dd[,match(factors, names(dd))])
  }
  else{
    if(length(factors) == 1){
      dd$interac <- dd[,match(factors, names(dd))]
    }
    else{
      stop('Please provide 1 or more factors to have the levels compared serialy!')
    }
  }
  
  
  if(type == 'two-sample' & is.null(H0)){
    e <- data.frame(c(sapply(1:(length(levels(dd$interac))-1), function(x){
      c(sapply(2:length(levels(dd$interac)), function (y){
        cohen.d(dd[dd[,'interac'] == levels(dd[,'interac'])[x], quant], dd[dd[,'interac'] == levels(dd[,'interac'])[y], quant])$estimate
      }))
    })))  
    
    rownames(e) <- c(sapply(1:(length(levels(dd$interac))-1), function(x){
      c(sapply(2:length(levels(dd$interac)), function (y){
        names(e) <- paste(levels(dd[,'interac'])[x], levels(dd[,'interac'])[y], sep = ' - ')
      }))
    }))  
    colnames(e) <- c('Cohen`s D')
  }
  else{
    if(type == 'one-sample' & !is.null(H0)){
      e <- data.frame(c(sapply(1:length(levels(dd$interac)), function (y){
        e <- abs((H0 - mean(dd[dd[,'interac'] == levels(dd[,'interac'])[y],quant]))/sd(dd[dd[,'interac'] == levels(dd[,'interac'])[y],quant]))
      })))
      
      rownames(e) <- c(sapply(1:length(levels(dd$interac)), function (y){
        names(e) <- levels(dd[,'interac'])[y]
      }))
      colnames(e) <- c('Cohen`s D')
    }
    else{
      stop('Error: for two-sample test, H0 must be zero. For one-sample tests H0 must be numeric')
    }
  }
  
  return(e)
}

factor_comp <- function(datalist, factor = 'hemisphere', isolated_factors = c('Exp', 'group', 'signal'), fdr_adjusted = F,
                        boot = F, R = 1000, seed = 1987, sim="balanced"){
  # INPUT: data_list as provided by mfiles_to_hierarchical_data() function, factor to be tested, 
  #and factors to be separated (not compared but not directly tested)
  #
  # OUTPUT: df with the isolated variables, regions, t and p values
  # NOTE: The factor and isolated_factors attributes should be factors in the bd$labels list
  # NOTE: The attributes R and sim are from the boot() and only matter if boot = T
  
  ld <- list()
  if(any(is.null(isolated_factors), is.na(isolated_factors))){
    ld[[1]] <- datalist
    names(ld) <- factor
  }else{
    datalist$labels['interac'] <- factor(with(datalist$labels, interaction(datalist$labels[,match(isolated_factors, names(datalist$labels))])))
    #separate the factors to be analyzed separately in lists
    
    ld <- lapply(seq_along(levels(datalist$labels$interac)), function(i){
      ld[[levels(datalist$labels$interac)[i]]] <- list(
        labels = data.frame(datalist$labels[datalist$labels$interac == levels(datalist$labels$interac)[i],]),
        counts = data.frame(datalist$counts[datalist$labels$interac == levels(datalist$labels$interac)[i],]))
    })
    names(ld) <- levels(datalist$labels$interac)
  }
  
  if(boot == F){
    R = NULL
    
    #performs t test in each region of each separate interac factor
    t <- lapply(seq_along(ld), function(i){
      data.frame(sapply(seq_along(ld[[i]]$counts), function(j){
        t <- t.test(counts ~ hems, data = cbind(hems = ld[[i]]$labels[,factor], counts = ld[[i]]$counts[,j]), var.equal = F)
      }))
    })
    
    
    results <- lapply(seq_along(t), function(i){
      data.frame(sapply(seq_along(t[[i]]), function(j){
        
        results <- cbind(t_value = t[[i]][[j]]$statistic, CI_lower = t[[i]][[j]]$conf.int[1], CI_upper = t[[i]][[j]]$conf.int[2], 
                         p_value = t[[i]][[j]]$p.value)
      }))
    })
    names(results) <- levels(datalist$labels$interac)
    for (i in seq_along(results)){
      colnames(results[[i]]) <- names(ld[[1]]$counts);
      rownames(results[[i]]) <- c('t_value', 'lower', 'upper', 'p_value')
    }
    
  }
  else{
    if(boot == TRUE){
      if(!is.numeric(R)){
        stop('R must be numeric! Number of resamplings!')
      }
      else{
        
        #Resampled data stats ditribution
        system.time(booted <- lapply(seq_along(ld), function(i){
          data.frame(sapply(seq_along(ld[[i]]$counts), function(j){
            set.seed(seed);
            boot(data = data.frame(hems = ld[[i]]$labels[,factor], counts = ld[[i]]$counts[,j]), statistic = factor_boot, R = R, sim = sim, factor = factor)
          }))
        }))
        names(booted) <- levels(datalist$labels$interac)
        for (i in seq_along(booted)){names(booted[[i]]) <- names(ld[[1]]$counts)}
        
        ## gather the results
        # get p-values
        ps <- list()
        ps <- lapply(seq_along(booted), function(i){
          data.frame(sapply(seq_along(booted[[i]]), function(j){
            ps <- p.value(booted[[i]][[j]])
          }))
        })
        names(ps) <- levels(datalist$labels$interac)
        for (i in seq_along(ps)){
          colnames(ps[[i]]) <- names(ld[[1]]$counts);
          rownames(ps[[i]]) <- c('stat', 'std_error', 'p_value')
        }
        
        #get CIs
        cis <- list()
        cis <- lapply(seq_along(booted), function(i){
          data.frame(sapply(seq_along(booted[[i]]), function(j){
            cis <- ifelse(sum(booted[[i]][[j]]$data$counts) > 0, data.frame(CI(booted[[i]][[j]])), data.frame(NaN, NaN))
          }))
        })
        names(cis) <- levels(datalist$labels$interac)
        for (i in seq_along(cis)){
          colnames(cis[[i]]) <- names(ld[[1]]$counts);
          rownames(cis[[i]]) <- c('lower', 'upper')
        }
        
        #Compiling results
        results <- lapply(seq_along(ps), function(i){
          data.frame(sapply(seq_along(ps[[i]]), function(j){
            results <- cbind(t_value = ps[[i]]['stat',j], std_error = ps[[i]]['std_error',j] , CI_lower = cis[[i]]['lower',j], 
                             CI_upper = cis[[i]]['upper',j], p_value = ps[[i]]['p_value',j])
          }))
        })
        names(results) <- levels(datalist$labels$interac)
        for (i in seq_along(results)){
          colnames(results[[i]]) <- names(ld[[1]]$counts);
          rownames(results[[i]]) <- c('t_value', 'std_error', 'lower', 'upper', 'p_value')
        }
        
      }
    }
  }
  
  #get effect_sizes
  #eff <- lapply(seq_along(ld), function(i){
  #  data.frame(sapply(seq_along(ld[[i]]$counts), function(j){
  #    eff <- eff_size(data.frame(hems = ld[[i]]$labels[,factor], counts = ld[[i]]$counts[,j]), 
  #                    factors = factor, quant = 'counts', type = 'two-sample', H0 = NULL)
  #  }))
  #})
  #names(eff) <- levels(datalist$labels$interac)
  #for (i in seq_along(eff)){
  #  colnames(eff[[i]]) <- names(ld[[1]]$counts);
  #  rownames(eff[[i]]) <- "Cohens_D"
  #}
  
  # FDR p_value adjustment if asked for
  if(fdr_adjusted == T){ 
    ps <- list()
    ps <- lapply(seq_along(results), function(i){
      ps[[i]] <- p.adjust(results[[i]]['p_value',], method = 'fdr')
    })
    #Joining effect_sizes AND adjusted p_values
    results <- lapply(seq_along(results), function(i){
      results[[i]] <- data.frame(rbind(results[[i]], adj_p_value = ps[[i]]))#, eff[[i]]))
    })
  }
  else{
    #Joining effect_sizes on results
    results <- lapply(seq_along(results), function(i){
      results[[i]] <- data.frame(rbind(results[[i]]))#, eff[[i]]))
    })
  }
  names(results) <- levels(datalist$labels$interac)
  for (i in seq_along(results)){ results[[i]] <- t(results[[i]]) }
  
  #OUTPUT
  return(as.data.frame(results))
}

level_vs_corr <- function(df){
  #INPUT: dataframe
  #
  #OUTPUT: a df with mean signal level, mean r^2 of signal, and the coefficient of variation of signal
  #OBS: This is used to verify whether the correlations coefficients are correlated with the signal level or to its variance
  
  mean_signal <- sapply(1:ncol(df), function(i){
    meanr <- mean(df[,i], na.rm = T)
  })
  
  m <- corr_matrix(df)$corr
  mean_r2 <- sapply(1:ncol(m), function(i){
    mean_r2 <- mean(m[,i]^2, na.rm = T)
  })
  
  co.var <- sapply(1:ncol(df), function(i){
    co.var <- (sd(df[,i],na.rm=TRUE)/mean(df[,i],na.rm=TRUE))
  })
  
  mean_r2[mean_r2 == 1] <- NA
  
  d <- data.frame(mean_r2, mean_signal, co.var)
  d <- d[!is.na(d$mean_r2) & !is.nan(d$mean_r2),]
  
  return(d)
}


model_boot <- function(data, i, mod = c('lm', 'glm')){
  
  #data = data[i,]
  data[,2] <- data[i,2]
  
  
  if(mod == 'lm'){
    model <- lm(data[,2] ~ data[,1])
    s <- c(summary(model)$coefficient[2], summary(model)$r.squared)
  }
  else{
    if(mod == 'glm'){
      model <- glm(data[,2] ~ data[,1])
      null <- glm(data[,2] ~ 1)
      r <- nagelkerke(model, null)[[2]][3]
      s <- c(summary(model)$coefficient[2], r)
    }
  }
  return(s)
}


activity_vs_connectivity <- function(datalist, quant = 'counts', factors = c('Exp', 'group', 'hemisphere', 'signal'), mod = c('lm', 'glm'), booted = F, R = 10000){
  
  datalist$labels['interac'] <- factor(with(datalist$labels, interaction(datalist$labels[,match(factors, names(datalist$labels))])))
  
  #separate the factors to be analyzed separately in lists
  ld <- list()
  ld <- lapply(seq_along(levels(datalist$labels$interac)), function(i){
    ld <- cbind(factor = rep(levels(datalist$labels$interac)[i], nrow(level_vs_corr(datalist$counts[datalist$labels$interac == levels(datalist$labels$interac)[i],]))),
                level_vs_corr(datalist[[quant]][datalist$labels$interac == levels(datalist$labels$interac)[i],]))
  })
  
  ldm <- data.frame()
  for(i in seq_along(ld)){
    ldm <- rbind(ldm, ld[[i]])
  }
  ldm$factor <- as.factor(ldm$factor)
  
  if(booted == F){
    if(mod == 'lm'){
      
      s <- list()
      s <- lapply(3:ncol(ldm), function(j){
        lapply(seq_along(levels(datalist$labels$interac)), function(i){
          s <- lm(ldm[ldm$factor == levels(ldm$factor)[i], 'mean_r2'] ~ ldm[ldm$factor == levels(ldm$factor)[i], j ])
        })
      })
      
      output <- list()
      output <- lapply(seq_along(s), function(j){
        data.frame(t(sapply(seq_along(levels(datalist$labels$interac)), function(i){
          
          output <- data.frame(
            sample = levels(datalist$labels$interac)[i],
            beta = summary(s[[j]][[i]])$coefficients[2,1],
            std_error = summary(s[[j]][[i]])$coefficients[2,2],
            t_value = summary(s[[j]][[i]])$coefficients[2,3],
            p_value = summary(s[[j]][[i]])$coefficients[2,4],
            upper.ci = confint(s[[j]][[i]])[2,1],
            lower.ci = confint(s[[j]][[i]])[2,2],
            R_squared = summary(s[[j]][[i]])$r.squared,
            F_value = summary(s[[j]][[i]])$fstatistic[1],
            model_p_value = pf(summary(s[[j]][[i]])$fstatistic[1], summary(s[[j]][[i]])$fstatistic[2], summary(s[[j]][[i]])$fstatistic[3], lower.tail = FALSE),
            stringsAsFactors = F)
        })))
      })
      names(output) <- c('mean_r2 ~ mean_signal','mean_r2 ~ coef_var')
    }
    else{
      if(mod == 'glm'){
        
        s <- list()
        s <- lapply(3:ncol(ldm), function(j){
          lapply(seq_along(levels(datalist$labels$interac)), function(i){  
            s <- glm(ldm[ldm$factor == levels(ldm$factor)[i], 'mean_r2'] ~ ldm[ldm$factor == levels(ldm$factor)[i], j ])
          })
        })
        
        null <- list()
        null <- lapply(3:ncol(ldm), function(j){
          null <- lapply(seq_along(levels(datalist$labels$interac)), function(i){
            null <- glm(ldm[ldm$factor == levels(ldm$factor)[i], 'mean_r2'] ~ 1)
          })
        })
        
        model <- list()
        model <- lapply(seq_along(s), function(j){
          lapply(seq_along(levels(datalist$labels$interac)), function(i){
          anova(s[[j]][[i]], null[[j]][[i]], test="Chisq") 
          })
        })
        
        output <- list()
        output <- lapply(seq_along(s), function(j){
        data.frame(t(sapply(seq_along(levels(datalist$labels$interac)), function(i){
          
          output <- data.frame(
            sample = levels(datalist$labels$interac)[i],
            beta = summary(s[[j]][[i]])$coefficients[2,1],
            std_error = summary(s[[j]][[i]])$coefficients[2,2],
            t_value = summary(s[[j]][[i]])$coefficients[2,3],
            p_value = summary(s[[j]][[i]])$coefficients[2,4],
            upper.ci = confint(s[[j]][[i]])[2,1],
            lower.ci = confint(s[[j]][[i]])[2,2],
            deviance = model[[j]][[i]]$Deviance[2],
            model_p_value = model[[j]][[i]]$`Pr(>Chi)`[2],
            pseudo_R_squared = nagelkerke(s[[j]][[i]],null[[j]][[i]])[[2]][3],
            likelyhood_ratio_Chiq = lrtest(s[[j]][[i]], null[[j]][[i]])$Chisq[2],
            likelyhood_ratio_p = lrtest(s[[j]][[i]], null[[j]][[i]])$`Pr(>Chisq)`[2],
            stringsAsFactors = F)
        })), stringsAsFactors = F)
        })
        names(output) <- c('mean_r2 ~ mean_signal','mean_r2 ~ coef_var')
      }
    }
  }
  else{
    if(booted == T){
      
      s <- list()
      s <- lapply(3:ncol(ldm), function(j){
        lapply(seq_along(levels(ldm$factor)), function(i){
        
        s <- boot(data = cbind(ldm[ldm$factor == levels(ldm$factor)[i], j ], ldm[ldm$factor == levels(ldm$factor)[i], 'mean_r2']), 
                       statistic = model_boot, R = R, mod = mod)
        })
      })
      
      ps <- list()
      ps <- lapply(seq_along(s), function(j){
        lapply(seq_along(levels(ldm$factor)), function(i){
        ps <- p.value(s[[j]][[i]])
        })
      })
      
      cis <- list()
      cis <- lapply(seq_along(s), function(j){
        lapply(seq_along(levels(ldm$factor)), function(i){
        cis <- CI(s[[j]][[i]])
        })
      })
      
      output <- list()
      output <- lapply(seq_along(s), function(j){
        data.frame(t(sapply(seq_along(levels(datalist$labels$interac)), function(i){
        
        output <- data.frame(
          sample = factor(levels(datalist$labels$interac)[i]),
          beta = ps[[j]][[i]][1,1],
          std_error = ps[[j]][[i]][2,1],
          p_value = ps[[j]][[i]][3,1],
          upper.ci = cis[[j]][[i]][1,1],
          lower.ci = cis[[j]][[i]][2,1],
          R_squared = ps[[j]][[i]][1,2],
          model_std_error = ps[[j]][[i]][2,2],
          model_p_value = ps[[j]][[i]][3,2],
          model_upper.ci = cis[[j]][[i]][1,2],
          model_lower.ci = cis[[j]][[i]][2,2],
          stringsAsFactors = F)
      })))
      })
      names(output) <- c('mean_r2 ~ mean_signal','mean_r2 ~ coef_var')
    }
    else{
      stop('only lm or glm are accepted as mod!')
    }
    if(mod == 'glm'){
      lapply(seq_along(s), function(j){
      names(output[[j]]) <- c('sample', 'beta', 'std_error', 'p_value', 'upper.ci', 'lower.ci', 'pseudo_R_squared', 'model_std_error', 
                    'likelyhood_ratio_p', 'model_upper.ci','model_lower.ci')
      })
    }
  }
  return(output)
}
