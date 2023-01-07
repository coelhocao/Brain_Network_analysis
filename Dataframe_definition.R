#############################################################################################
# These codes have been costumly designed for brain network analysis of 1-2 time point data
# generally acquired by immunolabelling of activity-dependent proteins (i.e. c-fos, Arc)
#
# Many functions here have been largely based on some of the codes shared by Justin Kenney,
# which can be found on this link https://github.com/jkenney9a/Networks
#
# NOTE: It is intended for undirected network data
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
# TODO
# In the trimmer() and remove_uncounted_child(), need to make the region/acronym fields neat across the functions
# Maybe a better fusion to deal with trimmer() and remove_uncounted_child()

if (!('dplyr' %in% installed.packages()[,'Package'])){install.packages("dplyr")}; require(dplyr)
if (!('Hmisc' %in% installed.packages()[,'Package'])){install.packages("Hmisc")}; require(Hmisc) #better functions for correlation matrices and p-values
if (!('data.table' %in% installed.packages()[,'Package'])){install.packages("data.table")}; require(data.table) #for rbindlist
if (!('data.tree' %in% installed.packages()[,'Package'])){install.packages("data.tree")}; require(data.tree) #for hierarchical organization of anatomical regions
if (!('lattice' %in% installed.packages()[,'Package'])){install.packages("lattice")}; require(lattice) #for generating graphs, colors and matrices
#################################################################
#IDs <- read.csv('C:/Users/CAOC/Dropbox/R/Network sufficiency/regions_IDs.csv', sep = ',', header = T, check.names = F, stringsAsFactors = F)
#require(rjson)
#Loads a file provided by Allen Brain Atlas with all the regions and their hierarchies
#abjson <- fromJSON(file = 'C:/Users/CAOC/Dropbox/R/Network sufficiency/allen_brian_atlas_hierarchy.json', simplify = T)
#hier_data <- FromListExplicit(abjson$msg[[1]][11], nameName = 'name') #this is the complete hierarchical data from Allen Brain Atlas
#################################################################

########################################################################################################################
###### Loading files and arraging Data frames ##########################################################################
########################################################################################################################

load_data <- function(file, header=FALSE, fill = T, balanced = T, from.ClearMap = T){
  #   Input: either a .txt or a .csv or a excel file name
  #   
  #   Output: Dataframe containing the file
  #   
  #   OBS: dataframes with missing data can't be loaded in the .txt format!!
  #   OBS2: data counts from ClearMap includes unbalanced columns with region names. For some reason, the last column,
  #   which has only a few lines with values, end up being put in the first column in the next raw.
  #   The balanced = F attribute deals with that. NOTE: simply not loading
  
  if(is.logical(from.ClearMap) & isTRUE(from.ClearMap)){
    if(balanced == T){
      if(tools::file_ext(file) == "txt"){
        d <- read.table(file=file, header=header, fill=fill, check.names=FALSE)[,1:2]
      }
      else{
        if(tools::file_ext(file) == "csv"){
          d <- read.csv(file=file, sep = ",", header=header, fill=fill, check.names=FALSE)[,1:2]
        }
      }
    }
    else{
      d <- readLines(file)
      d <- strsplit(d, ',')
      d <- lapply(d, function(x)as.data.frame(t(x), stringsAsFactors = F))
      d <- plyr::rbind.fill(d)
      d <- as.data.frame(d[,1:2])
      d[,1] <- as.numeric(d[,1])
      d[,2] <- as.numeric(d[,2])
    } 
  }
  else{
    if(is.logical(from.ClearMap) & isFALSE(from.ClearMap)){
      if(tools::file_ext(file) == "txt"){
        d <- read.table(file=file, header=header, fill=fill, check.names=FALSE)
      }
      else{
        if(tools::file_ext(file) == "csv"){
          d <- read.table(file=file, sep = ",", header=header, fill=fill, check.names=FALSE)
        }
      }
    }
  }
  return(d)
}

clean_data <- function(df, real_zeros=TRUE, missing_thresh= 0.5, fill_missing=TRUE){
  # Input: data as loaded
  # 
  #removes undesired zeros, remove columns with too many missing data, fill missing data with mean of column
  #
  # Output: cleaned data
  
  #if desired, Replace zeros in data with NA
  if (real_zeros==FALSE){
    df[mapply("==", df, 0)] <- NA
  }
  
  #Remove columns with more than missing_thresh missing values
  df_out <- df[,which(colSums(is.na(df))/dim(df)[1] < missing_thresh)]
  
  #Fills missing values with mean column value
  #(doesn't work if cleaned data doesn't have missing data to be filled)
  if (fill_missing==TRUE){
    col_avgs <- colMeans(df_out, na.rm=TRUE)
    index <- which(is.na(df_out), arr.ind=TRUE)
    if(!(length(index) > 1)){
      return(df_out)
      warning("Variables with missing data didn't pass the threshold and were excluded!")
    }
    else{
      df_out[index] <- col_avgs[index[,2]]
    }
  }
  return(df_out)
}


load_multiple_files <- function(IDpath, datapath, balanced = T, del_unlabeled = T){
  #   Input: path for the annotation file,a path where all the files are
  #   
  #   Output: one dataframe with all initial dataframes merged and separate into specified factors
  #   
  #   OBS: dataframes with missing data can't be loaded in the .txt format!!
  stopifnot(is.logical(balanced) |
              is.logical(del_unlabeled))
  
  temp <- list.files(path = datapath, pattern = '*Annotated_counts.csv')
  data_list <- lapply(temp, function(x){
    load_data(file= paste(datapath, x, sep="/"), header=FALSE, fill = T, balanced = balanced)
    })
  names(data_list) <- temp
  
  d <- read.csv(IDpath, sep = ',', header = T, check.names = F)
  b <- data.frame(sapply(1:length(data_list), function(y){
    a <- rbind(data_list[[y]], cbind(d[!d$id %in% data_list[[y]][,1],1], rep(0,length(d[!d$id %in% data_list[[y]][,1],1]))))
    b <- data.frame(a[match(d[,1],a[,1]),2])
  }))
  
  b <- data.frame(id = d[,1],parent_id = d[,8], region = as.character(d[,2]), acronym = gsub('-','.', as.character(d[,3])), 
                  parent_acronym = gsub('-', '.', as.character(d[,9])), b, stringsAsFactors = F)
  rownames(b) <- c(1:nrow(b))
  
  #Because I have not been able to identify these regions and I have to exclude them BEFORE summing up all levels so the dataframe marging works,
  # I am deleting it here on this step
  # from the six unlabelled regions This is what I'm almost sure of
  #182305696 - SSp-un  Primary Sensory area, unassigned
  #182305712 - SSp-un  Primary Sensory area, unassigned MAYBE NOT. BUT MUST CHECK
  #312782560 - Anterior area, under PTL
  #312782592 - VISrl1/3, under PTL  +  VISli, under VIS
  #312782624 - VISpor1/3, under VIS, + VISrl4/6, under PTLpVISli, under VIS
  #312782656 - VISpor3/4 , under VIS
  #
  if(del_unlabeled){
    b <- b[b[,'region'] != 'Unlabeled',]  
  }
  
  b <- list(annotation = b[,1:5], counts = b[,6:ncol(b)])
  
  return(b)
}

factors_from_names <- function(datapath){
  # INPUT: data files names
  # OUTPUT: factor dataframe with experiment (optional), subject, hemisphere, signal, group
  # NOTE: the file name has to have this info in the following order:
  # experiment (optional), subject, hemisphere, signal, group, separated by a "_".
  
  temp <- list.files(path = datapath, pattern = '*Annotated_counts.csv')
  temp <- gsub('_An.*', '',temp)
  labels <- t(data.frame(strsplit(temp, '_')))
  colnames(labels) <- c('experiment', 'subject', 'hemisphere', 'signal', 'group')

  return(data.frame(labels, stringsAsFactors = T))
}

hierarchical_sum <- function(df_tree){
  #INPUT: a dataframe with 4 columns: id, parent_id, region, signal (cell count)
  #
  #OUTPUT: Dataframe with the pathString of the datatree and the summed signal from children to root.
  #
  names(df_tree) <- c('id', 'parent', 'region', 'signal')
  
  tree <- FromDataFrameNetwork(df_tree[-1,])
  #tree$Do(function(node) node$levelName <- node$region)
  tree$Do(function(node) node$s <- ifelse(is.null(node$signal), 0, node$signal) + sum(Get(node$children, "s")), traversal = "post-order")
  tree <- ToDataFrameTree(tree,'name','s')
  return(tree)
}

drop_Layer_1 <- function(b, region = 'region'){
  # INPUT: dataframe and a var with region names
  #
  # OUTPUT: a vector with region names excluding the layer 1 from all cortices.
  
  tree <- FromDataFrameNetwork(b[-1,])
  # deleting fiber tracts, ventricular systems and grooves
  df <- ToDataFrameNetwork(tree$`8`, path = 'pathString', region = region)
  # Grab all names with the following elements
  d <- df[grepl('ayer 1', df$region),'region']
  d <- c(d, df[grepl('ayers 1', df$region),'region'])
  return(d)
}

remove_noncounted_child <- function(b){
  #
  # INPUT: DF with parents and childs IDs and counts data after hierarchical_sum()
  #
  # NOTE: This function was necessary because the segmentation script used by ClearMap
  # utilizes a level of segmentation different from the level in the complete data.tree
  # provided by the Allen Brain atlas. Ex., CEA is subdivided into CEAm, CEAl, CEAc, but
  # ClearMap quantifies them only at the level of CEA, and the regions above will have always counts = 0.
  # This is probably bc ClearMap uses an earlier script than the currently available.
  # Python users might be able to change/update the segmentation code on ClearMap
  # and not need this function.
  #
  # OPERATION: searches for the parents which sum of children counts = 0. Then, deletes its children.
  # this way, if any children was counted, the divisions remain. This serves for parents of parents.
  #
  # OUTPUT: DF on the higher (deeper) level counted by ClearMap
  
  #calculate sum of counts over subjs
  z <- sapply(seq_len(nrow(b[,7:ncol(b)])), function(x){
    z <- rowSums(b[x,7:ncol(b)])
  })
  z <- data.frame(acronym = factor(b$acronym),
                  parent = factor(b$parent_acronym), 
                  sum = as.numeric(z))
  
  # 
  e <- do.call(rbind.data.frame,
               sapply(levels(z$parent), FUN = function(x){
    if(sum(z[z$parent == x,'sum']) == 0){
      e <- z[z$parent == x,]
    }
  }))
  row.names(e) <- seq_len(nrow(e))

  z <- suppressMessages(anti_join(z,e))
  
  return(z$acronym)
}


trimmer <- function(b, z, region = 'region', output = 'df'){
  # Input: Data Tree with the attributes of interest (columns of subjects) summed from children to parents (with hierarchical_sum() function??)
  #
  # Output: A data frame trimmed without the Leafs not of interest.
  #
  # OBS: this function only trims children of the tree to the desired level. Does not trim parents
  # OBS2: The data tree to be inputed MUDT have a column (here called 'region') from which the names are taken
  
  # forms a tree with the dataframe inputed
  if(!('Node' %in% attr(b, 'class'))){
    tree_ <- FromDataFrameNetwork(b[-1,])
  }
  else{
    tree_ <- b
  }
  
  tree_$RemoveChild('1009') #Remove fiber tracts
  tree_$RemoveChild('1024') #Removes grooves
  tree_$RemoveChild('73') #Removes ventricular systems
  tree_$`8`$RemoveChild('512') #Removes cerebellum
  tree_$`8`$`343`$`1065`$RemoveChild('354') #Removes medulla
  
  #this function DOES NOT interpret anything that calls another value in the region for some reason
  df <- ToDataFrameNetwork(tree_$`8`, path = 'pathString', region = region) 
  
  # Moving Fields of Forel one level up so both zona incerta and Fields of Forel remain
  df[df$region == 'Fields of Forel', 'from'] <- 290
  df[df$region == 'Fields of Forel', 'path'] <- "997/8/343/1129/1097/290/804"
  
  
  # Grab all names with the following elements
  d <- df[grepl(',', df$region),'region']
  d <- c(d, df[grepl('/',df$region),'region'])
  d <- c(d, df[grepl('6',df$region),'region'])
  
  #Within the pool from the elements above, retain only the ones with the following elements
  a <- d[grepl('ayer', d)]
  a <- c(a, d[grepl('stratum', d)])
  a <- c(a, d[grepl('6', d)])
  
  #Grab parent and children with subnames, retain only children (to be deleted!)
  a <- c(a, df[grepl('Dentate gyrus',df$region),'region'])
  a <- a[a != "Dentate gyrus"]
  a <- c(a, df[grepl('stria terminalis',df$region),'region'])
  a <- a[a != "Bed nuclei of the stria terminalis"]
  a <- a[a != "Bed nuclei of the stria terminalis, anterior division"]
  a <- a[a != "Bed nuclei of the stria terminalis, posterior division"]
  a <- c(a, df[grepl('Mediodorsal nucleus of thalamus',df$region),'region'])
  a <- a[a != "Mediodorsal nucleus of thalamus"]  
  a <- c(a, df[grepl('Ventral part of the lateral geniculate complex',df$region),'region'])
  a <- a[a != "Ventral part of the lateral geniculate complex"]
  a <- c(a, df[grepl('Nucleus circularis',df$region),'region'])
  a <- c(a, df[grepl('Paraventricular hypothalamic nucleus',df$region),'region'])
  a <- a[a != "Paraventricular hypothalamic nucleus"]
  a <- a[a != "Paraventricular hypothalamic nucleus, magnocellular division"]
  a <- a[a != "Paraventricular hypothalamic nucleus, parvicellular division"]
  a <- a[a != "Paraventricular hypothalamic nucleus, descending division"]
  a <- c(a, df[grepl('Dorsomedial nucleus of the hypothalamus',df$region),'region'])
  a <- a[a != "Dorsomedial nucleus of the hypothalamus"]
  a <- c(a, df[grepl('Anterior hypothalamic nucleus',df$region),'region'])
  a <- a[a != "Anterior hypothalamic nucleus"]
  a <- c(a, df[grepl('Supramammilary nucleus',df$region),'region'])
  a <- a[a != "Supramammilary nucleus"]
  a <- c(a, df[grepl('Ventromedial hypothalamic nucleus',df$region),'region'])
  a <- a[a != "Ventromedial hypothalamic nucleus"]
  a <- c(a, df[grepl('Parabrachial nucleus',df$region),'region'])
  a <- a[a != "Parabrachial nucleus"]
  a <- a[a != "Parabrachial nucleus, lateral division"]
  a <- a[a != "Parabrachial nucleus, medial division"]
  
  df <- df[!df[,'region'] %in% a,]
  
  df <- FromDataFrameTable(df, pathName = 'path')

  if(output == 'df'){
    df <- ToDataFrameTable(df, region = 'region')
  }
  
  return(df)
}

color_annotation <- function(v, colcode = 'ggsci', seed = 1987){
  # INPUT: a factor vector
  # OUTPUT: a color code for each level of the factor (used to give a color to each anatomical group)
  # THIS FUCKING FUNCTION NEEDS SOME FUCKING FIXING !!!!!!!!!!!!!!!!!!!!!!!
  v <- factor(v)
  
  if(colcode == 'ggsci'){
    set.seed(seed)
    color <- sample(unique(c(pal_npg()(8), pal_aaas()(10), pal_gsea()(12), pal_tron()(7), pal_nejm()(8))), 
                    length(levels(v)))
  }
  else{
    if(all(is.character(colcode), (length(unique(v)) == length(levels(colcode))))){
      color <- colcode
    }
  }
  
  #HERE i WANT TO PUT SOME TOOL TO MEASURE DISTANCE BETWEEN 
  #THE COLORS AND BETTER SELECT THEM BY THEIR CONTRAST
  levels(v) <- color 
  return(as.character(v))
}

mfiles_to_dataframe <- function(IDpath, datapath, anato_path, region = 'region', balanced = T, 
                                drop_Layer_1 = T, anatomic_groups = T, del_unlabeled = T, 
                                colors = T, activity_zero = F, trim = c('trimmer', 'anatomical_Oh'), 
                                add_trim = NULL, by = 'acronym'){
  #
  # Input: Path to the file with regions IDs and to the data files.
  #
  # Output: data frame in which the data has been summed from children to parents across the whole tree and then trimmed until the levels of interest
  #
  # NOTE: the triiming function should be revised if your regions of interest are different or if you are NOT using the Allen Brain Atlas foradult mouse
  # NOTE2: Layer_1 == T will exclude Layer 1 of the whole isocortex. It can be a way to avoid ring effect in the imaging if they result in non-trivial 
  #amount of false positives
  

  d <- lapply(datapath, function(i){
    load_multiple_files(IDpath, i, balanced = balanced, del_unlabeled = del_unlabeled)  
  })
  
  b <- list(annotation = d[[1]]$annotation,
            counts = data.frame(rep(0, nrow(d[[1]]$counts))))
  for(i in d){
    b$counts <- cbind(b$counts, i$counts)
  }
  b$counts <- b$counts[,-1]
  #fazer essa funcao dar output de list separando o annotation
  
  if(isTRUE(anatomic_groups)){
    suppressWarnings(an <- anatomical_groups(IDpath, del_unlabeled = del_unlabeled, add_acronym = FALSE))
  }

  b <- data.frame(b$annotation, anatomical_groups = an, b$counts, stringsAsFactors = F)
  
  #This part excludes Cortical Layer1
  if(isTRUE(drop_Layer_1)){
    suppressWarnings(l <- drop_Layer_1(b, region = region))
    b <- b[!b$region %in% l,]
  }
  
  #summing the counts onto parent nodes
  for (i in 7:ncol(b)){
    a <- hierarchical_sum(data.frame(b[,1:3], b[,(i)], stringsAsFactors = F))
    b[,(i)] <- a[a$name == b$id,'s']
  }
  
  if(trim == 'trimmer'){
    # Trimming CHILDREN of the resulting Data Frame to the wanted level
    suppressWarnings(a <- trimmer(b, region = region, output = 'df'))
    # Deleting child nodes which are not the level counted by ClearMap
    z <- remove_noncounted_child(b)
    
    b <- b[b$acronym %in% z,]
    b <- b[b$region %in% a,]
  }
  else{
    if(tolower(trim == 'anatomical_Oh')){
      a = anatomical_oh_trim(b, anato_path, add_trim = add_trim, by = by, output = 'df')
      b <- b[b$acronym %in% a$acronym,]
    }
  }
  
  # Creating factors dataframe
  f <- lapply(datapath, function(i){
    factors_from_names(i)  
  })
  
  ffm <- data.frame(sapply(colnames(f[[1]]), function(i) { i = factor() }))
  for(i in f){
    ffm <- rbind(ffm, i)
  }
  
  # Creating list with 1) annotation, 2) Factors and 3) data
  j <- list()
  
  j$annotation <- b[,1:6]
  j$labels <- ffm
  j$counts <- b[,7:ncol(b)]
  j$counts <- data.frame(t(j$counts))
  colnames(j$counts) <- j$annotation$acronym
  rownames(j$counts) <- c(1:nrow(j$counts))
  
  if(isTRUE(colors)){
    j$annotation <- cbind(j$annotation, colors = as.character(color_annotation(j$annotation$anatomical_groups)))
  }
  j$labels$signal <- factor(j$labels$signal, levels(j$labels$signal)[c(2,1)])
  
  if(isFALSE(activity_zero)){
    j$counts = j$counts[colSums(j$counts) > 0]
    j$annotation <- inner_join(j$annotation, data.frame(acronym = colnames(j$counts)))
  }
  
  return(j)
}

anatomical_tree <- function(file, header = T, fill = F, check.names = F, set_name = 'region',
                            del_unlabeled = T){
  #
  # INPUT: .csv or .txt file
  #
  # OUTPUT: tree structure with a given var as levelNames if desired
  if(tools::file_ext(file) == "txt"){
    l <- read.table(file=file, header=header, fill=fill, check.names=FALSE)
  }
  else{
    if(tools::file_ext(file) == "csv"){
      l <- read.csv(file=file, sep = ",", header=header, fill=fill, check.names=FALSE)
    }
  }
  
  l <- data.frame(id = as.integer(l[,1]),parent_id = as.integer(l[,8]), region = as.character(l[,2]), 
                  acronym = as.character(l[,3]), parent_acronym = as.character(l[,9]), stringsAsFactors = F)
  if(del_unlabeled){ 
    l <- l[!l$region %in% 'Unlabeled',]
    }
  
  tree <- FromDataFrameNetwork(l[-1,])
  

  
  if(!is.null(set_name)){
    tree$Do(function(x){ x$name = x[[set_name]]}) 
  }
  return(tree)
}

rename_tree <- function(tree, newname = 'region'){
  #
  #
  return(tree$Do(function(x){ x$name = x[[newname]]}) )
}

collapse_hemisphere <- function(datalist){
  
  var <- data.frame(cbind(datalist$labels, datalist$counts))
  
  
  df <- data.frame(var %>%
                     group_by(experiment, subject, group, signal) %>%
                     summarise_if(is.numeric,mean))
  
  #alternative way
  #  df1 <- aggregate(.~ Exp + subject + group + signal, var, mean)
  
  datalist$labels <- data.frame(df[,1:4])
  datalist$counts <- data.frame(df[,5:ncol(df)])
  
  return(datalist)
}

z_score <- function(x){  
  z <- c((x - mean(x))/sd(x))
  return(z)
}
z_normalization <- function(datalist, factors = c('experiment', 'group', 'signal')){
  # INPUT: list of dataframes
  # 
  # OUTPUT: an additional dataframe in the list with the numeric data normalized by condition
  # OBS: gives out a warning that they ignored the groups, but they actually did not from my verification
  
  df <- datalist
  df$labels$interac <- factor(with(df$labels, interaction(df$labels[,match(factors, names(df$labels))])))
  df <- cbind(df$labels, df$counts)
  df <- df %>% dplyr::arrange(interac)
  datalist$labels <- df[,1:ncol(datalist$labels)]
  datalist$counts <- df[,(ncol(datalist$labels)+2):ncol(datalist$counts)]
  df <- aggregate(formula = .~ interac, data = df, FUN = z_score)
  df <- df[,7:ncol(df)]
  df <- data.frame(sapply(df, unlist))
  
  datalist$zcounts <- data.frame(df)
  
  return(datalist)
  
  ### THIS IS TO VERIFY THAT THE FUNCTION IS OUTPUTTING THE COORECT VALUES
  #  a <- sapply(seq_along((levels(datalist$labels$signal))), function(i){
  #    z_score(df[datalist$labels$signal == levels(datalist$labels$signal)[i] ,'VISp'])
  #  })
  #  a <- cbind(datalist$labels$signal,datalist$counts$VISp, df$VISp,c(t(a)))
  #  cbind(a[a[,1] == 1,], z_score(a[a[,1] == 1,2]))
}

##############################################################################################################################################
##############################################################################################################################################
##########       FUNCTIONS BELOW THIS LINE ARE NOT CODED PROPERLY. USE AT YOUR OWN RISK


anatomical_groups <- function(IDpath, del_unlabeled = T, add_acronym = FALSE){
  #
  # INPUT: path to regions ID file
  # OUTPUT: a tree where the children of the nodes selected as region groups have the selected node's name.
  #
  
  t <- anatomical_tree(IDpath, set_name = 'region', del_unlabeled = del_unlabeled)
  
  # THE CODE BELOW DOES NOT WORK ANYMORE.
  # NEED TO FIND ANOTHER TRIMMING METHOD
  
  #print(t$`Basic cell groups and regions`$Cerebrum$`Cerebral cortex`, pruneFun = function(x) x$level < 3)
  t$`Basic cell groups and regions`$Cerebrum$`Cerebral nuclei`$Do(function(x) x$name = 'Cerebral Nuclei')
  t$`Basic cell groups and regions`$Cerebrum$`Cerebral cortex`$`Cortical subplate`$Do(function(x){ x$name = 'Cortical subplate' })
  t$`Basic cell groups and regions`$Cerebrum$`Cerebral cortex`$`Cortical plate`$`Olfactory areas`$Do(function(x){ x$name = 'Olfactory regions' })
  t$`Basic cell groups and regions`$Cerebrum$`Cerebral cortex`$`Cortical plate`$`Hippocampal formation`$Do(function(x){ x$name = 'Hippocampal formation' })
  t$`Basic cell groups and regions`$Cerebrum$`Cerebral cortex`$`Cortical plate`$Isocortex$`Somatomotor areas`$Do(function(x){ x$name = 'Motor cortices' })
  t$`Basic cell groups and regions`$Cerebrum$`Cerebral cortex`$`Cortical plate`$Isocortex$`Somatosensory areas`$Do(function(x){ x$name = 'Sensory cortices' })
  t$`Basic cell groups and regions`$Cerebrum$`Cerebral cortex`$`Cortical plate`$Isocortex$`Frontal pole, cerebral cortex`$Do(function(x){ x$name = 'Polymodal association cortices' })
  t$`Basic cell groups and regions`$Cerebrum$`Cerebral cortex`$`Cortical plate`$Isocortex$`Infralimbic area`$Do(function(x){ x$name = 'Polymodal association cortices' })
  t$`Basic cell groups and regions`$Cerebrum$`Cerebral cortex`$`Cortical plate`$Isocortex$`Anterior cingulate area`$Do(function(x){ x$name = 'Polymodal association cortices' })
  t$`Basic cell groups and regions`$Cerebrum$`Cerebral cortex`$`Cortical plate`$Isocortex$`Prelimbic area`$Do(function(x){ x$name = 'Polymodal association cortices' })
  t$`Basic cell groups and regions`$Cerebrum$`Cerebral cortex`$`Cortical plate`$Isocortex$`Orbital area`$Do(function(x){ x$name = 'Polymodal association cortices' })
  t$`Basic cell groups and regions`$Cerebrum$`Cerebral cortex`$`Cortical plate`$Isocortex$`Agranular insular area`$Do(function(x){ x$name = 'Polymodal association cortices' })
  t$`Basic cell groups and regions`$Cerebrum$`Cerebral cortex`$`Cortical plate`$Isocortex$`Retrosplenial area`$Do(function(x){ x$name = 'Polymodal association cortices' })
  t$`Basic cell groups and regions`$Cerebrum$`Cerebral cortex`$`Cortical plate`$Isocortex$`Posterior parietal association areas`$Do(function(x){ x$name = 'Polymodal association cortices' })
  t$`Basic cell groups and regions`$Cerebrum$`Cerebral cortex`$`Cortical plate`$Isocortex$`Temporal association areas`$Do(function(x){ x$name = 'Polymodal association cortices' })
  t$`Basic cell groups and regions`$Cerebrum$`Cerebral cortex`$`Cortical plate`$Isocortex$`Ectorhinal area`$Do(function(x){ x$name = 'Polymodal association cortices' })
  t$`Basic cell groups and regions`$Cerebrum$`Cerebral cortex`$`Cortical plate`$Isocortex$`Perirhinal area`$Do(function(x){ x$name = 'Polymodal association cortices' })
  
  t$`Basic cell groups and regions`$Cerebrum$`Cerebral cortex`$`Cortical plate`$Isocortex$`Gustatory areas`$Do(function(x){ x$name = 'Sensory cortices' })
  t$`Basic cell groups and regions`$Cerebrum$`Cerebral cortex`$`Cortical plate`$Isocortex$`Auditory areas`$Do(function(x){ x$name = 'Sensory cortices' })
  t$`Basic cell groups and regions`$Cerebrum$`Cerebral cortex`$`Cortical plate`$Isocortex$`Visceral area`$Do(function(x){ x$name = 'Sensory cortices' })
  t$`Basic cell groups and regions`$Cerebrum$`Cerebral cortex`$`Cortical plate`$Isocortex$`Visual areas`$Do(function(x){ x$name = 'Sensory cortices' })
  
  t$`Basic cell groups and regions`$`Brain stem`$Interbrain$Thalamus$Do(function(x){ x$name = 'Thalamus' })
  #t$`Basic cell groups and regions`$`Brain stem`$Interbrain$Thalamus$`Thalamus, sensory-motor cortex related`$Do(function(x){ x$name = 'Thalamus, sensory-motor' })
  #t$`Basic cell groups and regions`$`Brain stem`$Interbrain$Thalamus$`Thalamus, polymodal association cortex related`$Do(function(x){ x$name = 'Thalamus, polymodal association' })
  t$`Basic cell groups and regions`$`Brain stem`$Interbrain$Hypothalamus$Do(function(x){ x$name = 'Hypothalamus' })
  t$`Basic cell groups and regions`$`Brain stem`$Midbrain$Do(function(x){ x$name = 'Midbrain' })
  t$`Basic cell groups and regions`$`Brain stem`$Hindbrain$Do(function(x){ x$name = 'Hindbrain' })
  
  if(add_acronym){
    b <- ToDataFrameTree(t, 'name', 'acronym')
    b <- b[,-1]
  }
  else{
    b <- ToDataFrameTree(t, 'name')
    b <- b[,-1]
  }

  return(b) 
}


