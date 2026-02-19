##########################
#
#
# data handling utils
#
#
##########################

library(geiger)
library(stringr)
library(phangorn)
library(readr)
library(dplyr)

motif_data_to_binary <- function(tree, data){
  # data assumed to have 'species' and 'gmm_cluster' columns
  # returns a list containing:
  # - tree
  # - binary data for each motif as data frame, with rows ordered like tree$tip.label
  
  if (!(name.check(tree, data.names = data$species) == 'OK')){
    warning('Name mismatch in data and tree')
    return(NULL)
  }
  
  for (i in 0:7){
    motif_i_present = (data$gmm_cluster == i)
    data[paste('motif_', i,'_present', sep='')] <- motif_i_present
  }
  
  motif_data <- data %>%
    group_by(species) %>%
    summarise(motif_0 = any(motif_0_present),
              motif_1 = any(motif_1_present),
              motif_2 = any(motif_2_present),
              motif_3 = any(motif_3_present),
              motif_4 = any(motif_4_present),
              motif_5 = any(motif_5_present),
              motif_6 = any(motif_6_present),
              motif_7 = any(motif_7_present))
  
  
  # sort motif_data so that entries are in the same order as tree labels
  motif_data <- data.frame(motif_data)
  
  species <- tree$tip.label
  
  perm <- match(species, motif_data$species)
  
  motif_data <- motif_data[perm,]
  motif_data <- motif_data[,2:9, drop=TRUE]
  
  # convert to 1/0 and set rownames
  motif_data <- data.frame(lapply(motif_data, as.integer))
  
  rownames(motif_data) <- species
  
  return(c(tree, motif_data, recursive=FALSE))
  }


write_to_nexus <- function(motif_data, motif_ind,
                          filename='data/motif_data_tmp.nex'){
  # expects motif data to have species as rownames, followed by binary columns
  # for each of 8 motifs in order
  # saves binary motif data for specified motif index as a nexus file
  
  if (length(motif_ind) > 1){
    warning('Use write_mv_data_to_nexus for multiple motifs')
    return(NULL)
  }
  
  write.nexus.data(setNames(motif_data[,motif_ind+1], rownames(motif_data)),
                   file=filename, format='standard')
  
  # replacing characters with '0 1'
  tmp_file <- read_file(filename)
  tmp_file <- str_replace(tmp_file, '0123456789', '0 1')
  write_file(tmp_file, filename)
}

write_mv_data_to_nexus <- function(motif_data, filename='data/motif_data_tmp.nex'){
  # write a data frame with binary columns to nexus file
  
  write.nexus.data(as.matrix(motif_data),
                   file=filename, format='standard')
  
  # replacing characters with '0 1'
  tmp_file <- read_file(filename)
  tmp_file <- str_replace(tmp_file, '0123456789', '0 1')
  write_file(tmp_file, filename)
}


match_ape_rb_nodes <- function(tree){
  # returns a permutation mapping the nodes labels from RevBayes
  # to those used in ape
  n <- length(tree$tip.label)
  
  write.tree(tree, file='data/tree_tmp_desc.tree')
  system('rb utils/get_descendants.Rev')
  
  desc <- read_file("output/descendants_tmp.file")

  desc_split <- str_split(desc, '],')
  desc <- list()

  
  # creating list of descendants from revbayes
  for (i in 1:length(desc_split[[1]])){
    str_i <- desc_split[[1]][i]
    str_i <- str_replace_all(str_i, '\\[', '')
    str_i <- str_replace_all(str_i, '\\]', '')
    str_i <- str_replace_all(str_i, ' ', '')
    desc[[i]] <- str_split(str_i, ',')[[1]]
  }
  
  # for ape tips
  desc_ape <- list()
  for (i in (n+1):(n+tree$Nnode)){
    tip_inds <- Descendants(tree, i)
    desc_ape[[i-n]] <- tree$tip.label[Descendants(tree, i, type='tips')[[1]]]
  }
  
  
  # convert to vectors after ordering tip names
  desc_vec <- vapply(desc, \(x) paste(sort(x), collapse='-'), FUN.VALUE=character(1))
  desc_ape_vec <- vapply(desc_ape, \(x) paste(sort(x), collapse='-'), FUN.VALUE=character(1))
  
  if (length(setdiff(desc_vec, desc_ape_vec)) > 0){
    warning('unequal sets of node descendants')
    return(NULL)
  }
  
  return(match(desc_ape_vec, desc_vec))
}

thresh_ace_to_binary <- function(liab_df, thresh_df){
  # converts threshold model liabilities to binary traits
  # expects inputs to be liability and threshold log from a run of RevBayes
  
  out_df <- liab_df
  
  if (nrow(out_df) != nrow(thresh_df)){
    warning('differing number of samples of liabilities and thresholds')
    return(NULL)
  }
  
  for (i in 1:nrow(out_df)){
    th_i <- (thresh_df$threshold)[i]
    out_df[i, 4:ncol(out_df)] <- ifelse(out_df[i, 4:ncol(out_df)] > th_i, 1, 0)
  }
  
  # drop metadata columns
  out_df <- out_df[, c(5:ncol(out_df)), drop=TRUE]
  
  return(out_df)
}

### functions to location and load run results from runs.csv

read_run_files <- function(filters, runs_df, join=FALSE,
                           raw=FALSE, burnin=0.2){
  # filters should be a list specifying values for the columns of runs.csv
  # function returns run files associated with those runs
  # if join = TRUE, attempts to join these files row-wise after discarding burnin
  # if raw = FALSE, turns liabilites from threshold model into binary motif presence
  
  
  # get filepaths
  keep_inds <- rep(TRUE, nrow(runs_df))
  
  for (col in names(filters)){
    keep_inds <- keep_inds & runs_df[[col]] == filters[[col]]
  }
  
  folders <- (runs_df$save_loc)[keep_inds]
  filepaths <- c()
  
  # TODO: standardise naming so this is not necessary
  if (filters[['analysis']] == 'MCMC'){
    if (filters[['model_type']] == 'Mk'){
      filename = 'states.txt'
    }
    if (filters[['model_type']] %in% c('threshold', 'threshold_alt')){
      filename= 'interior_liabilities.log'
    }
    if (filters[['model_type']] == 'MV_threshold'){
      filename='MV_thresh.log'
    }
  } else {
    filename = 'ml_res.file'
  }
  
  for (folder in folders){
    filepaths <- c(filepaths, paste(folder, filename, sep ='/'))
  }
  
  out = list()
  for (i in 1:length(filepaths)){
    file <- filepaths[i]
  
    out[[i]] <- read.csv(file, sep='\t')
  
    
    # if threshold model, combine to yield binary data
    if (filters[['model_type']] %in% c('threshold', 'threshold_alt')){
      if (!raw) {
        
        # if alternate threshold model, drop tip liabilities
        if (filters['model_type'] == 'threshold_alt') {
          n <- (ncol(out[[i]]) - 3)/2
          out[[i]] <- out[[i]][, append(1:4, (n+5):ncol(out[[i]]))]
        }
        
        threshold_i <- read.csv(paste(folders[i], 'threshold.log', sep='/'),
                                sep='\t')
        out[[i]] <- thresh_ace_to_binary(out[[i]], threshold_i)
        
      } else {
        # drop metadata columns and return raw liabilities
        out[[i]] <- out[[i]][, 5:ncol(out[[i]])]
        
      }
    }
  }
  
  
  if (join){
    out_i <- out[[1]][floor(burnin*nrow(out[[1]])):nrow(out[[1]]),]
    out_df <- out_i
    for (i in 2:length(out)){
      out_i <- out[[i]][floor(burnin*nrow(out[[i]])):nrow(out[[i]]),]
      out_df <- rbind(out_df, out_i)
    }
    return(out_df)
  } else { 
    return(out)} 

}

