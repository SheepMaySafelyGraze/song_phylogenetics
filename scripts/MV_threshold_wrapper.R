##########################
#
#
# run multivariate threshold model
#
#
##########################

library(ape)
library(stringr)

# read data and tree save to tmp files, choose motif to analyse
source('utils/data_handler.R')
data_path <- 'data/suboscine_data.csv'
tree_path <- 'data/suboscine_tree.tree'

analysis_scripts = data.frame(row.names=c('MCMC', 'ML'),
                              loc=c('rb-scripts/lines/MCMC.Rev',
                                    'rb-scripts/lines/ML.Rev'))

model_type = 'MV_threshold'
null_dist = FALSE
motif_inds = 0:1
n_runs = 1

for (analysis in rep(c('MCMC'), n_runs)){
  dat <- read.csv(data_path)
  tree<- read.tree(tree_path)
  
  res <- motif_data_to_binary(tree, dat)
  
  tree <- res[[1]]
  motif_data <- res[[2]]
  motif_data <- motif_data[,motif_inds+1]   # subset to motifs of interest
  
  write.tree(tree, 'data/tree_tmp.tree')
  write_mv_data_to_nexus(motif_data)
  
  
  # parameters
  n_chains = 1
  n_generations = 3e5
  burnin_generations = 300000 # for ML only
  tuning_interval = 1000 # for ML only
  printgen=500
  run_id = NULL ### populate with a run id, else date and time used
  
  if (is.null(run_id)) {
    run_id = date()
  }
  
  # folder of run results
  save_loc = paste('output/',
                   model_type,
                   '/',
                   analysis,
                   '_',
                   run_id,
                   sep=''
  )
  save_loc <- str_replace_all(save_loc, ':', '')
  
  # temporary Rev script
  tmp <- tempfile(fileext = '.Rev')
  
  con <- file(tmp, open='a')
  
  # write parameter values
  lines_0 <- c(
    paste('save_loc ="', save_loc, '"', sep=''),
    paste('n_chains =', n_chains, sep=''),
    paste('burnin_generations =', burnin_generations, sep=''),
    paste('n_generations =', n_generations, sep=''),
    paste('tuning_interval =', tuning_interval, sep=''),
    paste('printgen =', printgen, sep='')
  )
  
  writeLines(lines_0, con)
  
  # add model code
  lines_model <- readLines('rb-scripts/lines/MV_threshold.Rev')
  writeLines(lines_model, con)
  
  # add code to run analysis
  lines_analysis <- readLines('rb-scripts/lines/MCMC.Rev')
  writeLines(lines_analysis, con)
  
  
  close(con)
  
  
  # run script
  system(paste('rb', tmp))
  
  # wiriting motifs in single string
  # TODO: change this to support > 9 motifs
  motif_ind_str <- ""
  for (ind in sort(motif_inds)){
    motif_ind_str <- paste(motif_ind_str,as.character(ind),sep='')
  }
  
  # writing details of run to results list
  results <- read.csv('output/runs.csv', header=TRUE)
  results <- results[,2:ncol(results)]
  run_details <- c(date(),
                   model_type,
                   analysis,
                   run_id,
                   n_chains,
                   n_generations,
                   printgen,
                   tuning_interval,
                   save_loc,
                   data_path,
                   tree_path,
                   null_dist,
                   motif_ind_str)
  
  results[nrow(results)+1,] <- run_details
  write.csv(results, file='output/runs.csv')
  
  # remove script
  unlink(tmp)
}


