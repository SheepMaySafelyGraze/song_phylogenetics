### ------- resolving species synonmys and replacing with accepted names ------

library(phytools)
library(tidyverse)
library(taxize)

### load data and tree
trait_data <- read.csv('phylogenetic_analysis/traits_data_pc_gmm_8components_proba_100species.csv')
tree <- read.tree('phylogenetic_analysis/consensus_tree_100species.tre')

print(length(intersect(trait_data$species, tree$tip.label)))

unresolved_species_tree <- setdiff(tree$tip.label, trait_data$species)
unresolved_species_dat <- setdiff(trait_data$species, tree$tip.label)

### warning: VERY long runtime with manual input needed
# synonym_results <- synonyms(c(unresolved_species_dat, unresolved_species_tree), db='itis')

get_accepted_name_OLD <- function(species1, species2, syn_data){
  # either returns single, common accepted name
  # or raises an error if not synonymous
  
  # try to find each in database
  res1 <- syn_data[[species1]]
  res2 <- syn_data[[species2]]

  # get synonym names
  names1 = res1[, 'syn_name']
  names2 = res2[, 'syn_name']
  # if syn_name entry not present
  if (is.null(names1)|is.null(names2)){
    return(NA)
  }
  
  # if name in common, return (with underscore formatting)
  common_name <- intersect(names1, names2)[1]
  
  if (!is.na(common_name)){
    return(str_replace(common_name, ' ', '_'))
  } else {
    paste('No common name for:', species1, ' ', species2)
    return(NA)
  }
}

### resolving synonyms between tree and data and replacing with accepted names

replace_names <- function(tree, data, targets, replacement){
  data[data$species %in% targets, 'species'] <- replacement
  
  tree$tip.label[tree$tip.label %in% targets] <- replacement
  
  return(c(tree, data,recursive=FALSE))
}

get_accepted_name <- function(name1, name2) {
  # query gna 
  res <- gna_verifier(c(name1, name2))
  
  # if failed to resolve either
  if (any(is.na(res$currentCanonicalSimple))) {
    return(NA_character_)
  }
  
  # compare canonical names
  acc1 <- res$currentCanonicalSimple[res$submittedName == name1][1]
  acc2 <- res$currentCanonicalSimple[res$submittedName == name2][1]
  
  if (!is.na(acc1) && !is.na(acc2) && identical(acc1, acc2)) {
    return(str_replace(acc1, ' ', '_'))
  } else {
    return(NA_character_)
  }
}


# for each unmatched tip label, search unmatched data labels

for (tip.lab in unresolved_species_tree){
  print(paste('query: ', tip.lab))
  found <- 0
  for (dat.lab in unresolved_species_dat){
    if (found==0){
    tryCatch({acc_name <- get_accepted_name(tip.lab, dat.lab)},
             error = function(msg){acc_name=NA},
             warning = function(msg){acc_name=NA},
             finally = function(msg){acc_name=NA})

    if (!is.na(acc_name)){
      res <- replace_names(tree, trait_data, c(tip.lab, dat.lab), acc_name)
      tree <- res[[1]]
      trait_data <- res[[2]]
      found <- 1
      break
     }
    }
  }
  if(found==1){print(paste('accepted name:', acc_name))} else{
    print('No match found. Continuing.')
  }
  acc_name <- NA
}
  
# 9 species still not found

setdiff(tree$tip.label, trait_data$species)
setdiff(trait_data$species, tree$tip.label)

# manually matching
### WARNING: first pair SEEM TO BE TWO GENUINELY DIFFERENT SPECIES
res <- replace_names(tree, trait_data, c("Campylorhynchus_rufinucha", "Campylorhynchus_capistratus"),
              "Campylorhynchus_rufinucha")
tree <- res[[1]]
trait_data <- res[[2]]
res <- replace_names(tree, trait_data, c("Polioptila_plumbea", "Polioptila_bilineata"),
              "Polioptila_bilineata")
tree <- res[[1]]
trait_data <- res[[2]]
res <- replace_names(tree, trait_data, c("Trichastoma_celebense",  "Pellorneum_celebense"),
              "Pellorneum_celebense")
tree <- res[[1]]
trait_data <- res[[2]]
res <- replace_names(tree, trait_data, c("Rimator_malacoptilus", "Napothera_malacoptila"),
              "Napothera_malacoptila")
tree <- res[[1]]
trait_data <- res[[2]]
res <- replace_names(tree, trait_data, c("Stachyris_erythroptera", "Cyanoderma_erythropterum"),
              "Cyanoderma_erythropterum")
tree <- res[[1]]
trait_data <- res[[2]]
res <- replace_names(tree, trait_data, c("Zosterops_pallidus", "Zosterops_virens"),
              "Zosterops_virens")
tree <- res[[1]]
trait_data <- res[[2]]
res <- replace_names(tree, trait_data, c("Colluricincla_megarhyncha", "Colluricincla_rufogaster"),
              "Colluricincla_rufogaster")
tree <- res[[1]]
trait_data <- res[[2]]
res <- replace_names(tree, trait_data, c("Coracina_melanoptera", "Lalage_melanoptera"),
              "Lalage_melanoptera")
tree <- res[[1]]
trait_data <- res[[2]]
res <- replace_names(tree, trait_data, c("Sericornis_papuensis", "Aethomyias_papuensis"),
              "Aethomyias_papuensis")
tree <- res[[1]]
trait_data <- res[[2]]



setdiff(tree$tip.label, trait_data$species)
setdiff(trait_data$species, tree$tip.label)


# saving
write.csv(trait_data, 'phylogenetic_analysis/namefix_trait_data_100_1912.csv')
write.tree(tree, file='phylogenetic_analysis/namefix_tree_100_1912.tre')
  