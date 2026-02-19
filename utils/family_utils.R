##########################
#
#
# utils to determine suboscine species
#
#
##########################

library(dplyr)
library(stringr)

# list of suboscine families
suboscine_families <- scan('data/suboscine_families.csv', sep=",", what="")


check_suboscine <- function(family, suboscine_families){
  # checks whether family suboscine
  # allows for family variable to contain extra words (e.g. common family name)
  
  for (spec in suboscine_families){
    if (grepl(str_to_lower(spec), str_to_lower(family))) {
      return(TRUE)
    }
  }
  return(FALSE)
}

which_suboscine <- function(families, suboscine_families){
  # determine indices in families for which family is suboscine
  
  inds <- which(sapply(families, \(x) check_suboscine(x, suboscine_families)))
  
  return(inds)
}


species_to_genus <- function(species_binomial) {
  noms <- str_split(species_binomial, '_')
  return(noms[[1]][1])
}

