
library(plyr)
library(dplyr)
library(magrittr)
library(BacArena)
library(glpkAPI)
library(Rcpp)
library(igraph)
library(reshape2)
library(libSBML)
library(sybilSBML)
library(stringr)
library(parallel)

RenameModelReactions = function(model) {
  model@react_id <- gsub('__40__e__41__', '(e)', model@react_id) # for AGORA models
  model@react_id <- gsub('__40__c__41__', '(c)', model@react_id) # for AGORA models
  model@react_id <- gsub('__91__c__93__', '(c)', model@react_id) # for AGORA models
  model@met_id <- gsub('__91__', '(', model@met_id) # for AGORA models
  model@met_id <- gsub('__93__', ')', model@met_id) # for AGORA models
  
  model@react_id <- gsub('__', '_',model@react_id) # for models downloaded from Machado et al
  model@met_id <- gsub('__', '_',model@met_id) # for models downloaded from Machado et al
  model@met_id <- gsub('\\[', '(',model@met_id) # for models downloaded from Machado et al
  model@met_id <- gsub('\\]', ')',model@met_id) # for models downloaded from Machado et al
  return(model)
}