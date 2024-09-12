library(BacArena)
library(parallel)
library(glpkAPI)
library(Rcpp)
library(igraph)
library(reshape2)
library(libSBML)
library(sybilSBML)


myargs <- commandArgs(trailingOnly = TRUE)
print(myargs)

## ==========

diet_name = as.character(myargs[1])

##########################################################################################
###############   FOR GAPSEQ-ADAPTED (modified) MODELS  ##################################
###############         GAPFILLED USING 4 MEDIA         ##################################
######        INDIVID ASSEMBLY + COASSEMBLY (5 common species)       #####################
#############    + 116 pooled species from pooled reads           ########################
#############             + 5 Ecoli clusters from pooled reads       #####################
#############                      MAY 2024                       ########################
####################           RUN ON NIAGARA          ###################################
##########################################################################################

table_all_samples_abund <- read.table("Day7_50_10_bins_info_merged_111indAssebmly_116pooled_5coassembly_5EcoliClusters_min0.005.txt", sep='\t',header=T)

# save list of pre-loaded models for each sample and the table with taxonomic category for each specie
for (sample in sort(unique(table_all_samples_abund$sample))) {
  
  print(sample)
  model_list <- c()
  models <- table_all_samples_abund$binID_new[which(table_all_samples_abund$sample == sample)]
  
  for (model in models) {
    # if model is from individual assemblies:
    if (length(grep("SRR",model)) > 0) {
      print(paste0("SRR bin: ", model))
      bac_gapseq <- readSBMLmod(paste0("models_adapted_gapfilled_v2_4media_wMissingEXreacs_sbml_uniqueSpecies/",
                                       model,"_", diet_name, "_gapFilled-adapt_EXadded.xml"))
      # if model comes from pooled reads for common species from individual assemblies:
    } else if (length(grep("pooled",model)) > 0) {
      print(paste0("pooled specie: ", model))
      model_name <- gsub("_pooled","",model)
      bac_gapseq <- readSBMLmod(paste0("models_adapted_gapfilled_v2_4media_wMissingEXreacs_117pooledSpecies/",
                                       model_name,"_", diet_name, "_gapFilled-adapt_EXadded.xml"))
    } else if (length(grep("Ecoli_cluster",model)) > 0) {
      bac_gapseq <- readSBMLmod(paste0("models_adapted_gapfilled_v2_4media_5Ecoli_clusters_pooled/",
                                       model,"_", diet_name, "_gapFilled-adapt.xml"))
    } else {
      print(paste0("coassembly: ", model))
      bac_gapseq <- readSBMLmod(paste0("models_adapted_gapfilled_v2_4media_wMissingEXreacs_coassembly_5species/",
                                       model,"_", diet_name, "_gapFilled-adapt_EXadded.xml"))
    }
    
    
    bac_gapseq@mod_name <- unique(as.character(table_all_samples_abund$tax_specie[table_all_samples_abund$binID_new == model]))
    bac_gapseq@mod_desc <- unique(as.character(table_all_samples_abund$tax_specie[table_all_samples_abund$binID_new == model]))
    
    model_list <- c(model_list, bac_gapseq)
  }
  
  # save list of pre-loaded models for this compartment:
  saveRDS(model_list,
          paste0("ModelGeneration_Files/", sample,
                 "_day7_ind_assemblies_5coassembly_116pooled_5EcoliClusters_ModelList_gapfilled_", 
                 diet_name,"_gapseqAdapt_wEXreacs_May2024.rds"))
  
}
