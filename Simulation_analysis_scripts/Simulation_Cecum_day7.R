
myargs <- commandArgs(trailingOnly = TRUE)
print(myargs)

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
library(ggplot2)
library(stringr)
library(parallel)

# ===== Loading input parameters =============================

# proportion of concentration of essential me= as.numeric(myargs[1])
econc_proportion = as.numeric(myargs[1])

# proportion of diet being added to gizzard (i.e., 0, if running sensitivity analyses on averaged samples)
diet_proportion = as.numeric(myargs[2])

# sample name, i.e. SRR22360819
sample = as.character(myargs[3])

# number of timesteps
timesteps = as.numeric(myargs[4])

# diet_name, i.e. corn_absorbed_2fibre; wheat_absorbed; corn_absorbed; etc
diet_name = as.character(myargs[5])


#### ===== Load input files ================================

diet <- read.table(paste0("ModelGeneration_Files/", diet_name,
                        "_diet_grower_aggregated_2023.txt"),
                        sep = '\t', header = T)
table_abund <- read.table("Day7_50_10_bins_info_merged_111indAssebmly_116pooled_5coassembly_5EcoliClusters_min0.005.txt",
                          sep='\t',header=T)

scfa_total_abundances <- read.table("ModelGeneration_Files/D7SCFA_TotalAbundanceRescaled.txt",
                                    sep = '\t', header = T)

# load sbml models to construct a community
model_list_cecum <- readRDS(paste0("ModelGeneration_Files/", sample,
                 "_day7_ind_assemblies_5coassembly_116pooled_5EcoliClusters_ModelList_gapfilled_", 
                 diet_name,"_gapseqAdapt_wEXreacs_May2024.rds"))

##########################################################
############ UPLOAD DATA ON CLUSTER ######################
############# RUN MODELS WITH REPLICATES #################
##########################################################

cores <- 5
replicates <- 5
stir_TF = TRUE


## define initial commsize based on total abundance per g/ceca rescaled from 1000 to 5000:
initial_commsize <- scfa_total_abundances$logAbund_rescaled[scfa_total_abundances$SampleName == sample]

# increase Rbflv concentration - metabolite with largest negative shadow price for multiple models
diet$met_mM[diet$met_id == "EX_cpd00220_e0"] <- 2*diet$met_mM[diet$met_id == "EX_cpd00220_e0"]


now <- Sys.time()
today <- strsplit(as.character(as.POSIXlt(now)[1]), "\\ ")[[1]][1]

cl <- makeCluster(cores, type="PSOCK",
                  outfile=paste0("cluster_logs/day7_", timesteps, "_gapseq_",
                                 sample,"_diet_", diet_name,"_propEconc", econc_proportion,
                                 "_diet_proportion", diet_proportion, "_",
                                 today,".log"))


# load all the scripts with functions
clusterCall(cl, function() { source("/scratch/j/jparkin/utkinair/Modeling/day7_Chicken_Cecum_BenWilling/ModelGeneration_scripts/Function_ConstructArena_universal_CecumMG.R") })
clusterCall(cl, function() { source("/scratch/j/jparkin/utkinair/Modeling/day7_Chicken_Cecum_BenWilling/ModelGeneration_scripts/Functions_RemoveBacMetabolites_TypeSpecieAssoc.R") })
clusterCall(cl, function() { source("/scratch/j/jparkin/utkinair/Modeling/day7_Chicken_Cecum_BenWilling/ModelGeneration_scripts/Function_RenameModelReactions.R") })

# export all the parameters and input files
clusterExport(cl, c("diet_proportion", "econc_proportion", "initial_commsize"))
clusterExport(cl, "model_list_cecum") # export models to clusters
clusterExport(cl, c("diet","diet_name", "sample")) # export diet and sample name to clusters
clusterExport(cl, c("table_abund","timesteps", "stir_TF")) # export abundance table, timesteps and stir variable


print(system.time(simlist <- parLapply(cl, 1:replicates, function(i){ #do N replicate simulations
  
  # set seed for this replicate, it will be used in the creation of initial arena (CreateModelArena function)
  N <- sample(c(1:100), 1)
  seedN <- as.integer(239 + N)
  

  set_pH <- function(arena, pH) {
    arena <- BacArena::addSubs(arena, smax=1000/(10^pH), mediac="EX_cpd00067_e0", unit="mM", add=F) 
    return(arena)
  }
  
  set_o2 <- function(arena, max_conc) {
    arena <- BacArena::addSubs(arena, smax = max_conc, mediac="EX_cpd00007_e0", unit="mM", add=F) 
    return(arena)
  }

  arena_cecum_list <- CreateModelArena_Cecum(model_list = model_list_cecum,
                                             sampleName = sample,
                                             diet = diet, proportion = diet_proportion,
                                             initial_commsize = initial_commsize,
                                             econc = econc_proportion, 
                                             stir_var = stir_TF,  arena_size=100, speed=5, tstep=1,
                                             seed = seedN)
  arena_cecum <- arena_cecum_list[[1]]
  essential_metabolites <- arena_cecum_list[[2]]
  
  # remove acetate and propionate from essential metabolites, if present
  essential_metabolites <- essential_metabolites[which(!(essential_metabolites %in% c("EX_cpd00029_e0", "EX_cpd00141_e0")))]
  
  # remove oxygen from essential metabolites, if present
  essential_metabolites <- essential_metabolites[which(!(essential_metabolites %in% c("EX_cpd00007_e0")))]
  
  
  ## add nucleobases, gam and CoA to essential metabolites - they're missing (shadow prices)

  essential_metabolites <- unique(c(essential_metabolites,
                                "EX_cpd00307_e0",  "EX_cpd00092_e0","EX_cpd00367_e0","EX_cpd00182_e0",
                                "EX_cpd00311_e0", "EX_cpd00010_e0","EX_cpd00276_e0"))
  ## add mucins and urea:
  essential_metabolites <- unique(c(essential_metabolites,
                                    "EX_cpd00122_e0", "EX_cpd00232_e0","EX_cpd00832_e0","EX_cpd02992_e0",
                                    "EX_cpd11842_e0", "EX_cpd21520_e0", "EX_cpd00073_e0"))

  out_cecum = vector("list", length = timesteps + 1)
  
  for (timestep in 1:timesteps) {
      print(paste0('timestep ', as.character(timestep)))

      # set cecum pH and oxygen levels:
      arena_cecum <- set_pH(arena_cecum, 7) # ph=7
      arena_cecum <- set_o2(arena_cecum,1e-08)
      
      ### simulate cecum for 1 h
      sim_cecum <- BacArena::simEnv(arena_cecum, time = 1) #simulate for 1h
      # save outputs of this step
      out_cecum[[timestep]]$Population = BacArena::plotCurves(sim_cecum, retdata = T, graph = F)$Population
      out_cecum[[timestep]]$Substances = BacArena::plotCurves(sim_cecum, retdata = T, graph = F)$Substances
      out_cecum[[timestep]]$Orgdat <- sim_cecum@simlist
      out_cecum[[timestep]]$Fluxlist <- sim_cecum@mfluxlist
      
      # save only low negative shadow prices       
      shadowlist <- reshape2::melt(unlist(lapply(sim_cecum@shadowlist[[2]], function(x) x[x < -0.1])))       
      shadowlist$Metabolite <-  unlist(lapply(rownames(shadowlist), function(x) {strsplit(x, "EX_")[[1]][2]}))       
      shadowlist$Specie <- unlist(lapply(rownames(shadowlist), function(x) { gsub("\\.$", "",  strsplit(x, "EX_")[[1]][1]) } ))       
      out_cecum[[timestep]]$Shadowlist <- shadowlist
      
      if (timestep == 1) {
        type_specie_table_cecum <- CreateTypeSpecieAssoc_table(sim_cecum)
      } 
      
      # select % of the grid that will be randomly removed once in 4 hours (released to colon)
      if (timestep %% 4 == 0) {
        perc_to_remove_cecum_bac = 0.25
        perc_to_remove_cecum_mets = 0.25
        # get specs, orgdat and metabolites removed from cecum in removefunction_result_cecum list
        removefunction_result_cecum <- RemoveBacMetabolites(sim_cecum, 
                                                            perc_to_remove_bac=perc_to_remove_cecum_bac,
                                                            perc_to_remove_mets=perc_to_remove_cecum_mets,
                                                            type_specie_table=type_specie_table_cecum)
        arena2_cecum <- removefunction_result_cecum[[1]]
      } else {
        arena2_cecum <- BacArena::getArena(sim_cecum, 1)
      }
      
      # add metabolites (diet without metabolites that are absorbed upstream)
      arena2_cecum <- BacArena::addSubs(arena2_cecum,smax = diet_proportion * diet$met_mM,
                                                            unit="mM",
                                                            mediac = as.character(diet$met_id), add = T, addAnyway = T)
      # add essential metabolites that aren't in the diet:
      essential_metabolites_not_in_diet <- setdiff(essential_metabolites,as.character(diet$met_id))
      arena2_cecum <- BacArena::addSubs(arena2_cecum,smax = econc_proportion,
                                        unit="mM",
                                        mediac = as.character(essential_metabolites_not_in_diet), add = T, addAnyway = T)

      # update cecum arena 
      arena_cecum <- arena2_cecum
      
}
      
    
  
  simlist = list()
  simlist[[1]] = out_cecum
  simlist[[2]] = seedN
  simlist[[3]] = essential_metabolites
  # save the list with all args:
  simlist[[4]] = list('econc_proportion' = econc_proportion,'diet_proportion' = diet_proportion,
                       'timesteps' = timesteps, 'sample' = sample, 'diet_name'=diet_name)
  return(simlist)
} ))) 

stopCluster(cl)

# save with the timestamp in the name
now <- Sys.time()
today <- strsplit(as.character(as.POSIXlt(now)[1]), "\\ ")[[1]][1]

saveRDS(simlist, paste0("sim_outputs/", 
                        sample, "_",as.character(diet_name), "_diet_grower_",
                        as.character(timesteps), "h_proportion", as.character(diet_proportion),
                        "_econc", econc_proportion,"_", today,
                        "_gapseq.rds"))





