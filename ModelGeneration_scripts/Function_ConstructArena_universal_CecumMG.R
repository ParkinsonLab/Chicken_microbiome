library(dplyr)
library(magrittr)
library(BacArena)
library(glpkAPI)
library(Rcpp)
library(reshape2)
library(libSBML)
library(sybilSBML)
library(stringr)


##############################################################################################################
####################       FUNCTION TO CONSTRUCT ARENA FOR CECUM. MIXING, DIET   #############################
##############################################################################################################

CreateModelArena_Cecum = function(model_list, sampleName, 
                            diet, proportion, initial_commsize=500, econc, 
                            bacteroides_perc=NULL, bacteroides_model=NULL,
                            stir_var = T, arena_size=100, speed=5, tstep=1, seed=239) {
  
  table_abund_sample <- table_abund[which(table_abund$sample == sampleName),]
  
  arena <- BacArena::Arena(n=arena_size, m=arena_size, stir = stir_var, Lx=0.025*(arena_size/100), Ly=0.025*(arena_size/100), tstep = tstep,
                           seed = seed)
  
  essential_meds <- c()
  for (i in 1:length(model_list)) {
        # find necessary models in the ModelList and upload them through Bac function
        # and then add to arena in the amount defined above (abundance)
    model_list[[i]] <- RenameModelReactions(model_list[[i]])
    abundance = table_abund_sample$relative_abund[table_abund_sample$tax_specie == model_list[[i]]@mod_name]

    if (ceiling(abundance*initial_commsize)>0){

      bac = BacArena::Bac(model=model_list[[i]], growtype="exponential")
      arena = BacArena::addOrg(arena, bac, amount=ceiling(abundance*initial_commsize))
      
      # save essential metabolites for each bac-organism
      # 'only_return=T' if essential metabolites should only be returned but not added to arena
      essential_meds <- unique(c(essential_meds, BacArena::addEssentialMed(arena, bac, limit = 100, only_return = T)))
      
   }
  }
  
  if (!(is.null(bacteroides_perc))) {
      abundance = bacteroides_perc/(1-bacteroides_perc)
      bacteroides_model <- RenameModelReactions(bacteroides_model)
      bac =  BacArena::Bac(model=bacteroides_model, growtype="exponential")
      arena = BacArena::addOrg(arena, bac, amount=ceiling(abundance*initial_commsize))
  }
  
  # remove acetate and propionate from essential metabolites if present!
  essential_meds <- essential_meds[which(!(essential_meds %in% c("EX_cpd00029_e0", "EX_cpd00141_e0")))]
  
  # save only essential metabolites that aren't in the diet:
  essential_meds <- essential_meds[which(!(essential_meds %in% diet$met_id))]

  # add essential metabolites 
  if (length(essential_meds) > 0) {
    arena <- BacArena::addSubs(arena,smax=econc,unit="mM", mediac=essential_meds,add=F, addAnyway = T)
  }
  
  # add nucleobases:
  # gapseq: add cytosine, uracil, cytidine, adenosine, guanosine and CoA
  arena <- BacArena::addSubs(arena,smax=econc,unit="mM", 
                             mediac=c("EX_cpd00307_e0",  "EX_cpd00092_e0","EX_cpd00367_e0",
                                      "EX_cpd00182_e0", "EX_cpd00311_e0","EX_cpd00010_e0"), add=T, addAnyway = T)
  
  # add metabolite concentrations from input diet
  arena <- BacArena::addSubs(arena,
                               smax=proportion * diet$met_mM,
                               unit="mM", mediac=diet$met_id, 
                               add = F, addAnyway = T)
  
  
  return(list(arena, essential_meds))
}
