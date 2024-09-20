library(data.table)
library(dplyr)
library(reshape2)
library(stringr)
library(parallel)
library(Rcpp)

myargs <- commandArgs(trailingOnly = TRUE)
print(myargs)

###################### retrieve input arguments ###############################
pattern_arg =  as.character(myargs[1])
diet = as.character(myargs[2])

######################################################################################################
######################################################################################################
path = "sim_outputs/"

# helper function to cbind data tables with different number of rows
custom_cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

getRelAbundance_Commsize <- function(files) {
  final_relabund = data.table()
  final_relabund_mean = data.table()
  commsize_total = data.table()
  final_absabund = data.table()
  
  # Helper function to process relative abundance data
  process_relabund_data <- function(relabund_data) {
    relabund_data$essential_conc <- essential_conc
    relabund_data$proportionMedium <- proportion_medium
    relabund_data$medium <- medium
    relabund_data$Sample <- sample_name
    return(relabund_data)
  }
  
  for (file in files) {
    print(file)
    
    
    simlist <- readRDS(paste0(path, file))
    replicates_number <- length(simlist)
    
    essential_conc = simlist[[1]][[4]][[1]]
    proportion_medium  =  simlist[[1]][[4]][[2]]
    timesteps = simlist[[1]][[4]][[3]]
    sample_name = simlist[[1]][[4]][[4]]
    medium = simlist[[1]][[4]][[5]]
    
    tax_abund_total <- data.table(); tax_relabund <- data.table();  commsize <- data.table()
    
    for (rep in 1:replicates_number) { # number of replicates
        # simlist[[replicate_number]][[simlist itself - 1, seed -2, essential_mets - 3]][[hour_number(last one)]]
        poplist <- as.data.frame(simlist[[rep]][[1]][[1]]$Population[,2])
        
        for (hour in 2:timesteps) {
          # if there are new species added in the second iteration - beware other added rownames! quick fix:
          poplist2 <- as.data.frame(simlist[[rep]][[1]][[hour]]$Population[,2])
          # customly cbind data frames with different numbers of rows:
          if (nrow(poplist2) > nrow(poplist)) {
            poplist_new <- custom_cbind.fill(poplist, poplist2)
            # then add missing rownames (species names), if any, from poplist2:
            rownames(poplist_new)[(nrow(poplist)+1):nrow(poplist2)] <- rownames(poplist2)[(nrow(poplist)+1):nrow(poplist2)]
            poplist <- poplist_new
          } else {
            poplist <- cbind(poplist, poplist2)
          }
        }
        
        
        colnames(poplist) <- c(1:timesteps)
        poplist[is.na(poplist)] <- 0
        
        tax_abund_comp_replicate <- cbind(poplist,
                                          rep(rep, nrow(poplist)))
        tax_abund_comp_replicate <- as.data.frame(tax_abund_comp_replicate)
        
        commsize_comp_replicate <- apply(tax_abund_comp_replicate[,c(1:timesteps)], 2, function(x) sum(as.numeric(x)))
        commsize_comp_replicate$Replicate = rep
        
        tax_abund_comp_replicate$Specie <- rownames(poplist)
        
        tax_abund_total <- rbind(tax_abund_total, tax_abund_comp_replicate)
        tax_abund_total[is.na(tax_abund_total)] = 0
        
        # turn absolute abundances at each timepoint into relative abundances
        tax_abund_comp_replicate[,c(1:timesteps)] <- 
          apply(tax_abund_comp_replicate[, c(1:timesteps)], 2, function(x) { as.numeric(x)/sum(as.numeric(x)) } )
        tax_relabund <- rbind(tax_relabund, tax_abund_comp_replicate)
        commsize <- rbind(commsize, commsize_comp_replicate)
        colnames(commsize) <- c(1:timesteps, "Replicate")
      
      
    }
    
    colnames(tax_relabund) <- c(1:timesteps, "Replicate", "Specie")
    tax_relabund_mean <- aggregate(. ~  Specie, data = tax_relabund, FUN=mean)
    tax_relabund_mean$Replicate <- NULL
    
    tax_relabund <- process_relabund_data(relabund_data = tax_relabund)
    tax_relabund_mean  <-  process_relabund_data(relabund_data = tax_relabund_mean)
    
    # get mean and sd of 10 replicates for each hour for community size table
    setDT(commsize)
    commsize_stat  <- commsize[, as.list(unlist(lapply(.SD, function(x) list(mean = mean(x), sd = sd(x))))),
                               .SDcols = c(1:timesteps)]
    commsize_stat_melt = reshape2::melt(commsize_stat)
    commsize_stat_melt[, stat:=gsub("^[0-9]*\\.", "", variable)]
    commsize_stat_melt[, variable:=gsub("\\.mean|\\.sd", "", variable)]
    commsize_stat_melt_wide <- tidyr::spread(commsize_stat_melt, key = stat, value = value)
    colnames(commsize_stat_melt_wide) <- c("Hour", "Commsize_mean", "Commsize_sd")
    commsize_stat_melt_wide <- process_relabund_data(relabund_data = commsize_stat_melt_wide)
    
    tax_abund_total <- process_relabund_data(relabund_data = tax_abund_total)

    final_absabund <- rbind(final_absabund , tax_abund_total)
    final_relabund <- rbind(final_relabund , tax_relabund)
    final_relabund_mean <- rbind(final_relabund_mean, tax_relabund_mean)
    commsize_total <- rbind(commsize_total, commsize_stat_melt_wide)
    
    setDT(final_absabund); setDT(final_relabund); setDT(final_relabund_mean); setDT(commsize_total)
  }
  
  
  return(list(
    final_relabund,final_relabund_mean, commsize_total, final_absabund
  ))
}




#####################################################################


getFluxes_Conc_Shadows <- function(files) {
  bac_fluxes_total_df <- data.table()
  bac_fluxes_total_df_sum <- data.table()
  bac_fluxes_total_sum_mean <- data.table()
  substance_table_df_total <- data.table()
  substance_table_df_total_mean <- data.table()
  bac_shadowprices_total_df <- data.table()
  
  
  for (file in files) {
    start_time <- Sys.time()
    print(file)
    
    simlist <- readRDS(paste0(path, file))
    replicates_number <- length(simlist)
    
    essential_conc = simlist[[1]][[4]][[1]]
    proportion_medium  =  simlist[[1]][[4]][[2]]
    timesteps = simlist[[1]][[4]][[3]]
    sample_name = simlist[[1]][[4]][[4]]
    medium = simlist[[1]][[4]][[5]]
    
    # Helper function to process relative abundance data
    process_data <- function(fluxconc_data) {
      setDT(fluxconc_data)
      fluxconc_data[, c('essential_conc', 'proportion_medium', 'medium', 'SampleName') :=
                      .(essential_conc, proportion_medium, medium, sample_name)]
      return(fluxconc_data)
    }
    

    bac_fluxes_total_file = data.table()
    substance_table_total_file = data.table()
    bac_shadowprices_total_file = data.table()
    for (rep in 1:replicates_number) { # number of replicates
      for (hour in 2:timesteps) {
        ## Get fluxes and substances
        fluxlist <- simlist[[rep]][[1]][[hour]]$Fluxlist[[2]]
        subconclist <- simlist[[rep]][[1]][[hour]]$Substances[, 2]
        shadowlist_rep <- simlist[[rep]][[1]][[hour]]$Shadowlist
        
        
        # Combine exchange and added reactions (fluxes)
        bac_fluxes_total <- lapply(names(fluxlist), function(bac) {
          if (length(names(fluxlist[[bac]])) != 0) {
            # Get the reaction names for the current bacteria
            bac_reactions <- names(fluxlist[[bac]])
            
            # line added Feb'24: for some reason when simulated with sec_obj='mtf', 
            # fluxlist triples (exactly) in size - with 2x more fluxes with no names (NA). Remove them:
            bac_reactions <- bac_reactions[!(is.na(bac_reactions))]
            
            # Extract relevant information
            bac_name <- gsub("(-gapseqAdRm|_cobrapy_adapted)", "", bac)
            is_exchange_reaction <- startsWith(bac_reactions, "EX_")
            
             
            data.table(
              Reaction = c(bac_reactions[is_exchange_reaction]),
              Flux = c(fluxlist[[bac]][is_exchange_reaction]),
              Specie = rep(bac_name, sum(is_exchange_reaction)),
              Replicate = rep(rep, sum(is_exchange_reaction)),
              Hour = rep(hour, sum(is_exchange_reaction))
            )
            
          }
        })
        
        # Combine and optimize data.table rows
        bac_fluxes_total <- rbindlist(bac_fluxes_total)
        setDT(bac_fluxes_total)
        # to save time and memory, remove the rows with reaction fluxes=0
        bac_fluxes_total[Flux!=0,]
        
        # Append to the results data.table
        bac_fluxes_total_file <- rbind(bac_fluxes_total_file, bac_fluxes_total)
        
        # Shadows processing
        # get lowest 5 values (shadows) for each specie
        bac_shadowprices <- shadowlist_rep %>% group_by(Specie) %>% top_n(-5,  value)
        
        bac_shadowprices <- data.table(
          Metabolite = bac_shadowprices$Metabolite,
          Shadow = bac_shadowprices$value,
          Specie = bac_shadowprices$Specie,
          Replicate = rep,
          Hour = hour
        )
        
        
        # Combine and optimize data.table rows
        # bac_shadowprices <- rbindlist(bac_shadowprices)
        
        # Append to the results data.table
        bac_shadowprices_total_file <- rbind(bac_shadowprices_total_file, bac_shadowprices)
        
        # Substances processing
        substance_table <- data.table(
          Metabolite = names(subconclist),
          Conc = subconclist,
          Replicate = rep(rep, length(subconclist)),
          Hour = rep(hour, length(subconclist))
        )
        
        # Append to the results data.table
        substance_table_total_file <- rbindlist(list(substance_table_total_file, substance_table))
      }
      
    }
    end_time <- Sys.time()
    print(end_time - start_time)
    
    
    
    # Summarize results
    # get sum of produced/consumed metabolites over N h:
    bac_fluxes_total_file[, Flux:=as.numeric(as.character(Flux))]
    
    bac_fluxes_total_file_sum <- aggregate(Flux ~ Reaction + Specie + Replicate,
                                           data = bac_fluxes_total_file, FUN = sum)
    data.table::setDT(bac_fluxes_total_file_sum)
    
    # add all the parameters of the simulation (essential_conc, proportion medium, SampleName, ...)
    bac_fluxes_total_file_sum <- process_data(bac_fluxes_total_file_sum)
    bac_fluxes_total_file <- process_data(bac_fluxes_total_file)
    
    # get mean and sd of total (sum) fluxes by replicate:
    bac_fluxes_total_file_sum_mean <- bac_fluxes_total_file_sum[, list(Flux.mean=mean(Flux),Flux.sd=sd(Flux)),
                                                                by=c("Reaction", "Specie","essential_conc",
                                                                     "proportion_medium","medium","SampleName")]
    data.table::setDT(bac_fluxes_total_file_sum_mean)
    
    # add all the parameters of the simulation (essential_conc, proportion medium, SampleName, ...)
    substance_table_total_file <- process_data(substance_table_total_file)
    # get mean and sd of Conc by replicate:
    substance_table_total_file_mean  <- substance_table_total_file[, list(Conc.mean=mean(Conc),Conc.sd=sd(Conc)),
                                                                   by=c("Metabolite","Hour","essential_conc",
                                                                        "proportion_medium","medium","SampleName")]
    
    bac_shadowprices_total_file <- process_data(bac_shadowprices_total_file)
    
    bac_fluxes_total_df <- rbind(bac_fluxes_total_df, bac_fluxes_total_file)
    bac_fluxes_total_df_sum <- rbind(bac_fluxes_total_df_sum, bac_fluxes_total_file_sum)
    bac_fluxes_total_sum_mean <- rbind(bac_fluxes_total_sum_mean, bac_fluxes_total_file_sum_mean)
    
    substance_table_df_total <- rbind(substance_table_df_total, substance_table_total_file)
    substance_table_df_total_mean <- rbind(substance_table_df_total_mean, substance_table_total_file_mean)
    
    bac_shadowprices_total_df <- rbind(bac_shadowprices_total_df, bac_shadowprices_total_file)
    
  }
  end1_time <- Sys.time()
  print(end1_time - start_time)
  
  
  return(list(bac_fluxes_total_df,
              bac_fluxes_total_df_sum, bac_fluxes_total_sum_mean,
              substance_table_df_total, substance_table_df_total_mean,
              bac_shadowprices_total_df))
}

files <- list.files(path, pattern = paste0(diet,"_",pattern_arg))

date_of_simulation <- strsplit(strsplit(files[1], "_")[[1]][length(strsplit(files[1], "_")[[1]]) - 1], "\\.rds")[[1]][1]
timesteps <- as.numeric(strsplit(strsplit(files[1], 'h_')[[1]][1], "_")[[1]][length(strsplit(strsplit(files[1], 'h_')[[1]][1], "_")[[1]])])

# create a directory within OUT_files/ directory with a name of simulation date to save files there:
sim_date <- gsub("-","",date_of_simulation)
dir.create(paste0(path,"OUT_files/", sim_date))


final_relabund_commsize <- getRelAbundance_Commsize(files = files)
saveRDS(final_relabund_commsize, paste0(path,"OUT_files/", sim_date, "/OUTfile_RelAbund_CommSize_",
                                        "allSamples_", diet,
                                        "_sim", date_of_simulation,"_",
                                        timesteps,"_", pattern_arg, "_",
                                        "allReplicates_andMean.rds"))


outlist <- getFluxes_Conc_Shadows(files = files)


saveRDS(outlist , paste0(path,"OUT_files/", sim_date, "/OUTfile_Fluxes_Substances_Shadow_",
                         "allSamples_",diet,
                         "_sim", date_of_simulation,"_",
                         timesteps,"_", pattern_arg, "_",
                         "allReplicates_andMean.rds"))


path_to_tables <- "sim_outputs/OUT_files/Tables/"

essential_mets_table <- c()
for (file in files) {
  start_time <- Sys.time()
  print(file)
  
  simlist <- readRDS(paste0(path, file))
  replicates_number <- length(simlist)
  essential_mets = simlist[[1]][[3]]
  sample_name = simlist[[1]][[4]][[4]]
  medium = simlist[[1]][[4]][[5]]
  
  essential_mets_table <- rbind(essential_mets_table, cbind(essential_mets, 
                                                            rep(sample_name, length(essential_mets)),
                                                            rep(medium, length(essential_mets))))
  
}

write.table(essential_mets_table, paste0(path_to_tables,"EssentialMetsTable_", diet,"_", pattern_arg, "_",  sim_date, ".txt"),
                                        sep='\t',quote=F)
# }

# ================== compile abundance tables for plots ====================

meltTheRelabTable <- function(relab_table, absabund=F) {
  relab_table <- as.data.frame(relab_table)
  colnames(relab_table)[timesteps +1] <- "Replicate"
  relab_table <- relab_table[, -which(colnames(relab_table) == "1")]
  
  
  if (absabund == F) {
    relab_table <- reshape2::melt(relab_table,
                                  id.vars = c( "Specie", "Sample", "proportionMedium", 
                                               "essential_conc", "medium"))
    colnames(relab_table)[c(6,7)] <- c("Hour", "RelAbundance") 
  } else {
    relab_table <- reshape2::melt(relab_table,
                                  id.vars = c("Specie", "Sample", "proportionMedium", 
                                              "essential_conc", "medium", "Replicate"))
    colnames(relab_table)[c(7,8)] <- c("Hour", "AbsAbundance") 
  }
  
  relab_table$Hour <- as.numeric(gsub('X','', as.character(relab_table$Hour)))
  return(relab_table)
}

setwd("sim_outputs/OUT_files/")
relab_files <- list.files(path = sim_date, pattern = paste0("OUTfile_RelAbund_CommSize_.*", diet, "_sim"))
relab_files <- relab_files[grep(pattern_arg, relab_files)]

absab_tables <- c()
relab_tables <- c()
commsize_tables <- c()
for (file in relab_files) {
  out_tables <- readRDS(file.path(sim_date, file))
  relab_table_mean <- out_tables[[2]]
  commsize_table <- out_tables[[3]]
  absab_table <- out_tables[[4]]
  
  
  relab_table_mean_melt <- meltTheRelabTable(relab_table_mean)
  colnames(absab_table)[max(relab_table_mean_melt$Hour)] <- "Replicate"
  
  absab_table <- meltTheRelabTable(absab_table, absabund=T)
  
  absab_tables <- rbind(absab_tables, absab_table)
  relab_tables  <- rbind(relab_tables, relab_table_mean_melt)
  commsize_tables <- rbind(commsize_tables, commsize_table)
}

write.table(relab_tables, paste0(path_to_tables,
                                 "RelAbund_tables_allSamples_", diet,"_", pattern_arg, "_",  sim_date, "_meanByRep.txt"),
            sep='\t',quote=F)
write.table(absab_tables, paste0(path_to_tables,
                                 "AbsAbund_tables_allSamples_", diet,"_", pattern_arg, "_",  sim_date, "_allReps.txt"),
            sep='\t',quote=F)
write.table(commsize_tables, paste0(path_to_tables,
                                    "CommSize_tables_allSamples_", diet,"_", pattern_arg, "_",  sim_date, "_meanByRep.txt"),
            sep='\t',quote=F)




########################################################################################################################
#########                    GET FLUXES AND CONCENTRATIONS FOR SCFA FLUXES AND CONC PLOTS                   ############   
########################################################################################################################
setwd("sim_outputs/OUT_files/")

out_tables_33samples_diet <- list.files(path = file.path(getwd(), sim_date),
                                          pattern = paste0("OUTfile_Fluxes_Substances_Shadow_allSamples_", diet, "_sim"))

out_tables_33samples_diet <- out_tables_33samples_diet[grep(pattern_arg,out_tables_33samples_diet)]

path_to_tables <- "sim_outputs/OUT_files/Tables/"


flux_table_total_scfa_file <- data.table()
flux_table_total_scfa_file_meanByRep_bySpecie <- data.table()
conc_table_scfa_file <- data.table()

for (out_table_file in out_tables_33samples_diet) {
  start_time = Sys.time()
  print(out_table_file)
  flux_tables <- readRDS(file.path(sim_date, out_table_file))
  # the object is a list with 6 slots -  for each compartment, each compartment has also 6 slots:
  # 1) all fluxes 2) total (sum) of all fluxes 3) mean by replicate of sum of all fluxes
  # 4) all concentrations 5) mean of all concentrations by replicate
  # 6) shadow prices

  flux_table <- as.data.table(flux_tables[[1]])
  
  flux_table_sum_eachHour_byCommunity <- flux_table[, .(Flux = sum(Flux)), by = .(Hour, Reaction, Replicate, proportion_medium, medium, SampleName)]
  
  flux_table_sum <- as.data.table(flux_tables[[2]])
  flux_table_sum_meanByRep <- as.data.table(flux_tables[[3]])
  conc_table_allmets <- as.data.table(flux_tables[[4]])
  conc_table_meanByRep <- as.data.table(flux_tables[[5]])
  conc_table_meanByRep_lastHour <- conc_table_meanByRep[which(conc_table_meanByRep$Hour == max(conc_table_meanByRep$Hour)),]
  
  scfa_reactions <- c("EX_cpd00211_e0", "EX_cpd00029_e0", "EX_cpd00141_e0", "EX_cpd00159_e0")
  
  flux_table_sum_scfa <- flux_table_sum[Reaction %in% scfa_reactions, ]
  flux_table_scfa <- flux_table[Reaction %in% scfa_reactions, ]
  flux_table_sum_meanByRep_scfa <- flux_table_sum_meanByRep[Reaction %in% scfa_reactions, ]
  conc_table_meanByRep_scfa <- conc_table_meanByRep[Metabolite %in% scfa_reactions, ]
 
  flux_table_total_sum_scfa <- flux_table_sum_scfa[, .(Flux = sum(Flux)), by = .(Reaction, Replicate, proportion_medium, medium, SampleName)]
  
  flux_table_total_scfa_file <- rbindlist(list(flux_table_total_scfa_file, flux_table_total_sum_scfa))
  flux_table_total_scfa_file_meanByRep_bySpecie <- rbindlist(list(flux_table_total_scfa_file_meanByRep_bySpecie, flux_table_sum_meanByRep_scfa))
  conc_table_scfa_file <- rbindlist(list(conc_table_scfa_file, conc_table_meanByRep_scfa))
   
  end_time = Sys.time()
  print(end_time-start_time)

 }

shadows <- flux_tables[[6]]

write.table(shadows, 
            paste0(path_to_tables, "ShadowsTable_allSamples_allReps_",diet,"_", pattern_arg, "_sim", sim_date, ".txt"),
            sep = '\t', quote=F)
write.table(flux_table_sum_eachHour_byCommunity,
  paste0(path_to_tables,"flux_table_allMets_byCommunity_allReps_", diet, "_", pattern_arg, "_sim", sim_date, ".txt"),
  sep = '\t', quote=F)
write.table(conc_table_meanByRep_lastHour,
            paste0(path_to_tables,"conc_table_allMets_lastHour_meanByRep_", diet, "_", pattern_arg, "_sim", sim_date, ".txt"),
            sep = '\t', quote=F) 
write.table(flux_table_sum, 
            paste0(path_to_tables, "flux_table_sum_allMets_BySpecie_allReps_",diet,"_", pattern_arg, "_sim",sim_date, ".txt"),
            sep = '\t', quote=F)
write.table(flux_table_total_scfa_file, 
            paste0(path_to_tables, "flux_table_sum_SCFA_ByCommunity_",diet,"_", pattern_arg, "_sim",sim_date, ".txt"),
            sep = '\t', quote=F)
write.table(flux_table_total_scfa_file_meanByRep_bySpecie, 
            paste0(path_to_tables, "flux_table_sum_SCFA_bySpecie_meanByRep_",diet,"_", pattern_arg, "_sim",sim_date, ".txt"),
            sep = '\t', quote=F)
write.table(conc_table_scfa_file, 
            paste0(path_to_tables, "conc_table_48h_SCFA_",diet,"_", pattern_arg, "_sim",sim_date, ".txt"),
            sep = '\t', quote=F)

