library(ggplot2)
library(scales)
library(tidyverse)
library(dplyr)
library(tidyr)
library(vegan)
library(ComplexHeatmap)
library(circlize)



SCFAconc_experimental <- readxl::read_xlsx("D7SCFA_BioSampleIDs.xlsx", sheet = 1)
color_samples_byHB <- SCFAconc_experimental[, c("SampleName", "Bacteroides")]
color_samples_byHB$color <- ifelse(color_samples_byHB$Bacteroides == "HB", "gold", "forestgreen")
sample_categories <- color_samples_byHB[,1:2]

table_all_samples_abund <- read.table("ModelGeneration_Files/Day7_50_10_bins_info_merged_111indAssebmly_116pooled_5coassembly_5EcoliClusters_min0.005.txt",
                                      sep = '\t',header=T)

seed_metabolites <- readr::read_delim("gapseq_github_dat/all_seed_metabolites_edited.tsv", delim='\t',col_names=T)
seed_metabolites <- as.data.frame(seed_metabolites)


####################################################################################################
################             NUMBER OF ECs ACROSS MODELS AND SAMPLES           #####################
####################################################################################################
# table generated in python via script Find_uniqueECs_in_XMLmodels.py
ECs_models <- read.delim("Table_ECnumbers_inEachModel_gapfilled_corn_absorbed_2fibre_May2024.txt", sep='\t',header=T)
ECs_models$Specie <- unlist(lapply(ECs_models$SpecieName, 
                                   function(x) unique(table_all_samples_abund$tax_specie[table_all_samples_abund$binID_new == x])))

ECs_samples <- list()
for (sample in sort(unique(table_all_samples_abund$sample))) {
  ECs_samples[[sample]] <- c()
  models_sample <- table_all_samples_abund$binID_new[table_all_samples_abund$sample == sample]
  for (model in models_sample) {
    ECs <- unlist(strsplit(ECs_models$EC_numbers[ECs_models$SpecieName == model], ", ")[[1]])
    ECs_samples[[sample]] <- unique(c(ECs_samples[[sample]],ECs))
  }
}

ECs_samples_models <- list()
for (sample in sort(unique(table_all_samples_abund$sample))) {
  ECs_samples_models[[sample]] <- list()
  models_sample <- table_all_samples_abund$binID_new[table_all_samples_abund$sample == sample]
  for (model in models_sample) {
    ECs <- unlist(strsplit(ECs_models$EC_numbers[ECs_models$SpecieName == model], ", ")[[1]])
    model_name <-  unique(table_all_samples_abund$tax_specie[table_all_samples_abund$binID_new == model])
    ECs_samples_models[[sample]][[model_name]] <- ECs
  }
}


########################################################################################################################
###############   GET THE DISSIMILARITY SCORE FOR EACH COMMUNITY (EACH SAMPLE) BASED ON ECs OVERLAP   ##################
########################################################################################################################
# Use the Max Adjustment if you want to compare communities directly by their internal dissimilarity,
# assuming maximum difference as the baseline.
# Use the Log Adjustment if you are more interested in scaling the impact of additional species
# in larger communities, reflecting diminishing returns in potential interaction complexity.
cmds_coef_calculation <- function(ec_data, 
                                  adjust_for_size= c("max","log"),
                                  method = "jaccard") {
  # Calculate Bray-Curtis dissimilarity
  # Assuming ec_data is presence/absence data frame where columns are species and rows are EC numbers
  dissimilarity <- vegan::vegdist(t(ec_data), method = method)
  dissimilarity_matrix <- as.matrix(dissimilarity)
  
  # Calculate raw CMDS
  cmds_raw <- mean(dissimilarity_matrix[lower.tri(dissimilarity_matrix)])
  
  # Number of species (N)
  N <- ncol(ec_data)
  
  if (adjust_for_size == "max") {
    # Adjusted CMDS by maximum possible dissimilarity
    cmds_adjusted_max <- cmds_raw / (N * (N - 1) / 2)
    cmds_return <- cmds_adjusted_max
  }
  if (adjust_for_size == "log") {
    # Adjusted CMDS by logarithmic scaling
    cmds_adjusted_log <- cmds_raw / log(N)
    cmds_return <-cmds_adjusted_log
  }
  return(cmds_return)
}

ec_data_long <- lapply(names(ECs_samples_models), function(sample) {
  # Loop through each species in the sample
  lapply(names(ECs_samples_models[[sample]]), function(species) {
    ec_numbers <- ECs_samples_models[[sample]][[species]]
    # Create a data frame where 'species' is repeated to match the length of 'ec_numbers'
    data.frame(Species = rep(species, length(ec_numbers)), EC = ec_numbers, stringsAsFactors = FALSE)
  })
}) %>%
  unlist(recursive = FALSE) %>%
  bind_rows() %>%
  group_by(Species, EC) %>%
  summarise(Count = n(), .groups = 'drop')  # Count each EC number per species, remove grouping

# Pivot to wide format
ec_data_wide_mat <- tidyr::pivot_wider(ec_data_long, names_from = Species, values_from = Count, values_fill = list(Count = 0))

# Replace NA with 0 if any EC numbers are missing in some species
ec_data_wide_mat[is.na(ec_data_wide_mat)] <- 0
ec_data_wide_mat <- as.data.frame(ec_data_wide_mat)
rownames(ec_data_wide_mat) <- ec_data_wide_mat$EC
# If we want binary presence/absence data instead of counts:
ec_data_wide_mat_binary <- ec_data_wide_mat[,-1]
ec_data_wide_mat_binary <- ec_data_wide_mat_binary %>%
  mutate(across(.cols = everything(), .fns = ~ ifelse(. > 0, 1, 0)))

#### now calculate cmds score for each sample:
# subset the ec_data_wide_mat_binary to species that are in the sample, and remove ECs (rows) that aren't present in any species
sample_cmds_scores <- c()
for (sample in names(ECs_samples_models)) {
  specie_names <- names(ECs_samples_models[[sample]]) 
  ec_data_sample_binary <- ec_data_wide_mat_binary[, specie_names]
  ec_data_sample_binary_nonzero <- ec_data_sample_binary[which(rowSums(ec_data_sample_binary) > 0),]
  
  # calculate cmds score:
  cmds_log = cmds_coef_calculation(ec_data = ec_data_sample_binary_nonzero, adjust_for_size = "log")
  
  sample_cmds_scores <- rbind(sample_cmds_scores,
                              cbind(sample, cmds_max, cmds_log))
}

sample_cmds_scores <- as.data.frame(sample_cmds_scores)
colnames(sample_cmds_scores)[1] <- "SampleName"
sample_cmds_scores <- merge(sample_cmds_scores, sample_categories)
sample_cmds_scores$cmds_log <- as.numeric(sample_cmds_scores$cmds_log)

wilcox.test(sample_cmds_scores$cmds_log[sample_cmds_scores$Category == "NB"],
            sample_cmds_scores$cmds_log[sample_cmds_scores$Category == "HB"])  # p-value = 0.2436

# Calculate summary statistics
summary_stats <- sample_cmds_scores %>%
  group_by(Category) %>%
  summarize(
    Mean = mean(cmds_log),
    SD = sd(cmds_log),
    ymin = Mean - SD,
    ymax = Mean + SD )


# Create the plot cmds vs category (NB/HB)
ggplot(summary_stats, aes(x=Category, y=Mean)) + 
  geom_col(fill = "darkred", position = position_dodge()) + 
  geom_errorbar(aes(ymin=ymin, ymax=ymax), width=0.0, position = position_dodge(0.9)) + 
  theme_bw() + 
  coord_flip() +
  labs(y = "Community metabolic dissimilarity score 
(log-adjusted for community size)") +
  theme(axis.text.x = element_text(size=20),
        axis.title.y=element_blank(),
        axis.title = element_text(size=22),
        axis.text.y = element_text(size=20))

write.table(sample_cmds_scores, "CommunityMetDiversityScores_log_adjusted__33samples_corn2.txt",
            sep ='\t',quote=F)

#################################################################################
######### Save the table with average number of unique ECs specie has per sample

# Function to find unique EC numbers in a sample
find_unique_ec_within <- function(sample_name, current_model, ec_numbers) {
  # Extract the EC numbers for the current sample
  current_ec_numbers <- ec_numbers[[sample_name]][[current_model]]
  # Combine EC numbers from all other samples
  other_ec_numbers <- unique(unlist(ec_numbers[[sample_name]][names(ec_numbers[[sample_name]]) != current_model]))
  # Find EC numbers unique to the current sample
  unique_ec <- setdiff(current_ec_numbers, other_ec_numbers)
  # Return the unique EC numbers
  return(unique_ec)
}

# Create a new list to store unique EC numbers for each sample
unique_ec_sample <- list()
for (sample_name in names(ECs_samples_models)) {
  unique_ec_sample[[sample_name]] <- lapply(names(ECs_samples_models[[sample_name]]), 
                                            find_unique_ec_within, sample_name = sample_name, ec_numbers = ECs_samples_models)
  names(unique_ec_sample[[sample_name]]) <- names(ECs_samples_models[[sample_name]])
}

# get the number of unique ECs per specie per sample:
ec_data <- lapply(names(unique_ec_sample), function(sample) {
  sapply(names(unique_ec_sample[[sample]]), function(species) {
    data.frame(
      Sample = sample,
      Species = species,
      UniqueECs = length(unique(unique_ec_sample[[sample]][[species]])),
      stringsAsFactors = FALSE
    )
  }, simplify = FALSE)
}) %>%
  unlist(recursive = FALSE) %>%
  bind_rows()


# Pivot the dataframe to wide format
ec_matrix <- ec_data %>%
  tidyr::pivot_wider(names_from = Sample, values_from = UniqueECs, values_fill = list(UniqueECs = NA))

ec_matrix$rowSum <- rowSums(ec_matrix[,-1], na.rm = T)
ec_matrix$rowMean <- rowMeans(ec_matrix[,-c(1,35)], na.rm = T)
ec_matrix$rowMedian <- unlist(apply(ec_matrix[,-c(1,35,36)], 1, function(x){median(x, na.rm = T)}))

write.table(ec_matrix, "Table_numberUniqueECs_eachSpecie_perSample_corn2_5Ecoli.txt",
            sep='\t',quote=F)






##############################################################################
#################         CALCULATE MRO SCORE     ############################
##############################################################################
# load all the exchange fluxes by specie
flux_allMets_bySpecie <- read.table("flux_table_sum_allMets_bySpecie_allReps_corn_absorbed_2fibre_sim13052024.txt",
                                          sep = '\t', header=T)
# function to add names of the metabolites
addMetNames_flux <- function(df) {
  df$MetName <- NA
  for (met in unique(df$Reaction)) {
    if (is.na(met) | met == "EX_cpd11416_c0") {next}
    df$MetName[df$Reaction == met] <- seed_metabolites$name[seed_metabolites$id == gsub('EX_|_e0','',met)]
  }
  return(df)
}

flux_allMets_bySpecie <- addMetNames_flux(flux_allMets_bySpecie)
# remove total fluxes < 1e-06:
flux_allMets_bySpecie <- flux_allMets_bySpecie[which(abs(flux_allMets_bySpecie$Flux) > 1e-06),]
# subset only consumption data
flux_allMets_bySpecie_consumed <- flux_allMets_bySpecie[flux_allMets_bySpecie$Flux < 0,]

# create list with metabolites names that are being consumed by species in each sample
comm_consum_list <- list()

for (sample in unique(flux_allMets_bySpecie_consumed$SampleName)) {
  comm_consum_list[[sample]] <- list()
  for (replicate in unique(flux_allMets_bySpecie_consumed$Replicate[flux_allMets_bySpecie_consumed$SampleName == sample])) {
    fluxes_samples_rep = flux_allMets_bySpecie_consumed[which(flux_allMets_bySpecie_consumed$Replicate == replicate &
                                                                flux_allMets_bySpecie_consumed$SampleName == sample),]
    comm_consum_list[[sample]][[replicate]] <- list()
    
    for (bac in unique(fluxes_samples_rep$Specie)) {
      comm_consum_list[[sample]][[replicate]][[bac]] <- 
        unique(fluxes_samples_rep$MetName[fluxes_samples_rep$Specie == bac])
    }
  }
}

## generate MRO coefficient for each sample-replicate:
# helper function
C_n_k = function(n, k) {
  factorial(n) / factorial(n-k) / factorial(k)
}

mro_samples <- list()
for (sample in names(comm_consum_list)) {
  mro_reps <- c()
  for (rep in 1:length(comm_consum_list[[sample]])) {
    
    comm_size = length(names(comm_consum_list[[sample]][[rep]]))
    
    sum_intersection_mets_consum <- 0
    sum_all_mets_consum <- 0
    bacs <- names(comm_consum_list[[sample]][[rep]])
    for (i in 1:(length(bacs)-1)) {
      sum_all_mets_consum <- sum_all_mets_consum + length(comm_consum_list[[sample]][[rep]][[bacs[i]]])
      for (j in (i + 1):length(bacs)) {
        
        # intersection of the metabolites consumed:
        intersection_mets_consum <- intersect(comm_consum_list[[sample]][[rep]][[bacs[i]]],
                                              comm_consum_list[[sample]][[rep]][[bacs[j]]])
        
        sum_intersection_mets_consum <- sum_intersection_mets_consum + length(intersection_mets_consum)
      }
    }
    sum_all_mets_consum = sum_all_mets_consum + length(comm_consum_list[[sample]][[rep]][[bacs[length(bacs)]]])
    
    mro <- (comm_size/C_n_k(comm_size,2)) * (sum_intersection_mets_consum/sum_all_mets_consum)
    mro_reps <- c(mro_reps, mro)
  }
  mro_samples[[sample]] <- mro_reps
}

MRO_df <- data.frame()
for (diet in names(comm_consum_list)) {
  MRO_df_diet <- data.frame("MRO" = as.numeric(unlist(mro_diets_samples)), 
                            "SampleName" = rep(names(comm_consum_list), each=5), 
                            "Bacteroides" = rep(color_samples_byHB$Bacteroides[order(color_samples_byHB$SampleName)], each =5))
  MRO_df <- rbind(MRO_df, MRO_df_diet)
}


########## Visualize
MRO_df$Bacteroides <- factor(MRO_df$Bacteroides,
                             levels = c("LB","HB"))

MRO_df$comm_size <- NA
for (sample in unique(MRO_df$SampleName)) {
  comm_size = length(names(comm_consum_list[[sample]][[1]]))
  MRO_df$comm_size[MRO_df$SampleName == sample] <- comm_size
}

wilcox.test(MRO ~ Bacteroides, MRO_df, exact = F) # p-value = 1.203e-05 

MRO_df <- MRO_df[with(MRO_df, order(Bacteroides, -MRO)),]
MRO_df$SampleName <- factor(MRO_df$SampleName, levels = unique(MRO_df$SampleName ))
MRO_boxplots <- ggplot(MRO_df,
                            aes(x = SampleName, y = MRO)) + 
  geom_boxplot(aes(color = Bacteroides),outlier.size = 0.5, lwd=0.7) +
  scale_color_manual(values = c("forestgreen", "gold")) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank(),panel.grid.major = element_blank()
  )


MRO_df_size <- MRO_df[, c("SampleName", "comm_size","Bacteroides")]
MRO_df_size <- MRO_df_size[!duplicated(MRO_df_size),]

comm_size_bars <- ggplot(MRO_df_size,
                         aes(x = SampleName, y = comm_size)) + 
  geom_bar(stat="identity", fill = "darkblue")+
  # scale_color_manual(values = c("forestgreen", "gold")) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank(),
        panel.grid = element_blank())

mro_plot <- cowplot::plot_grid(comm_size_bars, MRO_boxplots,nrow=2, align = "v",rel_heights = c(1,4))









