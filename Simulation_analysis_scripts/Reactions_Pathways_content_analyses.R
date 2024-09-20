library(vegan)
library(igraph)
library(ggraph)
library(dplyr)
library(ggplot2)
library(tidyr)
library(cowplot)



################################################################################################
######   Get list of all reactions and metabolites in each model within each sample   ##########
################################################################################################
samples_list <- read.table("SRR_AccessionList_day7.txt")

# Initialize lists to store metabolites and reactions for each sample
samples_metabolites_list <- list()
samples_reacs_list <- list()
models_info_df <- data.frame()
samples_models_reactions_list <- list()

## loop through each sample and save the reaction_id and metabolite_id lists
for (sample in samples_list$V1) {
  sample_modelList_path <- paste0(path_to_models, sample, "_day7_ind_assemblies_5coassembly_116pooled_5EcoliClusters_ModelList_gapfilled_", 
                                  diet_name,"_gapseqAdapt_wEXreacs_May2024.rds")
  sample_modelList <- readRDS(sample_modelList_path)
  
  # Initialize lists to store metabolites and reactions for each model
  models_metabolites_list <- list()
  models_reacs_list <- list()
  samples_models_reactions_list[[sample]] <- list()
  
  for(model in sample_modelList){
    model_name <- model@mod_name
    models_metabolites_list[[model_name]] <- model@met_id
    models_reacs_list[[model_name]] <- model@react_id
    
    samples_models_reactions_list[[sample]][[model_name]] <- model@react_id
  }
  
  # ensure uniqueness across all models in a sample
  samples_metabolites_list[[sample]] <- unique(unlist(models_metabolites_list))
  samples_reacs_list[[sample]] <- unique(unlist(models_reacs_list))
  
  # remove "bio1" and "EX_cpd11416_c0" - biomass reacs
  samples_reacs_list[[sample]] <- samples_reacs_list[[sample]][-grep("bio1|EX_cpd11416_c0",samples_reacs_list[[sample]])]
  
}

# create a table for samples with their unique reaction and metabolites counts
samples_reacs_mets_info <- do.call(rbind, lapply(names(samples_reacs_list), function(sample_name) {
  data.frame(
    sampleName = sample_name,
    numberOfUniqueReactions = length(samples_reacs_list[[sample_name]]),
    numberOfUniqueMetabolites = length(samples_metabolites_list[[sample_name]])
  )
}))

################################################################################################
######          Get list of all pathways in each model within each sample             ##########
################################################################################################
# go through all models and save the pathways/subSystems in each model in the list
table_all_samples_abund <- read.table("Day7_50_10_bins_info_merged_111indAssebmly_116pooled_5coassembly_5EcoliClusters_min0.005.txt", sep='\t',header=T)

subSystems <- list()
for (model_name in sort(unique(table_all_samples_abund$binID_new))) {
  
  model <- readSBMLmod(paste0(path_to_models_folder,
                              model_name, "_gapFilled-adapt.xml"))
  # sum up the number of reactions within each subSystem
  subSys <- apply(model@subSys, 2, sum)
  
  # keep only subSystems with at least 3 reactions:
  subSys_min3 <- subSys[subSys>2]
  keep_subSystems <- names(subSys_min3)
  
  # save them to the list
  subSystems[[model_name]] <- keep_subSystems
}

# loop through all samples and summarize all the subSystems present in a sample:
subSystems_sample <- list()
for (sample in sort(unique(table_all_samples_abund$sample))) { 
  print(sample)
  for (model_name in table_all_samples_abund$binID_new[table_all_samples_abund$sample == sample]) {
    
    subSystems_sample[[sample]] <- union(subSystems_sample[[sample]],
                                         subSystems[[model_name]])
  }
}

########################################################################
########   PLOTS: reactions and pathways in HB vs NB samples  ##########
########################################################################
SCFAconc_experimental <- readxl::read_xlsx("/Users/irina/Desktop/Study/UofT/ParkinsonLab/Modeling_chickens/2023_ChickenCecalData_BenWilling/D7SCFA_BioSampleIDs.xlsx", sheet = 1)

color_samples_byHB <- SCFAconc_experimental[, c("SampleName", "Bacteroides","Func_capacity")]
color_samples_byHB$color <- "forestgreen"
color_samples_byHB$color[color_samples_byHB$Bacteroides == "HB"] <- "gold"
colors_samples <- setNames(as.vector(color_samples_byHB$color), color_samples_byHB$SampleName)


##################### density plots
## plot reacs numbers distributions
samples_reacs_mets_info$Bacteroides <- SCFAconc_experimental$Bacteroides[order(SCFAconc_experimental$SampleName)]

plot_reacs <- overlapping::overlap(list(samples_reacs_mets_info$numberOfUniqueReactions[samples_reacs_mets_info$Bacteroides == "NB"], 
                                        samples_reacs_mets_info$numberOfUniqueReactions[samples_reacs_mets_info$Bacteroides == "HB"]),plot=T)
reacs_density_plot <- ggplot(samples_reacs_mets_info) +
  geom_density(aes(x = numberOfUniqueReactions, fill = Bacteroides), alpha=0.6) +
  scale_fill_manual(values = c("NB" = "forestgreen","HB" = "gold")) + 
  theme_bw() +
  xlim(2000,3600) +
  labs(x = "Number of unique reactions") + 
  coord_flip() +
  scale_y_continuous(breaks = c(0,0.0015)) + 
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        legend.text = element_text(size=15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 19),
        axis.text = element_text(size = 17))

wilcox_reacs_NB_HB <- wilcox.test(numberOfUniqueReactions ~ Bacteroides, samples_reacs_mets_info, exact=F)
# p-value = 0.5035 

## plot pathways numbers distributions
samples_subsystems_info$Bacteroides <- SCFAconc_experimental$Bacteroides[order(SCFAconc_experimental$SampleName)]

pathways_density_plot <- ggplot(samples_subsystems_info) +
  geom_density(aes(x = numberOfUniquePathways, fill = Bacteroides), alpha=0.6) +
  scale_fill_manual(values = c("NB" = "forestgreen","HB" = "gold")) + 
  theme_bw() +
  xlim(0,700) +
  labs(x = "Number of unique pathways") + 
  coord_flip() +
  scale_y_continuous(breaks = c(0,0.006)) + 
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        legend.text = element_text(size=15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 19),
        axis.text = element_text(size = 17)
  )

wilcox_pathways_NB_HB <- wilcox.test(numberOfUniquePathways ~ Bacteroides, samples_subsystems_info, exact=F) 
# p-value = 0.004531 








############# ############# ############# 
#############   PCOA plots  ############# 
############# ############# ############# 
###### Presence/absence of reactions NB vs HB
## Create a binary matrix indicating presence/absence of reactions
# Unique list of all reactions across all samples
all_reactions <- unique(unlist(samples_reacs_list))

# Initialize a matrix of zeros
reaction_matrix <- matrix(0, nrow = length(all_reactions), ncol = length(samples_reacs_list), 
                          dimnames = list(all_reactions, names(samples_reacs_list)))

# Fill the matrix
for(i in names(samples_reacs_list)) {
  reactions_present <- samples_reacs_list[[i]]
  reaction_matrix[reactions_present, i] <- 1
}

reaction_matrix_presence <- reaction_matrix
# Calculate dissimilarity matrix using bray distance
diss_matrix <- vegdist(t(reaction_matrix), method = "bray")

# Perform PCoA
pcoa_results <- cmdscale(diss_matrix, eig = TRUE, k = 2) # k is the number of dimensions

eigenvalues <- pcoa_results$eig
# Calculate the percentages of variance explained by each axis
percent_variance_explained <- eigenvalues / sum(eigenvalues) * 100

# Convert to a dataframe for easier handling with ggplot2
pcoa_df <- as.data.frame(pcoa_results$points)
names(pcoa_df) <- c("PC1", "PC2")
pcoa_df$SampleName <- rownames(pcoa_df)

# Merge the PCoA scores with the group and color information
pcoa_data <- merge(pcoa_df, color_samples_byHB, by = "SampleName")

pcoa_plot <- ggplot(pcoa_data, aes(x = PC1, y = PC2, # color = Bacteroides, 
                                   fill = Bacteroides)) +
  geom_point(aes(shape = Bacteroides), size = 3, alpha = 0.8) +
  scale_shape_manual(values = c(21,24)) +
  scale_color_manual(values = c("NB" = "forestgreen","HB" = "gold")) +
  scale_fill_manual(values = c("NB" = "forestgreen","HB" = "gold")) +
  # scale_color_manual(values = c("blue","pink")) +
  ylim(-0.11,0.2) + 
  xlim(-0.3,0.34) + 
  theme_bw() +
  labs(title = "PCoA of samples by reaction content",
       x = paste("PCoA 1 (", format(percent_variance_explained[1], digits = 2), "%)", sep = ""),
       y = paste("PCoA 2 (", format(percent_variance_explained[2], digits = 2), "%)", sep = "")) +
  theme(legend.title = element_blank(),
        title = element_text(size=18),
        panel.grid.major = element_blank(),
        axis.text = element_text(size=18),
        legend.text = element_text(size=16),
        axis.title =  element_text(size=20)) +
  stat_ellipse(aes(group = Bacteroides, color = Bacteroides), type = "t", linetype = "dashed", level = 0.95)



######### Reactions in samples HB vs NB groups - weighted by relative abundances of species  
get_weightedByRelAbund_reacs_list <- function(samples_models_reactions_list) {
  # Initialize an empty list to store the weighted reaction sums for each sample
  weighted_reactions <- list()
  # Iterate through each sample in the samples_models_reactions_list
  for(sample_name in names(samples_models_reactions_list)) {
    # Get the species and their reactions in the sample
    species_reactions <- samples_models_reactions_list[[sample_name]]
    
    # Initialize an empty vector to store the weighted reactions
    weighted_reactions_sample <- c()
    
    # Iterate through each species in the sample
    for(model in names(species_reactions)) {
      # Get the relative abundance of the species in the sample
      rel_abundance <- table_all_samples_abund %>%
        filter(tax_specie == model, sample == sample_name) %>%
        pull(relative_abund)
      
      # Get the reactions for the species
      reactions <- species_reactions[[model]]
      
      # Weight the reactions by the relative abundance
      for(reaction in reactions) {
        if(reaction %in% names(weighted_reactions_sample)) {
          weighted_reactions_sample[reaction] <- weighted_reactions_sample[reaction] + rel_abundance
        } else {
          weighted_reactions_sample[reaction] <- rel_abundance
        }
      }
    }
    
    # Store the weighted reactions in the list
    weighted_reactions[[sample_name]] <- weighted_reactions_sample
  }
  return(weighted_reactions)
}

weighted_reactions <- get_weightedByRelAbund_reacs_list(samples_models_reactions_list)
all_reactions <- unique(unlist(lapply(weighted_reactions, names)))
reaction_matrix <- matrix(0, nrow = length(weighted_reactions), ncol = length(all_reactions), 
                          dimnames = list(names(weighted_reactions), all_reactions))

# Fill in the matrix with the weighted reaction counts
for(sample in names(weighted_reactions)) {
  for(reaction in names(weighted_reactions[[sample]])) {
    reaction_matrix[sample, reaction] <- weighted_reactions[[sample]][reaction]
  }
}

# Calculate dissimilarity using Bray-Curtis index
dissimilarity_matrix <- vegdist(reaction_matrix, method = "bray")

# Perform PCoA
pcoa_results <- cmdscale(dissimilarity_matrix, eig = TRUE, k = 2)  # for 2 dimensions
pcoa_df <- as.data.frame(pcoa_results$points)
names(pcoa_df) <- c("PC1", "PC2")
pcoa_df$SampleName <- rownames(pcoa_df)

eigenvalues <- pcoa_results$eig
# Calculate the percentages of variance explained by each axis
percent_variance_explained <- eigenvalues / sum(eigenvalues) * 100

# Merge the PCoA scores with the group and color information
pcoa_data <- merge(pcoa_df, color_samples_byHB, by = "SampleName")

pcoa_plot_weightedRelAb <- ggplot(pcoa_data, aes(x = PC1, y = PC2, 
                                                 #color = Bacteroides, 
                                                 fill = Bacteroides)) +
  geom_point(aes(shape = Bacteroides), size = 3, alpha = 0.8) +
  scale_shape_manual(values = c(21,24)) +
  scale_color_manual(values = c("NB" = "forestgreen","HB" = "gold")) +
  scale_fill_manual(values = c("NB" = "forestgreen","HB" = "gold")) +
  # scale_color_manual(values = c("blue","pink")) +
  theme_bw() +
  labs(title = "PCoA of samples by reaction content",
       x = paste("PCoA 1 (", format(percent_variance_explained[1], digits = 2), "%)", sep = ""),
       y = paste("PCoA 2 (", format(percent_variance_explained[2], digits = 2), "%)", sep = "")) +
  theme(legend.title = element_blank(),
        title = element_text(size=18),
        panel.grid.major = element_blank(),
        axis.text = element_text(size=14),
        legend.text = element_text(size=16),
        axis.title =  element_text(size=17)) +
  stat_ellipse(aes(group = Bacteroides, color = Bacteroides),type = "t", linetype = "dashed", level = 0.95)



#######  Repeat for pathways (subsystems):
get_weightedByRelAbund_SubSys_list <- function(SubSystems_samples_models) {
  # Initialize an empty list to store the weighted reaction sums for each sample
  weighted_subSys <- list()
  # Iterate through each sample in the samples_models_subSys_list
  for(sample_name in names(SubSystems_samples_models)) {
    print(sample_name)
    # Get the species and their subSys in the sample
    species_subSys <- SubSystems_samples_models[[sample_name]]
    
    # Initialize an empty vector to store the weighted subSys
    weighted_subSys_sample <- c()
    
    # Iterate through each species in the sample
    for(model in names(species_subSys)) {
      # Get the relative abundance of the species in the sample
      rel_abundance <- table_all_samples_abund %>%
        filter(gsub("_pooled","",binID_new) == model, sample == sample_name) %>%
        pull(relative_abund)
      
      # Get the subSys for the species
      subSys <- species_subSys[[model]]
      
      # Weight the subSys by the relative abundance
      for(pathway in subSys) {
        if(pathway %in% names(weighted_subSys_sample)) {
          weighted_subSys_sample[pathway] <- weighted_subSys_sample[pathway] + rel_abundance
        } else {
          weighted_subSys_sample[pathway] <- rel_abundance
        }
      }
    }
    
    # Store the weighted subSys in the list
    weighted_subSys[[sample_name]] <- weighted_subSys_sample
  }
  return(weighted_subSys)
}

weighted_subSys <- get_weightedByRelAbund_SubSys_list(SubSystems_samples_models)
all_subSys <- unique(unlist(lapply(weighted_subSys, names)))
subSys_matrix <- matrix(0, nrow = length(weighted_subSys), ncol = length(all_subSys), 
                        dimnames = list(names(weighted_subSys), all_subSys))

# Fill in the matrix with the weighted reaction counts
for(sample in names(weighted_subSys)) {
  for(pathway in names(weighted_subSys[[sample]])) {
    subSys_matrix[sample, pathway] <- weighted_subSys[[sample]][pathway]
  }
}

# Calculate dissimilarity using bray-curtis
dissimilarity_matrix <- vegdist(subSys_matrix, method = "bray")

# Perform PCoA
pcoa_results <- cmdscale(dissimilarity_matrix, eig = TRUE, k = 2)  # for 2 dimensions
pcoa_df <- as.data.frame(pcoa_results$points)
names(pcoa_df) <- c("PC1", "PC2")
pcoa_df$SampleName <- rownames(pcoa_df)

eigenvalues <- pcoa_results$eig

# Calculate the percentages of variance explained by each axis
percent_variance_explained <- eigenvalues / sum(eigenvalues) * 100

# Merge the PCoA scores with the group and color information
pcoa_data <- merge(pcoa_df, color_samples_byHB, by = "SampleName")

pcoa_plot_subSys_weightedRelAb <- ggplot(pcoa_data, aes(x = PC1, y = PC2, # color = Bacteroides, 
                                          fill = Bacteroides)) +
  geom_point(aes(shape = Bacteroides), size = 3, alpha = 0.8) +
  scale_shape_manual(values = c(21,24)) +
  scale_color_manual(values = c("NB" = "forestgreen","HB" = "gold")) +
  scale_fill_manual(values = c("NB" = "forestgreen","HB" = "gold")) +
  theme_bw() +
  labs(title = "PCoA of samples by pathway content",
       x = paste("PCoA 1 (", format(percent_variance_explained[1], digits = 2), "%)", sep = ""),
       y = paste("PCoA 2 (", format(percent_variance_explained[2], digits = 2), "%)", sep = "")) +
  theme(legend.title = element_blank(),
        title = element_text(size=18),
        panel.grid.major = element_blank(),
        axis.text = element_text(size=14),
        legend.text = element_text(size=16),
        axis.title =  element_text(size=17)) +
  stat_ellipse(aes(group = Bacteroides, color = Bacteroides),type = "t", linetype = "dashed", level = 0.95)




# create a matrix with only presence/absence of pathways (1 or 0, respectively)
subSys_matrix_presence <- subSys_matrix
subSys_matrix_presence[subSys_matrix_presence > 0] <- 1

# Calculate dissimilarity using bray-curtis
dissimilarity_matrix_presence <- vegdist(subSys_matrix_presence, method = "bray")

# Perform PCoA
pcoa_results <- cmdscale(dissimilarity_matrix_presence, eig = TRUE, k = 2)  # for 2 dimensions
pcoa_df <- as.data.frame(pcoa_results$points)
names(pcoa_df) <- c("PC1", "PC2")
pcoa_df$SampleName <- rownames(pcoa_df)

eigenvalues <- pcoa_results$eig

# Calculate the percentages of variance explained by each axis
percent_variance_explained <- eigenvalues / sum(eigenvalues) * 100

# Merge the PCoA scores with the group and color information
pcoa_data <- merge(pcoa_df, color_samples_byHB, by = "SampleName")

pcoa_plot_subSys <- ggplot(pcoa_data, aes(x = PC1, y = PC2, # color = Bacteroides, 
                                                        fill = Bacteroides)) +
  geom_point(aes(shape = Bacteroides), size = 3, alpha = 0.8) +
  scale_shape_manual(values = c(21,24)) +
  scale_color_manual(values = c("NB" = "forestgreen","HB" = "gold")) +
  scale_fill_manual(values = c("NB" = "forestgreen","HB" = "gold")) +
  theme_bw() +
  labs(title = "PCoA of samples by pathway content",
       x = paste("PCoA 1 (", format(percent_variance_explained[1], digits = 2), "%)", sep = ""),
       y = paste("PCoA 2 (", format(percent_variance_explained[2], digits = 2), "%)", sep = "")) +
  theme(legend.title = element_blank(),
        title = element_text(size=18),
        panel.grid.major = element_blank(),
        axis.text = element_text(size=14),
        legend.text = element_text(size=16),
        axis.title =  element_text(size=17)) +
  stat_ellipse(aes(group = Bacteroides, color = Bacteroides),type = "t", linetype = "dashed", level = 0.95)





#############################################################################################
############################### Run PERMANOVA using the adonis function and make plots
#############################################################################################
# Function to compute PERMANOVA and betadisper tests
perform_permanova_betadisper <- function(data_matrix, group_col) {
  # Compute dissimilarity matrix
  dist_matrix <- vegdist(as.matrix(data_matrix[, -which(colnames(data_matrix) == group_col)]))
  
  # Run PERMANOVA
  permanova_result <- adonis2(dist_matrix ~ Group, data = data_matrix, permutations = 999, method = "bray")
  print(permanova_result)
  
  # Run betadisper test
  dispersion <- betadisper(dist_matrix, group = data_matrix[[group_col]])
  TukeyHSD(dispersion)
  permutest(dispersion)
  
  return(dispersion)
}

# Function to generate dispersion plot
generate_dispersion_plot <- function(dispersion_obj, group_col, title) {
  dist_df <- as.data.frame(dispersion_obj$distances)
  dist_df$SampleName <- rownames(dist_df)
  dist_df$Bacteroides <- dispersion_obj$group
  
  plot <- ggplot(dist_df, aes(x = Bacteroides, y = dispersion_obj$distances)) + 
    geom_boxplot(aes(fill = factor(Bacteroides))) +
    labs(y = "Distance to centroid", x = "") +
    scale_fill_manual(values = c("NB" = "forestgreen", "HB" = "gold")) + 
    theme_classic() +
    coord_flip() +
    theme(
      legend.text = element_text(size = 15),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 16),
      axis.text.y = element_text(size = 15),
      axis.text.x = element_text(size = 16),
      legend.title = element_blank(),
      axis.title.x = element_blank(),
      plot.title = element_text(size = 20)
    ) +
    labs(title = title)
  
  return(plot)
}

# Perform analysis on pathways (presence/absence)
subsystems_matrix_presence_byGroup <- as.data.frame(subSys_matrix_presence)
subsystems_matrix_presence_byGroup$Group <- color_samples_byHB$Bacteroides[order(color_samples_byHB$SampleName)]
dispersion_pathways <- perform_permanova_betadisper(subsystems_matrix_presence_byGroup, "Group")
#           Df SumOfSqs      R2      F Pr(>F)    
# Group     1 0.033868 0.33846 15.861  0.001 ***
# Residual 31 0.066196 0.66154                  
# Total    32 0.100063 1.00000                                       
# TukeyHSD p=0.1303989 ; permutest p=0.122
# The non-significant results from the betadisper and permutation tests imply that 
# the within-group variability is not significantly different between the HB and NB groups.
# This means that while the overall compositions differ, the consistency of pathway presence 
# within each group is similar.

dispersion_pathways_presence <- generate_dispersion_plot(dispersion_pathways, "Group", "Pathway composition")

# Perform analysis on pathways (abundance)
subsystems_matrix_abund_byGroup <- as.data.frame(subSys_matrix)
subsystems_matrix_abund_byGroup$Group <- color_samples_byHB$Bacteroides[order(color_samples_byHB$SampleName)]
dispersion_pathways_abund <- perform_permanova_betadisper(subsystems_matrix_abund_byGroup, "Group")
#           Df SumOfSqs      R2      F Pr(>F)    
# Group     1  0.16745 0.21893 8.6891  0.001 ***
#  Residual 31  0.59741 0.78107                  
# Total    32  0.76486 1.00000     
# TukeyHSD p=p=0.015 ; permutest p=0.009
# This suggests that there are significant differences in the variability (dispersion) of pathway abundances
# within each group. Specifically, that the spread of the data points 
# from the centroid is significantly different between the HB and NB groups.

dispersion_pathways_abund_plot <- generate_dispersion_plot(dispersion_pathways_abund, "Group", "Pathway abundance")

# Combine plots
plot1 <- dispersion_pathways_presence + coord_flip() +
  theme(axis.title.x = element_text(size = 16))
plot2 <- dispersion_pathways_abund_plot + coord_flip() + 
  theme(axis.title.x = element_text(size = 16))
cowplot::plot_grid(plot1, plot2, ncol = 2)


# Perform analysis on reactions (presence/absence)
reactions_matrix_presence_byGroup <- as.data.frame(t(reaction_matrix_presence))
reactions_matrix_presence_byGroup$Group <- color_samples_byHB$Bacteroides[order(color_samples_byHB$SampleName)]
dispersion_reacs_presence <- perform_permanova_betadisper(reactions_matrix_presence_byGroup, "Group")
#           Df SumOfSqs      R2      F Pr(>F)    
# Group     1 0.015274 0.07642 2.5651  0.071 .
# Residual 31 0.184590 0.92358                
# Total    32 0.199864 1.00000                             
# Bacteroides presence/absence does not play a role in the overall reaction variety (based on presence/absence)
# TukeyHSD p=p=0.442 ; permutest p=0.452

dispersion_reacs_presence_plot <- generate_dispersion_plot(dispersion_reacs_presence, "Group", "Reaction composition")

# Perform analysis on reactions (abundance)
reactions_matrix_abund_byGroup <- as.data.frame(reaction_matrix_weightedAbund)
reactions_matrix_abund_byGroup$Group <- color_samples_byHB$Bacteroides[order(color_samples_byHB$SampleName)]
dispersion_reacs_abund <- perform_permanova_betadisper(reactions_matrix_abund_byGroup, "Group")
#           Df SumOfSqs      R2      F Pr(>F)    
# Group     1  0.16596 0.23117 9.3212  0.001 ***
# Residual 31  0.55194 0.76883                  
# Total    32  0.71790 1.00000  
# TukeyHSD p=p=0.0059 ; permutest p=0.004

dispersion_reacs_abund_plot <- generate_dispersion_plot(dispersion_reacs_abund, "Group", "Reaction abundance")

# Combine and save reaction plots
plot3 <- dispersion_reacs_abund_plot + coord_flip() +
  theme(axis.title.x = element_text(size = 16))
plot4 <- dispersion_reacs_presence_plot + coord_flip() +
  theme(axis.title.x = element_text(size = 16))
cowplot::plot_grid(plot3, plot4, ncol = 2)


# Perform analysis on taxonomic composition (at family level)
table_all_samples_abund_family <- table_all_samples_abund %>%
  group_by(sample, tax_family) %>%
  summarize(relative_abund_fam = sum(relative_abund), .groups = 'drop') %>%
  mutate(Group = color_samples_byHB$Bacteroides[match(sample, color_samples_byHB$SampleName)]) %>%
  pivot_wider(names_from = "tax_family", values_from = "relative_abund_fam", values_fill = 0)

table_all_samples_abund_family_wide <- as.data.frame(table_all_samples_abund_family)
rownames(table_all_samples_abund_family_wide) <- table_all_samples_abund_family_wide$sample

# Run PERMANOVA and dispersion analysis for taxonomic composition
dispersion_tax_family <- perform_permanova_betadisper(table_all_samples_abund_family_wide, "Group")
#           Df SumOfSqs      R2      F Pr(>F)    
# Group     1   1.6642 0.22584 9.0432  0.001 ***
# Residual 31   5.7047 0.77416                  
# Total    32   7.3688 1.00000                  
# TukeyHSD p=0.0041 ; permutest p=0.002

dispersion_taxfamily_plot <- generate_dispersion_plot(dispersion_tax_family, "Group", "Taxonomic composition at the family level")

