library(vegan)
library(igraph)
library(ggraph)
library(dplyr)
library(ggplot2)
library(tidyr)
library(cowplot)

# lists with unique reactions within samples and with reactions in each model within each sample 
# (generated in the "Reactions_Pathways_content_analyses.R script)
samples_reacs_list <- readRDS("Samples_day7_uniqueReacs_list_corn2_5Ecoli_May2024.RDS")            
samples_models_reactions_list <- readRDS("Samples_day7_allModelsReacs_list_corn2_5Ecoli_May2024.RDS")

table_all_samples_abund <- read.table("ModelGeneration_Files/Day7_50_10_bins_info_merged_111indAssebmly_116pooled_5coassembly_5EcoliClusters_min0.005.txt",
                                      sep = '\t',header=T)

sample_categories <- readxl::read_excel("~/Desktop/Study/UofT/ParkinsonLab/Modeling_chickens/Chicken_cecal_microbiome_project/ModelGeneration_Files/D7SCFA_BioSampleIDs.xlsx", sheet = "D7SCFA") %>%
  dplyr::select("SampleName","Bacteroides")

####################################################################################################
################           GENERATE NETWORKS FOR ALL SAMPLES                    ####################
####################################################################################################

########### First, create a circular network for each sample with species consolidated based on taxonomic families

# Helper function to create a graph for each sample - consolidated into tax families
create_model_graph_family <- function(sample_name, samples_models_reactions_list,
                                      keep_with_unique_function=F,
                                      ColorScheme, scaled_coef=NULL) {
  model_names <- names(samples_models_reactions_list[[sample_name]])
  
  # if keeping only models that have some unique ECs (at least 1) within the sample (cover some unique functionality in a community)
  if (keep_with_unique_function == T) {
    keep_species <- samples_species_uniqueECs[[sample_name]]
    model_names <- keep_species
  }
  
  model_names_families <- table_all_samples_abund[which(table_all_samples_abund$tax_specie %in% model_names),
                                                  c("tax_specie","tax_family")] 
  model_names_families <- model_names_families[!duplicated(model_names_families),]
  
  model_names_families$tax_family <- factor(model_names_families$tax_family,
                                            levels = families_ordered_to_plot[families_ordered_to_plot %in% model_names_families$tax_family])
  
  # consolidate the number of unique reactions within the family:
  model_names_families$N_unique_reacs_family <- NA
  for (tax_family in unique(model_names_families$tax_family)) {
    family_reacs <- c()
    for (specie in model_names_families$tax_specie[model_names_families$tax_family == tax_family]) {
      family_reacs <- unique(c(family_reacs, samples_models_reactions_list[[sample_name]][[specie]]))
    }
    
    model_names_families$N_unique_reacs_family[model_names_families$tax_family == tax_family] <-
      length(family_reacs)
  }
  
  model_families <- model_names_families %>% dplyr::select(tax_family,N_unique_reacs_family) %>% dplyr::distinct()
  
  matrix <- matrix(1, nrow = length(unique(model_families$tax_family)), ncol = length(unique(model_families$tax_family)))
  
  # Create a graph for the models
  g <- graph_from_adjacency_matrix(matrix, mode = "undirected", diag = FALSE)
  
  # Add model names and their colors to the graph
  V(g)$name <- model_families$tax_family
  V(g)$color <- ColorScheme[as.character(model_families$tax_family)]
  V(g)$order <- levels(model_names_families$tax_family)
  
  # Calculate the size for each model based on their reactions
  V(g)$size <- sapply(V(g)$name, function(fam) model_families$N_unique_reacs_family[model_families$tax_family == fam])
  
  # Use a circular layout
  plot <- ggraph(g, layout = 'circle',order = levels(model_names_families$tax_family)) +
    # geom_edge_link(aes(edge_width = weight), color = "gray80") +
    geom_node_point(aes(size = size, color = color)) +
    scale_size_continuous(range = c(30, 70)) + # Adjust the range based on your data
    # scale_edge_width(range = scaled_coef*c(0.1, 2)) + # Adjust the range based on Jaccard index
    scale_color_identity() +
    # scale_color_manual(values = ColorScheme[model_names]) +
    theme_void() +
    coord_fixed() +
    theme(legend.position = "none")
  
  return(plot)
}


# Calculate the size for each sample-node based on unique reactions within a sample
sample_sizes <- sapply(sample_names, function(s) length(samples_reacs_list[[s]]))
sample_sizes_rescaled <- scales::rescale(sample_sizes, to = c(8,20))

# make plots with families (species consolidated)
for (sample_name in sample_names) {
  pdf(paste0("Networks_circular_33samples/Network_",sample_name,
             "_consolidatedFamilies_colorByOrder_noEdges_scaled30_60_corn2_5Ecoli.pdf"),
      width = 0.5*sample_sizes_rescaled[sample_name],height = 0.5*sample_sizes_rescaled[sample_name])
  
  plot <- create_model_graph_family(sample_name,
                                    samples_models_reactions_list,
                                    ColorScheme = families_ordered_to_plot_colors) 
  print(plot)
  dev.off()
}





########################################################################################################################
###########  GENERATE NETWORK WITH 33 SAMPLES (AS NODES) CONNECTED BASED ON MUTUAL INFORMATION SCORE  ##################
########################################################################################################################

##### MUTUAL INFORMATION SCORE   
# Count reactions in each sample and aggregate 

# Combine all samples into a single matrix for analysis
all_reactions <- unique(unlist(samples_reacs_list))
sample_names <- names(samples_reacs_list)
combined_reaction_matrix <- matrix(0, ncol = length(sample_names), nrow = length(all_reactions),
                                   dimnames = list(all_reactions,sample_names))

combined_reaction_matrix_presence <- matrix(0, ncol = length(sample_names), nrow = length(all_reactions), 
                                            dimnames = list(all_reactions,sample_names))

# Fill the matrix with 1/0 for presence/absence of reactions
for (sample in sample_names) {
  reactions <- samples_reacs_list[[sample]]
  combined_reaction_matrix_presence[reactions, sample] <- 1
}

library(infotheo)
# Function to calculate mutual information between columns
calculate_mi <- function(matrix) {
  n <- ncol(matrix)
  mi_matrix <- matrix(0, n, n)
  colnames(mi_matrix) <- colnames(matrix)
  rownames(mi_matrix) <- colnames(matrix)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      mi <- mutinformation(discretize(matrix[,i]), discretize(matrix[,j]))
      mi_matrix[i, j] <- mi
      mi_matrix[j, i] <- mi
    }
  }
  return(mi_matrix)
}

# Calculate mutual information matrix
# mi_matrix <- calculate_mi(combined_reaction_matrix)
mi_matrix_presence <- calculate_mi(combined_reaction_matrix_presence)

# Create graph from mutual information matrix
threshold <- 0.1  # Define a threshold to filter weak connections
MIscore_matrix_reactions_0.1 <- mi_matrix_presence

MIscore_matrix_reactions_0.1[MIscore_matrix_reactions_0.1 < 0.1] <- 0

# Calculate the size for each sample based on unique reactions
sample_sizes <- sapply(sample_names, function(s) length(samples_reacs_list[[s]]))

# Create a graph for the models
g <- graph_from_adjacency_matrix(MIscore_matrix_reactions_0.1, mode = "undirected", weighted = TRUE, diag = FALSE)

# Add model names and their colors to the graph
V(g)$name <- sample_names
V(g)$color <- colors_samples[sample_names]

# Calculate the size for each model based on their reactions
V(g)$size <- sample_sizes

colfunc <- colorRampPalette(c("darkblue", "gray94"))
my.palette <- rev(colfunc(length(E(g)$weight)))

E(g)$color = my.palette[rev(order(E(g)$weight))]

range01 <- function(x){(x-min(x))/(max(x)-min(x))} 

# Plot network
plot_samples_network <- ggraph(g, layout = 'stress') + 
  geom_edge_link(aes(edge_alpha = range01(weight), edge_width = range01(weight), color=weight)) + # Adjust the color based on Jaccard index
  geom_node_point(aes(size = size,alpha = 0.2, color = color)) +
  scale_size_continuous(range = c(10, 50)) + # Adjust the range based on your data
  scale_edge_width(range = c(0.1, 3)) + 
  scale_edge_colour_gradient2(
    low = "gray94",
    mid = "lightblue",
    high = "darkblue",
    midpoint = 0.17,
    space = "Lab",
    guide = "edge_colourbar",
    aesthetics = "edge_colour"
  ) +
  scale_color_identity() +
  geom_node_text(aes(label = name), repel = TRUE) + # Add labels to the nodes  
  theme_void() 

library(RCy3)
cytoscapePing()
RCy3::createNetworkFromIgraph(g,"Network_33samples_MIscore")
#### continue in Cytoscape and then manually overlay the circular sample-specific networks
### from previous step onto the sample-nodes in the MI network




# Function to perform Wilcoxon test between two sets of values
perform_wilcox_test <- function(values1, values2) {
  test_result <- wilcox.test(values1, values2, exact = FALSE)
  return(test_result$p.value)
}

# Melt the MI matrix
mi_df <- melt(mi_matrix_presence, varnames = c("Var1","Var2"), value.name = "value")

# Build sample‐level categories
sample_categories <- sample_categories %>%
  mutate(
    Bacteroides = ifelse(Bacteroides == "NB", "NB", "HB"),
    Ecoli       = ifelse(SampleName %in% table_all_samples_abund$sample[
        table_all_samples_abund$tax_specie == "Escherichia_coli"
      ], "Ecoli", "noEcoli"
    )
  ) %>%
  select(SampleName, Bacteroides, Ecoli)

# Join in both labels
mi_labeled <- mi_df %>%
  left_join(sample_categories, by = c("Var1" = "SampleName")) %>%
  rename(Bact1 = Bacteroides, Ecoli1 = Ecoli) %>%
  left_join(sample_categories, by = c("Var2" = "SampleName")) %>%
  rename(Bact2 = Bacteroides, Ecoli2 = Ecoli)

# Define the categories of interest:
mi_flagged <- mi_labeled %>%
  mutate(
    # broad groups:
    is_HB_HB = (Bact1 == "HB" & Bact2 == "HB"),
    is_NB_NB = (Bact1 == "NB" & Bact2 == "NB"),
    # E. coli subgroups:
    is_HB_Ecoli_HB_Ecoli       = is_HB_HB & Ecoli1 == "Ecoli"   & Ecoli2 == "Ecoli",
    is_HB_noEcoli_HB_noEcoli   = is_HB_HB & Ecoli1 == "noEcoli" & Ecoli2 == "noEcoli",
    is_NB_Ecoli_NB_Ecoli       = is_NB_NB & Ecoli1 == "Ecoli"   & Ecoli2 == "Ecoli",
    is_NB_noEcoli_NB_noEcoli   = is_NB_NB & Ecoli1 == "noEcoli" & Ecoli2 == "noEcoli"
  )

# calculate mean ± SD for each category (flag)
stats_list <- list(
  HB_HB        = filter(mi_flagged, is_HB_HB)$value,
  NB_NB        = filter(mi_flagged, is_NB_NB)$value,
  HB_Ecoli     = filter(mi_flagged, is_HB_Ecoli_HB_Ecoli)$value,
  HB_noEcoli   = filter(mi_flagged, is_HB_noEcoli_HB_noEcoli)$value,
  NB_Ecoli     = filter(mi_flagged, is_NB_Ecoli_NB_Ecoli)$value,
  NB_noEcoli   = filter(mi_flagged, is_NB_noEcoli_NB_noEcoli)$value
) %>%
  lapply(function(x) list(mean = mean(x), sd = sd(x), values = x))


# Perform Wilcoxon tests between categories
wilcox_results <- list(
  HB_NB = perform_wilcox_test(stats_list$HB_HB$values, stats_list$NB_NB$values),
  HB_Ecoli_vs_NB_Ecoli = perform_wilcox_test(stats_list$HB_Ecoli$values, stats_list$NB_Ecoli$values),
  HB_noEcoli_vs_NB_noEcoli = perform_wilcox_test(stats_list$HB_noEcoli$values, stats_list$NB_noEcoli$values),
  HB_Ecoli_vs_NB_noEcoli = perform_wilcox_test(stats_list$HB_Ecoli$values, stats_list$NB_noEcoli$values)
)

# Collect raw p-values into a data.frame
wilcox_df <- data.frame(
  comparison = names(wilcox_results),
  pvalue     = unlist(wilcox_results),
  stringsAsFactors = FALSE
)

# Apply BH correction across all four tests
wilcox_df$p.adjusted <- p.adjust(wilcox_df$pvalue, method = "BH")

# Print a summary
print(wilcox_df)
# comparison                                        pvalue       p.adjusted
# HB_NB                                       HB_NB 1.306710e-01 1.742280e-01
# HB_Ecoli_vs_NB_Ecoli         HB_Ecoli_vs_NB_Ecoli 2.212948e-02 4.425896e-02
# HB_noEcoli_vs_NB_noEcoli HB_noEcoli_vs_NB_noEcoli 1.000000e+00 1.000000e+00
# HB_Ecoli_vs_NB_noEcoli     HB_Ecoli_vs_NB_noEcoli 1.014021e-05 4.056085e-05



