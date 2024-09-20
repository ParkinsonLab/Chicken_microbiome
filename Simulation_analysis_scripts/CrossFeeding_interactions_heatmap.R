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

table_all_samples_abund <- read.table("Day7_50_10_bins_info_merged_111indAssebmly_116pooled_5coassembly_5EcoliClusters_min0.005.txt",
                                      sep = '\t',header=T)
### add family/order to plot annotations:
table_all_samples_abund$tax_order <- unlist(lapply(table_all_samples_abund$classification,
                                                     function(x) { strsplit(strsplit(x, "o__")[[1]][2],";")[[1]][1] }))

families_ordered_to_plot <- c()
for (order in unique(table_all_samples_abund$tax_order)) {
  families_ordered_to_plot <- c(families_ordered_to_plot,
                                unique(table_all_samples_abund$tax_family[table_all_samples_abund$tax_order == order] ))
}

families_ordered_to_plot_colors <- c("#7f3667","#a53e76","#bb437e","#d24787", # Oscillospirales - dusty pink shades
                                     "#0f4b00","#146d00","#1aa000","#4ecb61", # Christensenellales - green shades
                                     "#9966cc","#564787","#351c75", # Lachnospirales - purple shades
                                     "#909090", # Actinomycetales - grey
                                     "black", # RFN20
                                     "#754f12", "#c49869", # TANB77 - brown shades
                                     "#700f6b", # RF 39 (UBA660) - deep violet purple
                                     "#ff8200", # orange - Enterobac
                                     "#616b22", # CAJFEE01 - dark olive
                                     "#c0f1ab", # Campylobacterales - light ugly green
                                     "#ff0000", "#bd0000", # Lactobacillales - bright red
                                     "#fff897","#ffc900", # Erysipelotrichales - yellow/gold shades
                                     "#ff8c7f", # Acholeplasmatales- salmon/peach
                                     "#73919c", # Coriobacteriales - dustry grey-ish lightblue
                                     "#cce2ff", "#a5abed", "#6b74d0","#3742b5", # Bacteroidales - light sky blue shades
                                     "#9a0000", # Clostridiales - dark red
                                     "#b99b73")   # Mycoplasmatales - olive beige
names(families_ordered_to_plot_colors) <- families_ordered_to_plot

seed_metabolites <- readr::read_delim("gapseq_github_dat/all_seed_metabolites_edited.tsv", delim='\t',col_names=T)
seed_metabolites <- as.data.frame(seed_metabolites)


#################################################################################
# -------------------------------------------------------------------------------
# ---------------            CROSS-FEEDING HEATMAP           --------------------
# -------------------------------------------------------------------------------
#################################################################################
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

# filter out by total sum of fluxes over 16h > 0.1 (< -0.1)
flux_allMets_bySpecie <- flux_allMets_bySpecie %>%
  filter(abs(Flux) >= 0.1)

# Identify producers and consumers
producers <- flux_allMets_bySpecie_4diets %>%
  filter(Flux > 0) %>%
  dplyr::select(SampleName, Replicate, Specie, Flux, MetName) %>%
  dplyr::rename(Producer = Specie, Producer_Flux = Flux)

consumers <- flux_allMets_bySpecie_4diets %>%
  filter(Flux < 0) %>%
  dplyr::select(SampleName, Replicate, Specie, Flux, MetName) %>%
  dplyr::rename(Consumer = Specie, Consumer_Flux = Flux)

# Merge the producers and consumers dataframes on metabolite and sample information
cross_feeding <- merge(producers, consumers, by = c("SampleName", "Replicate","MetName"))

common_interactions_byRep <- cross_feeding %>%
  group_by(SampleName,  MetName, Producer, Consumer) %>%
  summarise(Interaction_Count = n()) %>%
  ungroup() %>%
  arrange(desc(Interaction_Count)) %>% 
  mutate(Interaction_Count_byRep = Interaction_Count/5)

sample_categories <- data.frame(SampleName = color_samples_byHB$SampleName, Category = color_samples_byHB$Bacteroides)

# Join this with common_interactions dataframe
common_interactions <- merge(common_interactions_byRep, sample_categories, by = "SampleName")
write.table(common_interactions,"Table_common_interactions_byCategory_sim130524.txt", sep='\t',header=T)

# one-directional count of mets exchanged
producer_consumer_summary <- common_interactions %>%
  group_by(Category,  Producer, Consumer) %>%
  summarise(Num_Metabolites_Exchanged = n_distinct(MetName), .groups = 'drop') %>%
  arrange(Category, Producer, Consumer)

# treat producer and consumer as an unordered pair (e.g. specieA->specieB is equal to specieB->specieA to count both directions of cross-feeding)
producer_consumer_summary_2directions <- common_interactions %>%
  mutate(Species_Pair = pmin(Producer, Consumer), # Get the minimum of Producer and Consumer lexicographically
         Species_Pair_Other = pmax(Producer, Consumer)) %>% # Get the maximum of Producer and Consumer lexicographically
  group_by(Category, Species_Pair, Species_Pair_Other, SampleName) %>%
  summarise(Num_Metabolites_Exchanged = n_distinct(MetName), .groups = 'drop') %>%
  arrange(Category, Species_Pair, Species_Pair_Other)


### normalize # of metabolites exchanged by the number of samples where each cross-feeding pair is formed:
producer_consumer_summary_2directions_normByNSamples <- 
  producer_consumer_summary_2directions %>% 
  group_by(Category, Species_Pair, Species_Pair_Other) %>%
  summarise(Num_Metabolites_Exchanged_meanBySample = mean(Num_Metabolites_Exchanged))

producer_consumer_summary_2directions_normByNSamples_HB <-
  producer_consumer_summary_2directions_normByNSamples[producer_consumer_summary_2directions_normByNSamples$Category =="HB",]
producer_consumer_summary_2directions_normByNSamples_NB <-
  producer_consumer_summary_2directions_normByNSamples[producer_consumer_summary_2directions_normByNSamples$Category =="NB",]

### get the top cross-feeding species - that have at least 10 interactions with 12+ distinct metabolites:
species_min12x10mets_HB <- c()
for (specie in unique(producer_consumer_summary_2directions_normByNSamples_HB$Species_Pair)) {
  if (length(which(producer_consumer_summary_2directions_normByNSamples_HB$Num_Metabolites_Exchanged_meanBySample[
    producer_consumer_summary_2directions_normByNSamples_HB$Species_Pair == specie
  ] >= 12)) >= 10) {
    species_min12x10mets_HB <- c(species_min12x10mets_HB,specie)
  }
}

species_min12x10mets_NB <- c()
for (specie in unique(producer_consumer_summary_2directions_normByNSamples_NB$Species_Pair)) {
  if (length(which(producer_consumer_summary_2directions_normByNSamples_NB$Num_Metabolites_Exchanged_meanBySample[
    producer_consumer_summary_2directions_normByNSamples_NB$Species_Pair == specie
  ] >= 12)) >= 10) {
    species_min12x10mets_NB <- c(species_min12x10mets_NB,specie)
  }
}

species_min12x10mets_both <- unique(c(species_min12x10mets_HB,species_min12x10mets_NB))
write.table("87species_min10interactionsXmin12mets_corn2_sim130524.txt", sep='\t', quote=F)

producer_consumer_summary_2directions_normByNSamples_min12x10mets <-
  producer_consumer_summary_2directions_normByNSamples[producer_consumer_summary_2directions_normByNSamples$Species_Pair %in% species_min12x10mets_both &
                                                               producer_consumer_summary_2directions_normByNSamples$Species_Pair_Other %in% species_min12x10mets_both, ]

interaction_summary_top_crossfeeders <- producer_consumer_summary_2directions_normByNSamples_min12x10mets %>%
  group_by(Category, Species_Pair = pmin(Species_Pair, Species_Pair_Other), 
           Species_Pair_Other = pmax(Species_Pair, Species_Pair_Other)) %>%
  summarise(Total_Metabolites_Exchanged = sum(Num_Metabolites_Exchanged_meanBySample), .groups = 'drop')



##########################################################################################
#######     GENERATE A HEATMAP WITH INTERACTIONS OF TOP CROSS-FEEDING TAXA     ###########
##########################################################################################

# Create a comprehensive list of all species appearing in any role
all_species <- unique(c(interaction_summary_top_crossfeeders$Species_Pair, interaction_summary_top_crossfeeders$Species_Pair_Other))

# Determine presence in categories
presence_in_NB <- unique(interaction_summary_top_crossfeeders %>% 
                           filter(Category == "NB") %>%
                           dplyr::select(Species_Pair, Species_Pair_Other) %>%
                           unlist() %>%
                           as.character())

presence_in_HB <- unique(interaction_summary_top_crossfeeders %>% 
                           filter(Category == "HB") %>%
                           dplyr::select(Species_Pair, Species_Pair_Other) %>%
                           unlist() %>%
                           as.character())

# Calculate the total interactions (sum of both roles) for each species
total_interactions <- interaction_summary_top_crossfeeders %>%
  mutate(Species = Species_Pair) %>%
  bind_rows(interaction_summary_top_crossfeeders %>%
              mutate(Species = Species_Pair_Other)) %>%
  group_by(Species) %>%
  summarise(Total_Interactions = sum(Total_Metabolites_Exchanged)) %>%
  arrange(desc(Total_Interactions))  # arrange to prioritize species with more interactions


# Step 2: Create a sorted list of species based on total interactions
sorted_species <- total_interactions$Species

# Create a data frame for presence categories
species_presence <- data.frame(
  Species = sorted_species,
  Presence = sapply(sorted_species, function(x) {
    if (x %in% presence_in_NB & x %in% presence_in_HB) {
      "Both"
    } else if (x %in% presence_in_NB) {
      "NB Only"
    } else if (x %in% presence_in_HB) {
      "HB Only"
    } else {
      "None"
    }
  })
)

# Order species by presence for better visualization
species_presence <- species_presence %>%
  mutate(Presence = factor(Presence, levels = c("NB Only","Both",  "HB Only"))) %>%
  arrange(Presence)


# Standardize the interaction matrix based on sorted species
standardize_matrix_sorted <- function(data, category, sorted_species) {
  data <- data %>% filter(Category == category)
  matrix <- matrix(0, nrow = length(sorted_species), ncol = length(sorted_species), dimnames = list(sorted_species, sorted_species))
  
  for (i in seq_along(sorted_species)) {
    for (j in seq_along(sorted_species)) {
      if (i != j) {  # Avoid counting self-interactions if not relevant
        species_i <- sorted_species[i]
        species_j <- sorted_species[j]
        sum_exchanges <- sum(data$Total_Metabolites_Exchanged[(data$Species_Pair == species_i & data$Species_Pair_Other == species_j) |
                                                                (data$Species_Pair == species_j & data$Species_Pair_Other == species_i)], na.rm = TRUE)
        matrix[i, j] <- sum_exchanges
      }
    }
  }
  return(matrix)
}

# Reorder the interaction matrix according to this new order
sorted_species <- species_presence$Species

# Create interaction matrices for each category using the new sorted approach
interaction_matrix_NB <- standardize_matrix_sorted(interaction_summary_top_crossfeeders, "NB", sorted_species)
interaction_matrix_HB <- standardize_matrix_sorted(interaction_summary_top_crossfeeders, "HB", sorted_species)

# add tax family, order and colors to species_presence:
species_presence$tax_family <- 
  unlist(lapply(species_presence$Species, function(x) unique(Day7_out_table_pooled_all$tax_family[Day7_out_table_pooled_all$tax_specie == x])))
species_presence$tax_order <- 
  unlist(lapply(species_presence$Species, function(x) unique(Day7_out_table_pooled_all$tax_order[Day7_out_table_pooled_all$tax_specie == x])))
species_presence$color_family <- unlist(lapply(species_presence$tax_family, function(x) families_ordered_to_plot_colors[x] ))

# Combine the matrices
combined_matrix <- cbind(interaction_matrix_NB, interaction_matrix_HB)

# Create the heatmap annotation based on columns
annotation_cols <- data.frame(
  Category = c(rep("NB", ncol(interaction_matrix_NB)), rep("HB", ncol(interaction_matrix_HB))) 
)

ha <- HeatmapAnnotation(df = annotation_cols, col = list(Category = c("NB" = "forestgreen", "HB" = "gold")))

# Identify and filter out columns with all zero values
filtered_combined_matrix <- combined_matrix[, colSums(combined_matrix != 0) > 0]

# Update the annotation for columns accordingly
filtered_annotation_cols <- annotation_cols$Category[colSums(combined_matrix != 0) > 0]

# Create the top annotation for the filtered heatmap
ha_filtered <- HeatmapAnnotation(df = data.frame(Category = filtered_annotation_cols),
                                 col = list(Category = c("NB" = "forestgreen", "HB" = "gold")), 
                                 show_annotation_name = F, show_legend = F)

# Create annotations based on the presence data
presence_annotation <- HeatmapAnnotation(df = data.frame(Presence = species_presence$Presence),
                                         col = list(Presence = c("NB Only" = "darkblue", "HB Only" = "cyan", "Both" = "salmon")),
                                         which = "row",show_annotation_name = F)

# Create annotations based on the tax_family data
tax_family_annotation_row <- HeatmapAnnotation(df = data.frame(Presence = species_presence$Presence,
                                                               tax_family = species_presence$tax_family),
                                               col = list(Presence = c("NB Only" = "darkblue", "HB Only" = "cyan", "Both" = "salmon"),
                                                          tax_family = families_ordered_to_plot_colors),
                                               which = "row",show_annotation_name = F)

specie_tax_family_colnames <- unlist(lapply(colnames(filtered_combined_matrix), function(x) {
  unique(species_presence$tax_family[species_presence$Species == x]) } ))

tax_family_annotation_col <- HeatmapAnnotation(df = data.frame(tax_family = specie_tax_family_colnames),
                                               col = list(tax_family = families_ordered_to_plot_colors),
                                               which = "col",show_annotation_name = F, show_legend = F)

# Generate the heatmap
heatmap_filtered <- Heatmap(filtered_combined_matrix, 
                            name = "Total interactions",
                            top_annotation = ha_filtered,
                            # right_annotation = presence_annotation,
                            right_annotation = tax_family_annotation_row,
                            bottom_annotation = tax_family_annotation_col,
                            show_row_names = TRUE,
                            show_column_names = TRUE,
                            column_split = filtered_annotation_cols,
                            row_split = species_presence$Presence,
                            cluster_row_slices = T,
                            # cluster_rows = F, 
                            cluster_columns = TRUE,
                            # col = circlize::colorRamp2(c(0, 15, 35), c("white", "blue", "red"))) # for 87 species 10x12 mets
                            col = circlize::colorRamp2(c(0, 10, 30), c("white", "blue", "red"))) # for all non-dietary mets

## save the heatmap as excel table for supplementary data:
ht = draw(heatmap_filtered) ## this draws the heatmap
row_order = row_order(ht)
column_order = column_order(ht)
# since we applied splitting on rows and columns, `row_order` and column_order are lists - need to unlist it
row_order = unlist(row_order)
column_order = unlist(column_order)

# reorder `mat`
mat2 = filtered_combined_matrix[row_order, column_order]
mat2 <- as.data.frame(mat2)
# add presence in category column
mat2$Presence = species_presence$Presence[row_order]
# add family column
mat2$Taxonomic_family <- species_presence$tax_family[row_order]
write.csv(mat2, "InteractionHeatMap_normByNSamples_fromClusteredHeatmap_4A.csv")



############################################
#### Correlate interactions with unique ECs in top cross-feeding species
# (load this table generated in the script CMDS_MRO_uniqueECs_calculation.R)
n_unique_ec_perSample_matrix <- read.table("Table_numberUniqueECs_eachSpecie_perSample_5Ecoli.txt",
             sep='\t',header=T)

# subset to only species with min 10 cross-feeding pairs formed of at least 12 metabolites exchanged in each:
n_unique_ec_perSample_matrix_10_12 <- n_unique_ec_perSample_matrix[which(
  n_unique_ec_perSample_matrix$Species %in% species_min10x12mets_both), c("Species","rowMean")]

# Get the number of total interactions for each species in each role (consumer of producer) within each category
interaction_summary_both_all <- producer_consumer_summary_2directions %>%
  group_by(Category, Species_Pair = pmin(Species_Pair, Species_Pair_Other), 
           Species_Pair_Other = pmax(Species_Pair, Species_Pair_Other)) %>%
  summarise(Total_Metabolites_Exchanged = sum(Num_Metabolites_Exchanged), .groups = 'drop')

# Calculate the total metabolites exchanged (in both roles) for each species
total_interactions_all <- interaction_summary_both_all %>%
  mutate(Species = Species_Pair) %>%
  bind_rows(interaction_summary_both_all %>%
              mutate(Species = Species_Pair_Other)) %>%
  group_by(Species) %>%
  summarise(Total_Interactions = sum(Total_Metabolites_Exchanged)) %>%
  arrange(desc(Total_Interactions))  # arrange to prioritize species with more interactions

# normalize by the number of samples where the specie is present:
# first make a table with specie and N samples where it's present:
specie_numSamples <- Day7_out_table_pooled_all %>% 
  group_by(tax_specie) %>% 
  summarise(numSamples = length(unique(sample))) %>%
  mutate(Species = tax_specie)

total_interactions_all_normBy_numSamples <- 
  merge(total_interactions_all,specie_numSamples, by = "Species") %>%
  mutate(Total_Interactions_norm = Total_Interactions/numSamples)

# merge with sum of interactions formed on average by this specie across samples (normalized by N samples where specie is present):
n_unique_ec_perSample_matrix_10_12_interactions_normBy_numSample <- merge(n_unique_ec_perSample_matrix_10_12, 
                                                                          total_interactions_all_normBy_numSamples)


## plot interactions vs unique ECs 
ggpubr::ggscatter(n_unique_ec_perSample_matrix_10_12_interactions_normBy_numSample,
                  x = "rowMean", y = "Total_Interactions_norm", size = 1.5,
                  add = "reg.line", add.params = list(color = "darkred", fill = "gray84"),
                  conf.int = TRUE,  palette = "jco" ) +
  ggpubr::stat_cor(aes(color = "darkblue"), method = "pearson") + 
  ggrepel::geom_text_repel(data = n_unique_ec_perSample_matrix_10_12_interactions_normBy_numSample[which(n_unique_ec_perSample_matrix_10_12_interactions_normBy_numSample$rowMean > 5 &
                                                                                                           n_unique_ec_perSample_matrix_10_12_interactions_normBy_numSample$Total_Interactions_norm > 230),],
                           aes(label=n_unique_ec_perSample_matrix_10_12_interactions_normBy_numSample$Specie[which(n_unique_ec_perSample_matrix_10_12_interactions_normBy_numSample$rowMean > 5 &
                                                                                                                     n_unique_ec_perSample_matrix_10_12_interactions_normBy_numSample$Total_Interactions_norm > 230)],
                               # size= 0.1,
                               vjust=1), 
                           max.overlaps = 40,
                           segment.color="#b3b5b3",
                           seed=14,  
                           force=10,
                           position = ggpp::position_nudge_center(x = 0.3,
                                                                  y = 0.4 ,
                                                                  direction = "split"),
                           size = 3.5) +
  theme_classic()+
  scale_x_log10() +
  labs(x="Average number of unique ECs in a sample", y = "Total interactions (normalized by sample)") +
  theme(axis.title = element_text(size=17),
        axis.text = element_text(size=16),
        legend.position = "none")






summary_stats_CrossFeeding <- total_interactions_samples_all_cmds_mro_score %>%
  group_by(Category) %>%
  summarize(
    Mean = mean(CrossFeeding_coef),
    SD = sd(CrossFeeding_coef),
    ymin = Mean - SD,
    ymax = Mean + SD )

wilcox.test(total_interactions_samples_all_cmds_mro_score$CrossFeeding_coef[total_interactions_samples_all_cmds_mro_score$Category == "NB"],
            total_interactions_samples_all_cmds_mro_score$CrossFeeding_coef[total_interactions_samples_all_cmds_mro_score$Category == "HB"]) 
# p-value = 0.02729 

# Create the plot CrossFeeding_coef vs category (NB/HB)
ggplot(summary_stats_CrossFeeding, aes(x=Category, y=Mean)) + 
  geom_col(fill = "#330273", position = position_dodge()) + 
  geom_errorbar(aes(ymin=ymin, ymax=ymax), width=0.0, position = position_dodge(0.9)) + 
  theme_bw() + 
  coord_flip() +
  labs(y = "Cross-feeding coefficient") +
  theme(axis.text.y = element_text(size=20),
        axis.title.y=element_blank(),
        axis.title = element_text(size=20),
        axis.text.x = element_text(size=18))













