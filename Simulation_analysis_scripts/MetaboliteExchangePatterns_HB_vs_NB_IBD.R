library(tidyr)
library(broom)
library(dplyr)
library(ggplot2)
library(reshape2)


seed_metabolites <- readr::read_delim("gapseq_github_dat/all_seed_metabolites_edited.tsv", delim='\t',col_names=T)
seed_metabolites <- as.data.frame(seed_metabolites)

SCFAconc_experimental <- readxl::read_xlsx("D7SCFA_BioSampleIDs.xlsx", sheet = 1)
color_samples_byHB <- SCFAconc_experimental[, c("SampleName", "Bacteroides")]
color_samples_byHB$color <- ifelse(color_samples_byHB$Bacteroides == "HB", "gold", "forestgreen")
sample_categories <- color_samples_byHB[,1:2]


# ==================================================================================================================
#####################################################################################################################
########                        METABOLITE EXCHANGE PATTERNS IN HB VS NB COMMUNITIE                 #################
#####################################################################################################################
# ==================================================================================================================
# function to add names of the metabolites
addMetNames_flux <- function(df) {
  df$MetName <- NA
  for (met in unique(df$Reaction)) {
    if (is.na(met) | met == "EX_cpd11416_c0") {next}
    df$MetName[df$Reaction == met] <- seed_metabolites$name[seed_metabolites$id == gsub('EX_|_e0','',met)]
  }
  return(df)
}


# load absolute abundances of species at each hour of simulation (for all sample-replicates)
AbsAbund <-  read.table(paste0("AbsAbund_tables_allSamples_corn_absorbed_2fibre_allReps.txt"),
                        sep = '\t', header=T)
# summarize to get community sizes at each hour
Commsize_corn_allReps <- AbsAbund %>% group_by(Hour,Sample,Replicate) %>% summarise(CommSize = sum(AbsAbundance))
colnames(Commsize_corn_allReps)[2] <- "SampleName"

# load total fluxes by community at each hour
flux_allMets_byCommunity_allReps_16hours <- 
  read.table(paste0("flux_table_allMets_byCommunity_allReps_corn_absorbed_2fibre_sim13062024.txt"),
             sep = '\t', header=T) %>%
  addMetNames_flux()

# Combine with community sizes data and normalize flux by community size at each hour
flux_total_allMets_16hours_normalizedByCommSize <- flux_allMets_byCommunity_allReps_16hours %>%
  left_join(Commsize_corn_allReps, by = c("Hour", "SampleName", "Replicate")) %>%
  mutate(Flux_norm = Flux / CommSize)


# Sum all the fluxes over 16h
flux_total_allMets_normalizedByCommSize <- flux_total_allMets_16hours_normalizedByCommSize %>% 
  group_by(SampleName, MetName, Replicate) %>% 
  summarise(TotalFlux_norm = sum(Flux_norm)) %>%
  ungroup() %>%
  # get mean by replicates:
  group_by(SampleName,MetName) %>%
  summarise(TotalFlux_norm_mean = mean(TotalFlux_norm))

# Assign NB/HB category and prepare the data with MetName and Category as factors
df_flux <- flux_total_allMets_normalizedByCommSize %>%
  left_join(sample_categories, by = "SampleName") %>%
  mutate(
    MetName = as.factor(MetName),
    Category = as.factor(Category)
  )

# Calculate mean fluxes and pivot to wide format
mean_fluxes <- df_flux %>%
  group_by(MetName, Category) %>%
  summarise(MeanFlux = mean(TotalFlux_norm_mean, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Category, values_from = MeanFlux, names_prefix = "Flux_mean_") %>%
  # remove rows where both Flux_mean_HB and Flux_mean_NB = 0:
  filter(!(Flux_mean_HB == 0 & Flux_mean_NB == 0))


# Running the Wilcoxon rank-sum test for each MetName
results <- df_flux %>%
  group_by(MetName) %>%
  # Ensure that there are at least two different categories to compare
  filter(length(unique(Category)) > 1) %>%
  # Perform the test comparing categories within each MetName group
  do(broom::tidy(wilcox.test(TotalFlux_norm_mean ~ Category, data = ., exact = FALSE))) %>%
  ungroup() %>%
  # Select only the columns of interest
  select(MetName, statistic, p.value)

results$p.adjusted <- p.adjust(results$p.value, method = "BH")

# Combine Wilcoxon results with mean fluxes and calculate fold change
final_results <- results %>%
  left_join(mean_fluxes, by = "MetName") %>%
  filter(!is.na(p.adjusted)) %>%
  mutate(
    Flux_FoldChange = Flux_mean_HB / Flux_mean_NB,
    Flux_log2FC = log2(abs(Flux_FoldChange)),
    # Adjust the sign of log2FC based on fold change direction (for plotting, to make consistent)
    Flux_log2FC = case_when(
      Flux_FoldChange < -1 ~ -Flux_log2FC,
      Flux_FoldChange > -1 & Flux_FoldChange < 0 & Flux_mean_HB > 0 ~ abs(Flux_log2FC),
      TRUE ~ Flux_log2FC
    )
  )

# add categories (to assign colors for plotting)
final_results <- final_results %>%
  mutate(ChangeType = case_when(
    Flux_mean_HB > 0 & Flux_mean_NB > 0 & Flux_log2FC > 0  ~ "Increased production in HB",
    Flux_mean_HB > 0 & Flux_mean_NB > 0 & Flux_log2FC < 0  ~ "Increased production in NB",
    Flux_mean_HB < 0 & Flux_mean_NB < 0 & Flux_log2FC > 0  ~ "Increased consumption in NB",
    Flux_mean_HB < 0 & Flux_mean_NB < 0 & Flux_log2FC < 0  ~ "Increased consumption in HB",
    Flux_mean_HB > 0 & Flux_mean_NB < 0                   ~ "Switch from consumption in NB to production in HB",
    Flux_mean_HB < 0 & Flux_mean_NB > 0                   ~ "Switch from production in NB to consumption in HB",
    TRUE                                                   ~ "Other"
  ))

final_results$ChangeType[final_results$p.adjusted > 0.05] <- "non-significant"

# Plotting
Vlc_plot <- ggplot(final_results, aes(x = Flux_log2FC, y = -log10(p.adjusted), label = MetName)) +
  geom_point(aes(color = ChangeType), alpha = 0.75, size = 1.2) +
  ggrepel::geom_text_repel(aes(label = ifelse(p.adjusted < 0.05 & abs(Flux_log2FC) > 0.65, as.character(MetName), ""),
                               color = ChangeType), 
                           size = 6, 
                           max.overlaps = Inf) +
  scale_x_continuous(limits = c(-6,6), breaks = c(-5,-3,-1,1,3,5)) +
  # xlim(-6,6) +
  scale_color_manual(values = c("Switch from consumption in NB to production in HB" = "darkred",
                                "Switch from production in NB to consumption in HB" = "magenta",
                                "Increased production in HB" = "forestgreen",
                                "Increased production in NB" = "darkblue",
                                "Increased consumption in HB" = "purple",
                                "Increased consumption in NB" = "orange",
                                "non-significant" = "black")) +
  theme_bw() +
  labs(x = "Log2(Flux FoldChange)", y = "-Log10(p-value adjusted)",
       # title = "Volcano plot of metabolic fluxes",
       color = "Change Type") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  theme(legend.position = "right",
        axis.text = element_text(size=15),
        axis.title = element_text(size=17))



