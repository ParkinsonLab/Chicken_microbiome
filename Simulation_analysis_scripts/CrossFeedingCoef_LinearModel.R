library(ggplot2)
library(scales)
library(tidyverse)
library(dplyr)
library(tidyr)
library(vegan)
library(circlize)


SCFAconc_experimental <- readxl::read_xlsx("D7SCFA_BioSampleIDs.xlsx", sheet = 1)
color_samples_byHB <- SCFAconc_experimental[, c("SampleName", "Bacteroides")]
color_samples_byHB$color <- ifelse(color_samples_byHB$Bacteroides == "HB", "gold", "forestgreen")
sample_categories <- color_samples_byHB[,1:2]

table_all_samples_abund <- read.table("Day7_50_10_bins_info_merged_111indAssebmly_116pooled_5coassembly_5EcoliClusters_min0.005.txt",
                                      sep = '\t',header=T)


########################################################################################
################     CROSS-FEEDING COEFFICIENT AND LINEAR MODEL      ###################
########################################################################################

# load the dataframe with all the cross-feeding pairs and metabolites exchanged in each sample 
# (generated in the script CrossFeeding_interactions_heatmap.R)
common_interactions <- read.table("Table_common_interactions_byCategory_sim130524.txt", sep='\t',header=T)

### recalculate number of total interactions but for each sample
# treat producer and consumer as an unordered pair 
# (e.g. specieA->specieB is equal to specieB->specieA to count both directions of cross-feeding)
producer_consumer_summary_2directions_samples <- common_interactions %>%
  mutate(Species_Pair = pmin(Producer, Consumer), # Get the minimum of Producer and Consumer lexicographically
         Species_Pair_Other = pmax(Producer, Consumer)) %>% # Get the maximum of Producer and Consumer lexicographically
  group_by(SampleName, Category, Species_Pair, Species_Pair_Other) %>%
  summarise(Num_Metabolites_Exchanged = n_distinct(MetName), .groups = 'drop') %>%
  arrange(Category, Species_Pair, Species_Pair_Other)

interaction_summary_samples_all <- producer_consumer_summary_2directions_samples %>%
  group_by(SampleName, Category, Species_Pair = pmin(Species_Pair, Species_Pair_Other), 
           Species_Pair_Other = pmax(Species_Pair, Species_Pair_Other)) %>%
  summarise(Total_Metabolites_Exchanged = sum(Num_Metabolites_Exchanged), .groups = 'drop')

# Calculate the total interactions (sum of exchanged metabolites, in both roles) for each sample
total_interactions_samples_all <- interaction_summary_samples_all %>%
  mutate(Species = Species_Pair) %>%
  bind_rows(interaction_summary_samples_all %>%
              mutate(Species = Species_Pair_Other)) %>%
  group_by(SampleName) %>%
  summarise(Total_Interactions = sum(Total_Metabolites_Exchanged)) %>%
  arrange(desc(Total_Interactions))  # arrange to prioritize samples with more interactions

## merge with CMD scores:
# load CMD scores per sample (generated in the script CMDS_MRO_uniqueECs_calculation.R)
sample_cmds_scores <- read.table("CommunityMetDiversityScores_log_adjusted_33samples_corn2.txt",
                                 sep ='\t',header=T)
# combine tables
total_interactions_samples_cmds <- merge(total_interactions_samples_all, sample_cmds_scores)

C_n_k = function(n, k) {
  factorial(n) / factorial(n-k) / factorial(k)
}

## add cross-feeding coeff = total cross-feeding pairs / C(n,2)
total_interactions_samples_cmds$CrossFeeding_coef <- 
  unlist(lapply(1:nrow(total_interactions_samples_cmds), function(row) {
    total_interactions_samples_cmds$Total_Interactions[row] / 
      C_n_k(n=length(ECs_samples_models[[total_interactions_samples_cmds$SampleName[row]]]), k=2)
  }))

summary_stats_CrossFeeding <- total_interactions_samples_cmds %>%
  group_by(Category) %>%
  summarize(
    Mean = mean(CrossFeeding_coef),
    SD = sd(CrossFeeding_coef),
    ymin = Mean - SD,
    ymax = Mean + SD )

wilcox.test(total_interactions_samples_cmds$CrossFeeding_coef[total_interactions_samples_cmds$Category == "NB"],
            total_interactions_samples_cmds$CrossFeeding_coef[total_interactions_samples_cmds$Category == "HB"]) 
# p-value = 0.02729 

# Plot CrossFeeding_coef vs category (NB/HB)
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



####################################################
################ BUILD LINEAR MODEL
####################################################

# Function to check the presence of a specific species in the sample models
check_presence <- function(sample_names, species) {
  sapply(sample_names, function(sample) {
    species %in% names(ECs_samples_models[[sample]])
  }) %>% as.integer()
}

# Add variables presence/absence of  L. crispatus, and A. butyraticus
total_interactions_samples_cmds <- total_interactions_samples_cmds %>%
  mutate(
    Lcrispatus_presence = check_presence(SampleName, "Lactobacillus_crispatus"),
    Abutyraticus_presence = check_presence(SampleName, "Anaerostipes_butyraticus")
  )

# Create E. coli category based on its abundance
Ecoli_category <- table_all_samples_abund %>%
  filter(tax_specie == "Escherichia_coli") %>%
  mutate(
    Ecoli_category = case_when(
      relative_abund < 0.05 ~ "LowEcoli",
      relative_abund >= 0.05 ~ "HighEcoli",
      TRUE ~ "NoEcoli"
    )
  ) %>%
  select(sample, Ecoli_category)

# Merge E. coli categories into the main dataframe
total_interactions_samples_cmds <- total_interactions_samples_cmds %>%
  left_join(Ecoli_category, by = c("SampleName" = "sample")) %>%
  mutate(Ecoli_presence_category = ifelse(is.na(Ecoli_category), "NoEcoli", Ecoli_category)) %>%
  select(-Ecoli_category)

# Fit the linear model
model_cf3 <- lm(
  CrossFeeding_coef ~ Ecoli_presence_category:(Category + Abutyraticus_presence + Lcrispatus_presence) + comm_size,
  data = total_interactions_samples_cmds)
summary(model_cf3)

#  Interaction plot for Category (HB or NB) and Ecoli_presence
ggplot(total_interactions_samples_cmds, 
       aes(x=Category, y=CrossFeeding_coef, color=factor(Ecoli_presence_category,
                                                         levels = c("NoEcoli","LowEcoli","HighEcoli"))))+
  geom_point() +
  scale_color_manual(values=c("gray84", "lightblue","red"), 
                     labels = c("E.coli absent", "E.coli low abundance (<5%)", "E.coli high abundance (>5%)")) +
  stat_summary(fun=mean, geom="line", aes(group=Ecoli_presence_category)) +
  labs(title="Effect of E. coli Presence Across Categories
on cross-feeding activity") + 
  theme_classic() + 
  theme(axis.text = element_text(size=14) ,
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=16),
        legend.title = element_blank())

# add variable describing simultaneous presence/absence of E. coli and A. butyraticus
total_interactions_samples_cmds$Interaction_Ecoli_Abutyraticus <- 
  interaction(total_interactions_samples_cmds$Ecoli_presence_category, total_interactions_samples_cmds$Abutyraticus_presence, sep = " & ")

# Interaction plot showing the effect of E. coli and A. butyraticus on CrossFeeding_coef
interactionPlot <- ggplot(total_interactions_samples_cmds, 
                          aes(x = Category, y = CrossFeeding_coef, color = Interaction_Ecoli_Abutyraticus,
                              group = Interaction_Ecoli_Abutyraticus)) +
  geom_line(aes(linetype = Interaction_Ecoli_Abutyraticus), 
            stat = "summary", fun = "mean") +
  geom_point(stat = "summary", fun = "mean", size = 3, shape = 21, fill = "white") +
  labs(x = "Community Category (NB vs. HB)", 
       y = "CrossFeeding Coefficient", 
       title = "Interaction Effects of E.coli and A.butyraticus
on cross-feeding activity") +
  scale_color_manual(values = c("red", "blue", "gray", "purple","orange","cyan"),
                     labels = c("High E.coli & A.butyraticus Absent", 
                                "Low E.coli & A.butyraticus Absent", 
                                "No E.coli & A.butyraticus Absent", 
                                "High E.coli & A.butyraticus Present",
                                "Low E.coli & A.butyraticus Present",
                                "No E.coli & A.butyraticus Present")) +
  theme_bw() + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 11),
        axis.text = element_text(size=14) ,
        axis.title.y = element_text(size=16),
        axis.title.x = element_blank())


