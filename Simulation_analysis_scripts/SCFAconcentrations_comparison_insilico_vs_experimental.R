library(dplyr)
library(reshape2)
library(data.table)
library(ggplot2)

conc_SCFA_corn_absorbed_2fibre <- read.table("conc_table_48h_SCFA_corn_absorbed_2fibre_diet_grower_sim20240513.txt",
                                             sep = '\t', header=T)
# add names for SCFA:
conc_SCFA_corn_absorbed_2fibre <- conc_SCFA_corn_absorbed_2fibre %>% 
  filter(Metabolite != "EX_cpd00159_e0") %>% # remove lactate
  mutate(MetName = case_when(
      Metabolite == "EX_cpd00141_e0" ~ "Propionate",
      Metabolite == "EX_cpd00211_e0" ~ "Butyrate",
      TRUE ~ "Acetate"  # Default to "Acetate" for remaining metabolite
    )
  )

# subset to hour=16
conc_SCFA_corn_absorbed_2fibre_h16 <- conc_SCFA_corn_absorbed_2fibre[which(conc_SCFA_corn_absorbed_2fibre$Hour == 16),]

#load experimental SCFA concentration data and reshape
SCFAconc_experimental <- readxl::read_xlsx("D7SCFA_BioSampleIDs.xlsx", sheet = 1) %>%
  select(Acetate, Propionate, Butyrate, SampleName) %>%
  melt(id.vars = "SampleName", variable.name = "MetName", value.name = "Conc.mean") %>%
  # add columns to merge with the in silico concentrations data
  mutate(
    medium = "experiment",
    Conc.sd = 0,
    Hour = 16
  )


# merge in silico concentrations with experimental
conc_SCFA_corn_insilico_exp <- rbind(conc_SCFA_corn_absorbed_2fibre_h16[,c("Hour","medium","SampleName","Conc.mean","Conc.sd","MetName")], 
                                     SCFAconc_experimental_subset_melt[,c("Hour","medium","SampleName","Conc.mean","Conc.sd","MetName")])


conc_SCFA_corn_insilico_exp_wide <- conc_SCFA_corn_insilico_exp %>% 
  reshape2::dcast(formula = SampleName + MetName ~ medium, value.var = "Conc.mean")
conc_SCFA_corn_insilico_exp_wide$MetName <- factor(conc_SCFA_corn_insilico_exp_wide$MetName,
                                                   levels = c("Acetate","Propionate","Butyrate"))
# plot scatter plot with fitted line and correlation
SCFA_corr_scatter <- ggpubr::ggscatter(conc_SCFA_corn_insilico_exp_wide,
                                       x = "experiment", y = "corn_absorbed_2fibre", 
                                       add = "reg.line", # add.params = list(color = "darkred", fill = "gray84"),
                                       conf.int = TRUE, 
                                       color = "MetName", palette = "lancet") +
  ggpubr::stat_cor(aes(color = MetName), method = "spearman") + 
  ggpubr::theme_classic2()+
  labs(y = "Predicted concentration", x = "Measured concentration") + 
  facet_wrap(MetName ~ . , scales = "free") + 
  theme(axis.title = element_text(size=17),
        axis.text = element_text(size=13),
        strip.text = element_text(size=15),
        legend.position = "none")



# Update Bacteroides category directly within the dataframe
color_samples_byHB <- SCFAconc_experimental[, c("SampleName", "Bacteroides", "Func_capacity")]
sample_categories <- data.frame(SampleName = color_samples_byHB$SampleName, Category = color_samples_byHB$Bacteroides)

conc_SCFA_corn_insilico_exp <- conc_SCFA_corn_insilico_exp %>%
  mutate(Bacteroides = ifelse(SampleName %in% sample_categories$SampleName[sample_categories$Category == "NB"], "NB", "HB"))

# Plot barplots with measured and experimental HB vs NB 
scfa_plots <- list()
for (met in c("Acetate","Propionate","Butyrate")) {
  df <- conc_SCFA_corn_insilico_exp[which(conc_SCFA_corn_insilico_exp$MetName == met),]
  
  # Scaling factor
  if (met %in% c("Propionate")) { 
    sf <- 500
    color_values = c("#ed0000","#f27963")
  } 
  if (met %in% c("Butyrate")) {
    sf <- 90
    color_values = c("#42b540", "#7ee6ae")
  }
  if (met %in% c("Acetate")) {
    sf <- 50
    color_values = c("#00468b", "#8b74d6")
  }
  
  
  # multiply the experimental value by sf:
  df$Conc.mean[df$medium == "experiment"]  <- sf*df$Conc.mean[df$medium == "experiment"]
  
  # Calculate summary statistics for error bars
  df_summary <- df %>%
    group_by(Bacteroides, medium) %>%
    summarise(Conc_mean = mean(Conc.mean), Conc_sd = sd(Conc.mean)) %>%
    mutate(ymin = Conc_mean - Conc_sd, ymax = Conc_mean + Conc_sd)
  
  scfa_plots[[met]] <- ggplot(df_summary, aes(x = factor(Bacteroides, levels = c("NB","HB")), y = Conc_mean, fill = medium, group = interaction(Bacteroides, medium))) + 
    geom_col(position = position_dodge(width = 0.9), width = 0.85) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), position = position_dodge(width = 0.9), width = 0.1) +
    scale_fill_manual(values = color_values) +
    scale_y_continuous(name = "Predicted concentration", labels = scales::comma,
                       sec.axis = sec_axis(~./sf, name = "Measured concentration", labels = scales::comma)) +
    labs(title = met) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 16),
          axis.title.x = element_blank(),
          legend.title = element_blank(),
          legend.position = "bottom",
          title = element_text(size = 19))
  
}

## run wilcox.test:
wilcox_LBvsHB_scfa <- conc_SCFA_corn_insilico_exp %>%
  # run the Wilcoxon for each metabolite within each 'medium' (experiment and simulation)
  group_by(medium, MetName) %>%
  summarise(
    pvalue = wilcox.test(Conc.mean[Bacteroides == "LB"], Conc.mean[Bacteroides == "HB"],
                         exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  # adjust p-values across the three MetNames *within* each medium
  group_by(medium) %>%
  mutate(
    qvalue = p.adjust(pvalue, method = "BH")
  ) %>%
  ungroup()

print(wilcox_LBvsHB_scfa)
#  medium               MetName        pvalue     qvalue
# corn_absorbed_2fibre Acetate    0.0197     0.0296    
# corn_absorbed_2fibre Butyrate   0.255      0.255     
# corn_absorbed_2fibre Propionate 0.00000116 0.00000347
# experiment           Acetate    0.00787    0.00787   
# experiment           Butyrate   0.00361    0.00787   
# experiment           Propionate 0.00633    0.00787 
