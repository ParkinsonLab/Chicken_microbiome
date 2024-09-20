library(ggplot2)
library(scales)
library(tidyverse)
library(dplyr)
library(tidyr)
library(vegan)
library(ComplexHeatmap)
library(circlize)

seed_metabolites <- readr::read_delim("gapseq_github_dat/all_seed_metabolites_edited.tsv", delim='\t',col_names=T)
seed_metabolites <- as.data.frame(seed_metabolites)

SCFAconc_experimental <- readxl::read_xlsx("D7SCFA_BioSampleIDs.xlsx", sheet = 1)
color_samples_byHB <- SCFAconc_experimental[, c("SampleName", "Bacteroides")]
color_samples_byHB$color <- ifelse(color_samples_byHB$Bacteroides == "HB", "gold", "forestgreen")
sample_categories <- color_samples_byHB[,1:2]


####################################################################################
#######        IDENTIFY KEY METABOLITES EXCHANGED BY SPECIES       #################
###############    WITH CONTRIBUTION BY KEYSTONE TAXA        #######################
####################################################################################
# function to add names of the metabolites
addMetNames_flux <- function(df) {
  df$MetName <- NA
  for (met in unique(df$Reaction)) {
    if (is.na(met) | met == "EX_cpd11416_c0") {next}
    df$MetName[df$Reaction == met] <- seed_metabolites$name[seed_metabolites$id == gsub('EX_|_e0','',met)]
  }
  return(df)
}

# load the medium formulation
medium_corn2 <- read.table("ModelGeneration_Files/corn_absorbed_2fibre_diet_grower_aggregated_2023.txt",
                           sep='\t',header=T)  %>%
  select(Reaction = met_id, Flux = met_mM) %>%
  addMetNames_flux() %>%
  mutate(
    Specie = "Medium",
    Replicate = 1
  )

# multiply by 2e-04 (the proportion of the medium added per hour) and 16 hours :
medium_corn2$Flux <- 2e-04 * 16 * medium_corn2$Flux

# and convert to fmol/arena - as total influx of the metabolite with medium:
# *625 - convert mM (in diet table) to mM/cell (will be fmol/cell in arena - 10^12 times higher) because:
# "arena grid size [cm]: 0.025 x 0.025";  "1 mM in arena correspons to mmol/grid_cell: 6.25e-10"
# 'mM'    = { conv <- 10^12 * 0.01 * object@scale}, # conversion of mMol in fmol/grid_cell = 10^12*0.01*6.25e-08 = 625
# fmol/arena = fmol/grid cell * (object@n*object@m) = 625*10^4

medium_corn2$Flux <- 625*10^4*medium_corn2$Flux

# Load essential metabolites and adjust columns similarly to the medium
essential_mets <- read.table("ModelGeneration_Files/EssentialMetsTable_corn_absorbed_2.txt", sep = '\t') %>%
  dplyr::select(Reaction = V1, SampleName = V2) %>%
  dplyr::mutate(
    Flux = 3e-06 * 625 * 10^4, # 3e-06 - concentration at which essential metabolites are added at timepoint=1
    Specie = "Medium", 
    Replicate = 1
  ) %>%
  addMetNames_flux()

# Function to repeat dataframe with new SampleName column
repeat_dataframe_with_new_col <- function(df, new_col_values) {
  df[rep(seq_len(nrow(df)), each = length(new_col_values)), ] %>%
    mutate(SampleName = rep(new_col_values, times = nrow(df)))
}

# Apply the function to repeat medium data for each sample
medium_corn2_ext <- repeat_dataframe_with_new_col(medium_corn2, sample_categories$SampleName)

# Combine medium data with essential metabolites
medium_corn2_ext_ess <- bind_rows(medium_corn2_ext, essential_mets)

# Load all the exchange fluxes by specie
flux_allMets_bySpecie <- read.table("flux_table_sum_allMets_bySpecie_allReps_corn_absorbed_2fibre_sim13052024.txt",
                                    sep = '\t', header=T) %>% 
  addMetNames_flux()

# Add medium to fluxes and join with sample categories
flux_allMets_bySpecie_corn2_wMedium <- bind_rows(flux_allMets_bySpecie_corn2, medium_corn2_ext_ess) %>%
  left_join(sample_categories, by = "SampleName")

# Calculate mean flux by replicate and then total fluxes by category
total_flux_allMets_corn2_wMedium_byCategory_meanByRep <- flux_allMets_bySpecie_corn2_wMedium %>%
  group_by(SampleName, Specie, MetName, Category) %>%
  summarise(Flux = mean(Flux), .groups = "drop") %>%
  group_by(MetName, Category) %>%
  summarise(
    Total_Consumption = sum(Flux[Flux < 0]),
    Total_Production = sum(Flux[Flux > 0]),
    .groups = "drop")

# Identify producers and consumers
producers_corn2_meanByRep <- flux_allMets_bySpecie_corn2_wMedium_meanByRep %>%
  filter(Flux > 0) %>%
  dplyr::rename(Producer = Specie, Producer_Flux = Flux)

consumers_corn2_meanByRep <- flux_allMets_bySpecie_corn2_wMedium_meanByRep %>%
  filter(Flux < 0) %>%
  dplyr::rename(Consumer = Specie, Consumer_Flux = Flux)

cross_feeding_corn2_meanByRep <- merge(producers_corn2_meanByRep, consumers_corn2_meanByRep, 
                                       by = c("SampleName","MetName", "Category"))

# calculate average contribution within category (all samples summarised)
contribution_producers_byCategory <- producers_corn2_meanByRep %>%
  group_by(MetName, Producer, Category) %>%
  summarise(Producer_Flux_total = sum(Producer_Flux)) %>%
  ungroup() %>%
  left_join(total_flux_allMets_corn2_wMedium_byCategory_meanByRep, by = c("MetName", "Category")) %>%
  group_by(MetName, Producer, Category) %>%
  summarise(Production_Contribution = Producer_Flux_total / Total_Production) 

contribution_consumers_byCategory <- consumers_corn2_meanByRep %>%
  group_by(MetName, Consumer, Category) %>%
  summarise(Consumer_Flux_total = sum(Consumer_Flux)) %>%
  ungroup() %>%
  left_join(total_flux_allMets_corn2_wMedium_byCategory_meanByRep, by = c("MetName", "Category")) %>%
  group_by(MetName, Consumer, Category) %>%
  summarise(Consumption_Contribution = Consumer_Flux_total / Total_Consumption,
            Consumption_Contribution_ofProduction = abs(Consumer_Flux_total) / Total_Production) 


# Define key taxa involved in cross-feeding
key_hubs <- c("Escherichia_coli", "Bacteroides_fragilis", "Anaerostipes_butyraticus", "Lactobacillus_crispatus")

# Filter metabolites produced by bacteria (not in the medium) and consumed at least at 10%
contribution_producers_byCategory_bac0.3_key <- contribution_producers_byCategory %>%
  # Filter out unwanted metabolites and join with consumers data
  filter(!(MetName %in% c("O2", medium_corn2_ext_ess$MetName))) %>%
  left_join(contribution_consumers_byCategory, by = c("MetName", "Category")) %>%
  # Further filter by conditions for production and consumption contributions:
  # 1) there is at least one producer per metabolite per Category that contributes > 33% of production
  # 2) at least 1% of the produced metabolite is being consumed by other members
  # 3) or one of the major consumers is keystone specie
  filter(
    (Production_Contribution > 0.3 & sum(Consumption_Contribution_ofProduction) > 0.01 & Producer %in% key_hubs) |
      (Consumption_Contribution_ofProduction > 0.3 & Consumer %in% key_hubs)
  ) %>%
  group_by(MetName, Category)

# Retrieve these metabolites:
mets_bac0.3_keystone <- unique(contribution_producers_byCategory_bac0.3_key$MetName)

### Select the key metabolites:
# 1. Rank the producers within each MetName and Category, select first top3 producers
ranked_producers <- contribution_producers_byCategory[which(contribution_producers_byCategory$MetName %in% mets_bac0.3_keystone),] %>%
  group_by(MetName, Category) %>%
  arrange(desc(Production_Contribution)) %>%
  mutate(rank = row_number()) %>%
  filter(rank < 4)

# 2. Extract the top 3 producers for each MetName and Category
top_producers <- ranked_producers %>%
  select(MetName, Category, Producer, rank)

# 3. Identify metabolites where the top1 producers in NB and HB are different
# as we are interested in the niche adjustments when the keystone species are absent
same_top_producers <- top_producers %>%
  tidyr::spread(Category, Producer) %>%
  filter(NB[rank==1] == HB[rank==1]) %>%
  select(MetName)

# also filter out mets where top1 producer is also top1 consumer:
top_producers_are_top_consumers <- contribution_consumers_byCategory[which(
  contribution_consumers_byCategory$MetName %in% mets_bac0.3_keystone),] %>%
  group_by(MetName, Category) %>%
  arrange(desc(Consumption_Contribution)) %>%
  mutate(rank = row_number()) %>% #  rank consumers
  filter(rank < 4) %>% 
  left_join(top_producers,by = c("MetName", "Category","rank")) %>%  # join with top_producers data
  filter(Consumer == Producer) %>%
  select(MetName)

# filter out metabolites with the same top producer across categories or if producer and consumer is the same specie
mets_bac0.3_keystone_filtered <- 
  mets_bac0.3_keystone[-which(mets_bac0.3_keystone %in% 
                                c(unique(same_top_producers$MetName), top_producers_are_top_consumers$MetName))]

# remove sugars (produced in very low quantities, by 1-2 species in up to 3 samples only): "D-Gluconate", "D-Ribose",  "D-Arabinose":
# and leave only D-Lactate (remove L-Lactate)
mets_bac0.3_keystone_filtered <- 
  mets_bac0.3_keystone_filtered[-which(mets_bac0.3_keystone_filtered %in% c("D-Gluconate", "D-Ribose", "L-Lactate", "D-Arabinose"))]

# Filter producers and consumers for selected metabolites
contribution_producers_byCategory_selmets <- contribution_producers_byCategory %>%
  filter(MetName %in% mets_bac0.3_keystone_filtered)

contribution_consumers_byCategory_selmets <- contribution_consumers_byCategory %>%
  filter(MetName %in% mets_bac0.3_keystone_filtered)

# Get names of top producers and consumers by contribution
top_contributing_producers_names <- contribution_producers_byCategory_selmets %>%
  group_by(MetName, Category) %>%
  arrange(desc(Production_Contribution)) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  distinct(Producer)

top_contributing_consumers_names <- contribution_consumers_byCategory_selmets %>%
  group_by(MetName, Category) %>%
  arrange(desc(Consumption_Contribution)) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  distinct(Consumer)

# Create dataframe for Sankey plot generation with all the top consumer-producer fluxes
top_contributing_crossfeeders_selmets <- contribution_producers_byCategory_selmets %>%
  filter(Producer %in% top_contributing_producers_names$Producer) %>%
  inner_join(
    contribution_consumers_byCategory_selmets %>%
      filter(Consumer %in% top_contributing_consumers_names$Consumer),
    by = c("MetName", "Category")
  ) %>%
  # Remove rows with low production and consumption contributions
  filter(
    Production_Contribution > 0.05,
    Consumption_Contribution > 0.03
  ) %>%
  # Manually remove species with only small contributions that are not critical for display
  filter(
    !Producer %in% c("UBA1417_sp900552925", "Clostridium_Q_saccharolyticum_A", 
                     "SRR22360819_Lachnospiraceae_bin.20", "Ruthenibacterium_avium"),
    !Consumer %in% c("Borkfalkia_sp944368795", "Mediterraneibacter_stercoravium", 
                     "Mediterraneibacter_caccogallinarum")
  )


########## GENERATE SANKEY PLOT
library(networkD3)

# Define unique nodes including metabolites
nodes_split <- top_contributing_crossfeeders_selmets %>%
  select(Producer, Consumer, MetName, Category) %>%
  pivot_longer(cols = c(Producer, Consumer, MetName), names_to = "Role", values_to = "Name") %>%
  mutate(Role = case_when(
    Role == "Producer" ~ "Producer",
    Role == "Consumer" ~ "Consumer",
    TRUE ~ "Metabolite"
  ),
  Unique_Name = paste(Name, Category, Role, sep = "_")) %>%
  select(Unique_Name, Category) %>%
  distinct() %>%
  arrange(Unique_Name) %>%
  mutate(ID = row_number() - 1)  # Ensure IDs are zero-indexed

# Prepare links for producers to metabolites and metabolites to consumers
links_split <- bind_rows(
  # Links from producers to metabolites
  top_contributing_crossfeeders_selmets %>%
    mutate(
      Producer_Name = paste(Producer, Category, "Producer", sep = "_"),
      Metabolite_Name = paste(MetName, Category, "Metabolite", sep = "_")
    ) %>%
    group_by(Producer_Name, Metabolite_Name, Category) %>%
    summarise(Value = first(Production_Contribution), .groups = 'drop') %>%
    mutate(
      Source = nodes_split$ID[match(Producer_Name, nodes_split$Unique_Name)],
      Target = nodes_split$ID[match(Metabolite_Name, nodes_split$Unique_Name)]
    ) %>%
    drop_na(Source, Target),
  
  # Links from metabolites to consumers
  top_contributing_crossfeeders_selmets %>%
    mutate(
      Metabolite_Name = paste(MetName, Category, "Metabolite", sep = "_"),
      Consumer_Name = paste(Consumer, Category, "Consumer", sep = "_")
    ) %>%
    group_by(Metabolite_Name, Consumer_Name, Category) %>%
    summarise(Value = first(Consumption_Contribution), .groups = 'drop') %>%
    mutate(
      Source = nodes_split$ID[match(Metabolite_Name, nodes_split$Unique_Name)],
      Target = nodes_split$ID[match(Consumer_Name, nodes_split$Unique_Name)]
    ) %>%
    drop_na(Source, Target)
)

# Define colors for categories
categories <- unique(nodes_split$Category)
colors <- setNames(c("gold", "forestgreen"), categories)

# Generating the Sankey plot
sankey <- sankeyNetwork(
  Links = links_split,
  Nodes = nodes_split,
  Source = "Source",
  Target = "Target",
  Value = "Value",
  NodeID = "Unique_Name",
  NodeGroup = "Category",
  LinkGroup = "Category",
  colourScale = JS(sprintf(
    'd3.scaleOrdinal().domain(["%s"]).range(["%s"])',
    paste(categories, collapse='","'), 
    paste(colors[categories], collapse='","')
  ))
)

# Print the plot
print(sankey)




