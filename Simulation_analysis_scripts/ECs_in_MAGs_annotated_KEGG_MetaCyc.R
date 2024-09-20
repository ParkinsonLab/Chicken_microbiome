library(dplyr)
library(reshape2)
library(ggplot2)


table_all_samples_abund <- read.table("ModelGeneration_Files/Day7_50_10_bins_info_merged_111indAssebmly_116pooled_5coassembly_5EcoliClusters_min0.005.txt",
                                      sep = '\t',header=T)

################################################################################################################
###########     FIGURE SUPPLEMENTARY: ECs in MAGs annotated by KEGG/MetaCyc pathways          ##################
################################################################################################################
# table generated in python via script Find_uniqueECs_in_XMLmodels.py
ECs_models <- read.delim("Table_ECnumbers_inEachModel_gapfilled_corn_absorbed_2fibre_May2024.txt", sep='\t',header=T)
ECs_models$Specie <- unlist(lapply(ECs_models$SpecieName, function(x) unique(Day7_out_table_pooled_all$tax_specie[Day7_out_table_pooled_all$binID_new == x])))
ECs_models$tax_family <- unlist(lapply(ECs_models$SpecieName, function(x) unique(Day7_out_table_pooled_all$tax_family[Day7_out_table_pooled_all$binID_new == x]))) 

# Split the EC_numbers into separate rows
expanded_df <- ECs_models %>%
  separate_rows(EC_numbers, sep = ", ") %>%
  mutate(EC_numbers = trimws(EC_numbers)) # Remove any leading/trailing whitespace

# Create a binary presence/absence matrix
binary_matrix <- expanded_df %>%
  pivot_wider(names_from = EC_numbers, values_from = EC_numbers, values_fill = list(EC_numbers = 0), values_fn = list(EC_numbers = length)) %>%
  mutate(across(-c(Specie, SpecieName, tax_family), ~ ifelse(. > 0, 1, 0)))

# Set the row names to the Specie column and remove the Specie column
binary_matrix <- as.data.frame(binary_matrix)
rownames(binary_matrix) <- binary_matrix$SpecieName

# replace with Specie for all except "Ecoli_cluster", and fix duplicated Faecenecus_gallistercoris:
rownames2 <- binary_matrix$Specie
rownames2[grep('Ecoli_cluster',rownames(binary_matrix))] <- rownames(binary_matrix)[grep('Ecoli_cluster',rownames(binary_matrix))]
rownames2[grep('Faecenecus_gallistercoris',rownames2)] <- paste0(rownames2[grep('Faecenecus_gallistercoris',rownames2)],"_sp",
                                                                 c(1:length(rownames2[grep('Faecenecus_gallistercoris',rownames2)] )))
rownames(binary_matrix)  <- rownames2

### now read the keeg pwy info table:
kegg_pwy <- read.delim("kegg_pwy_info.tbl")

### now read the metacyc pwy info table:
metacyc_pwy <- read.delim("meta_pwy.tbl.txt")


# Extract EC numbers from the binary matrix column names
ec_numbers <- colnames(binary_matrix)[-c(1:3)] # excluding the first column which is Specie, SPecie_name and tax_family

# Function to get name and hierarchy for EC numbers
extract_ec_info <- function(data, ec_numbers) {
  ec_info <- data %>%
    separate_rows(reaEc, sep = ",") %>%
    mutate(reaEc = trimws(reaEc)) %>%
    filter(reaEc %in% ec_numbers) %>%
    group_by(reaEc) %>%
    summarise(
      name = list(name),
      hierarchy = list(hierarchy)
    )
  return(ec_info)
}

# Apply the function to the pathway data
ec_info_kegg <- extract_ec_info(kegg_pwy, ec_numbers)
ec_info_meta <- extract_ec_info(metacyc_pwy, ec_numbers)


# Unnest the lists into separate rows
ec_info_kegg_unnested <- ec_info_kegg %>%
  unnest(cols = c(name, hierarchy))
ec_info_meta_unnested <- ec_info_meta %>%
  unnest(cols = c(name, hierarchy))

# Collapse the table to keep only unique rows based on reaEc and hierarchy
unique_ec_info_kegg <- ec_info_kegg_unnested %>%
  distinct(reaEc, hierarchy, .keep_all = TRUE)

unique_ec_info_meta <- ec_info_meta_unnested %>%
  distinct(reaEc, hierarchy, .keep_all = TRUE)

# add ECs that are in metaCyc but not in Kegg:
meta_ecs <- setdiff(unique(ec_info_meta$reaEc), unique(ec_info_kegg$reaEc))
unique_ec_info_meta <- unique_ec_info_meta[which(unique_ec_info_meta$reaEc %in% meta_ecs),]

# fix meta names
unique_ec_info_meta$Pathway <- NA
unique_ec_info_meta$Pathway[grep("Lipid-biosynthesis|Lipid-Biosynthesis|Lipid-Degradation",unique_ec_info_meta$hierarchy)] <- "Lipid metabolism"
unique_ec_info_meta$Pathway[grep("Amino-Acid-Degradation",unique_ec_info_meta$hierarchy)] <- "Amino acid metabolism"
unique_ec_info_meta$Pathway[grep("Amino-Acid-Biosynthesis",unique_ec_info_meta$hierarchy)] <- "Amino acid metabolism"
unique_ec_info_meta$Pathway[grep("Glycan",unique_ec_info_meta$hierarchy)] <- "Glycan biosynthesis and metabolism"
unique_ec_info_meta$Pathway[grep("Energy-Metabolism",unique_ec_info_meta$hierarchy)] <- "Energy metabolism"
unique_ec_info_meta$Pathway[grep("Carbohydrates-Degradation",unique_ec_info_meta$hierarchy)] <- "Carbohydrate metabolism"
unique_ec_info_meta$Pathway[grep("Carbohydrates-Biosynthesis",unique_ec_info_meta$hierarchy)] <- "Carbohydrate metabolism"
unique_ec_info_meta$Pathway[grep("CARBOXYLATES-DEG",unique_ec_info_meta$hierarchy)] <- "Carbohydrate metabolism"
unique_ec_info_meta$Pathway[grep("Cofactor",unique_ec_info_meta$hierarchy)] <- "Metabolism of cofactors and vitamins"
unique_ec_info_meta$Pathway[grep("AROMATIC",unique_ec_info_meta$hierarchy)] <- "Xenobiotics biodegradation and metabolism"

# remove remaining NA's
unique_ec_info_meta <- unique_ec_info_meta[which(!(is.na(unique_ec_info_meta$Pathway))),]

unique_ec_info_meta <- unique_ec_info_meta %>%
  distinct(reaEc, Pathway, .keep_all = TRUE)

# add kegg names
unique_ec_info_kegg$Pathway <- gsub("kegg;Metabolism;|kegg;Genetic Information Processing;|kegg;Environmental Information Processing;","", unique_ec_info_kegg$hierarchy)
# rbind kegg and metacyc:
unique_ec_info_kegg_meta <- rbind(unique_ec_info_kegg, unique_ec_info_meta)

# Identify EC numbers with multiple annotations
multiple_annotations <- unique_ec_info_kegg_meta %>%
  group_by(reaEc) %>%
  filter(n() > 1) %>%
  ungroup()

# Remove rows with "kegg;Metabolism;Global and overview maps" if there are other annotations for the same EC number
filtered_ec_info <- multiple_annotations %>%
  filter(!(reaEc %in% reaEc[hierarchy == "kegg;Metabolism;Global and overview maps"] & hierarchy == "kegg;Metabolism;Global and overview maps"))

# Combine the filtered multiple annotations with single annotations
unique_ec_info_woGlobal <- unique_ec_info_kegg_meta %>%
  filter(!reaEc %in% multiple_annotations$reaEc) %>%
  bind_rows(filtered_ec_info)

# Function to create Pathway_category1 and Pathway_category2
create_pathway_category <- function(pathways) {
  unique_pathways <- unique(pathways)
  if (length(unique_pathways) == 1) {
    return(data.frame(Pathway_category1 = unique_pathways, Pathway_category2 = NA))
  } else if (length(unique_pathways) == 2) {
    sorted_pathways <- sort(unique_pathways)
    return(data.frame(Pathway_category1 = sorted_pathways[1], Pathway_category2 = sorted_pathways[2]))
  } else {
    return(data.frame(Pathway_category1 = "Multiple superpathways", Pathway_category2 = NA))
  }
}

# Add the Pathway_category columns
unique_ec_info_woGlobal <- unique_ec_info_woGlobal %>%
  group_by(reaEc) %>%
  do(create_pathway_category(.$Pathway)) %>%
  right_join(unique_ec_info_woGlobal, by = "reaEc") %>%
  ungroup()



# Ensure all ECs are accounted for in unique_ec_info_woGlobal
all_ecs <- colnames(binary_matrix)[-c(1:3)]  # Exclude tax_family & Specie & SpecieName column

# Create a data frame for all ECs, filling missing ones with "Unknown"
all_ec_info <- tibble(
  reaEc = all_ecs
) %>%
  left_join(unique_ec_info_woGlobal, by = "reaEc") %>%
  mutate(
    Pathway_category1 = ifelse(is.na(Pathway_category1), "Unknown", Pathway_category1),
    Pathway_category2 = ifelse(is.na(Pathway_category2), "Unknown", Pathway_category2)
  )

all_ec_info_unduplicated <- all_ec_info[,c(1:3)]
all_ec_info_unduplicated <- all_ec_info_unduplicated %>%
  distinct()

# Define the order of pathways
pathway_order <- c("Multiple superpathways", "Carbohydrate metabolism", "Energy metabolism", "Amino acid metabolism",
                   "Metabolism of other amino acids"  ,"Lipid metabolism","Metabolism of terpenoids and polyketides",
                   "Metabolism of cofactors and vitamins" ,"Glycan biosynthesis and metabolism" ,"Xenobiotics biodegradation and metabolism",
                   "Nucleotide metabolism", "Biosynthesis of other secondary metabolites", "Global and overview maps",
                   "Translation","Signal transduction" ,"Unknown")



# Convert Pathway_category1 to a factor with the specified order
all_ec_info_unduplicated <- all_ec_info_unduplicated %>%
  mutate(
    Pathway_category1 = factor(Pathway_category1, levels = pathway_order)
  )

# Reorder the all_ec_info to match the column names of the binary_matrix and reorder columns based on Pathway_category1
all_ec_info_unduplicated <- all_ec_info_unduplicated %>%
  arrange(Pathway_category1, match(reaEc, colnames(binary_matrix)[-c(1:3)]))

# Reorder the columns in binary_matrix based on Pathway_category1
ordered_ecs <- all_ec_info_unduplicated$reaEc

binary_matrix <- binary_matrix %>%
  select(tax_family, SpecieName, all_of(ordered_ecs))

all_ec_info_unduplicated$Pathway_category1 <- as.character(all_ec_info_unduplicated$Pathway_category1)

# Combine all unique pathway categories from Pathway_category1
all_unique_categories <- unique(c(all_ec_info_unduplicated$Pathway_category1, all_ec_info_unduplicated$Pathway_category2))
all_unique_categories <- all_unique_categories[!is.na(all_unique_categories)]  # Remove NA


# make color palette for pathway categories:
# Use Spectral palette and adjust the number of colors if there are more unique categories
# pathway_colors <- setNames(brewer.pal(n = min(length(all_unique_categories), 11), "Spectral"), all_unique_categories)

# Ensure colors are not repeated if more than 11 unique categories
if(length(all_unique_categories) > 11) {
  pathway_colors <- setNames(colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(length(all_unique_categories)), all_unique_categories)
}

pathway_colors["Multiple superpathways"] <- "black"
pathway_colors["Unknown"] <- "white"

# Create column annotation for Pathway_category1 and Pathway_category2 with unified color mapping
column_anno <- HeatmapAnnotation(
  Pathway_category1 = all_ec_info_unduplicated$Pathway_category1,
  Pathway_category2 = all_ec_info_unduplicated$Pathway_category2,
  col = list(
    Pathway_category1 = pathway_colors,
    Pathway_category2 = pathway_colors
  ),
  annotation_legend_param = list(
    Pathway_category1 = list(title = "Pathway Category 1", at = pathway_order, labels = pathway_order),
    Pathway_category2 = list(title = "Pathway Category 2", at = pathway_order, labels = pathway_order)
  ), show_annotation_name = F
)

# Create row annotation for tax_family
row_anno <- rowAnnotation(
  tax_family = binary_matrix$tax_family,
  col = list(tax_family = families_ordered_to_plot_colors), show_annotation_name = F
)

presence_in_NB <- unique(Day7_out_table_pooled_all$binID_new[Day7_out_table_pooled_all$sample %in% 
                                                               sample_categories$SampleName[sample_categories$Category == "LB"]])
presence_in_HB <- unique(Day7_out_table_pooled_all$binID_new[Day7_out_table_pooled_all$sample %in% 
                                                               sample_categories$SampleName[sample_categories$Category == "HB"]])

# Create a data frame for presence categories
species_presence <- data.frame(
  Species = unique(Day7_out_table_pooled_all$binID_new),
  Presence = sapply(unique(Day7_out_table_pooled_all$binID_new), function(x) {
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


# Remove the tax_family column from binary_matrix for heatmap plotting
heatmap_data <- binary_matrix %>% select(-tax_family)

# reorder rows in species_presence to match heatmap_data order:
species_presence <- species_presence %>% arrange(match(Species, heatmap_data$SpecieName))
heatmap_data <- heatmap_data[,-1] # remove specieName now

# Create the heatmap
EC_bins <- Heatmap(
  as.matrix(heatmap_data),
  name = "ECs",
  cluster_columns = F,
  show_column_dend = F,
  top_annotation = column_anno,
  right_annotation = row_anno,
  show_row_names = TRUE,
  show_column_names = F,
  row_split = species_presence$Presence,
  cluster_row_slices = T,
  show_heatmap_legend = T,
  col = circlize::colorRamp2(c(0, 1), c("#0f075e",  "#ed1561"))
)


## save the heatmap as excel table for supplementary data:
ht = draw(EC_bins) ## this draws the heatmap
row_order = row_order(ht)
column_order = column_order(ht)
# since weapplied splitting on rows and columns, `row_order` and column_order are lists - need to unlist it
row_order = unlist(row_order)
column_order = unlist(column_order)
# reorder `mat`
mat2 = heatmap_data[row_order, column_order]
# add top annotation:
mat2_topannot <- rbind(as.character(all_ec_info_unduplicated$Pathway_category1),
                       as.character(gsub("Unknown","", all_ec_info_unduplicated$Pathway_category2)),
                       mat2)

write.csv(mat2_topannot, "Table_HeatMap_allECs_KEGGMetaCycannotated_allMAGs.csv")
