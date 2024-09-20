library(dplyr)
library(reshape2)

# load the chicken feed corn diet
diet_corn_grower <- read.table("ModelGeneration_Files/corn_diet_grower_aggregated_2023.txt", sep='\t',header=T)
# Create a copy and remove the 'value' column
diet_corn_grower_cecum <- diet_corn_grower
diet_corn_grower_cecum$value <- NULL

##########################################################################################################################
##########     Adjust concentrations to imitate absorption/breakdown of metabolites before reaching the cecum   ##########
##########################################################################################################################

### 0) Minerals (ref: https://www.animbiosci.org/journal/view.php?doi=10.5713/ajas.2010.90129) 
# Adjust percentages of intake amounts for specific minerals
# Ca - 40%; Phosphorus - 47%; Potassium - 34%; Sodium - 66%; Magnesium - 26%; Iron - 21%; Manganese - 11%; Zinc - 10%, Copper - 8%
mineral_names <- c("Calcium", "Phosphorus", "Potassium", "Sodium", "Magnesium", "Fe3+", "Manganese", "Zinc", "Copper")
mineral_factors <- c(0.4, 0.47, 0.34, 0.66, 0.26, 0.21, 0.11, 0.1, 0.08)
diet_corn_grower_cecum$met_mM[diet_corn_grower_cecum$met_name %in% mineral_names] <- 
  mineral_factors * diet_corn_grower_cecum$met_mM[diet_corn_grower_cecum$met_name %in% mineral_names]

###  1) Carbohydrates  (ref - https://doi.org/10.1093/jn/110.1.117, https://doi.org/10.1016/j.psj.2020.08.072
###   and https://doi.org/10.3382/ps/pev244)
# retain 3% of starches; add the products of their degradation products (malto-..-oses) as 2% of inital concentration:
starch_ids <- c("EX_cpd90003_e0", "EX_cpd90004_e0")
diet_corn_grower_cecum$met_mM[diet_corn_grower_cecum$met_id %in% starch_ids] <- 
  0.03 * diet_corn_grower_cecum$met_mM[diet_corn_grower_cecum$met_id %in% starch_ids]

# Calculate 2% of starches for degradation products
conc_2perc_starch19 <- 0.02 * diet_corn_grower$met_mM[diet_corn_grower$met_id == "EX_cpd90004_e0"]
conc_2perc_starch27 <- 0.02 * diet_corn_grower$met_mM[diet_corn_grower$met_id == "EX_cpd90003_e0"]

# Add malto-oses (products of starch degradation)
malto_oses_corn <- cbind(c("Maltohexaose", "Maltooctaose", "Maltoundecaose", "Maltotridecaose"),
                          c(rep("g_per100g",4)),
                          c("Maltohexaose", "Maltooctaose", "Maltoundecaose", "Maltotridecaose"),
                          c("EX_cpd90007_e0","EX_cpd90005_e0","EX_cpd90006_e0","EX_cpd90008_e0"),
                          c(conc_2perc_starch19_corn,conc_2perc_starch27_corn,conc_2perc_starch19_corn,conc_2perc_starch27_corn),
                          c(rep("Carbohydrates",4)))
colnames(malto_oses_corn) <- colnames(diet_wheat_grower_cecum)

diet_corn_grower_cecum <- rbind(diet_corn_grower_cecum, malto_oses_corn)

# Remove glucose and fructose completely, retain 2% of galactose, raffinose, sucrose ,maltose and stachyose:
diet_corn_grower_cecum$met_mM <- as.numeric(diet_corn_grower_cecum$met_mM)
diet_corn_grower_cecum$met_mM[diet_corn_grower_cecum$met_name %in% c("D-Glucose", "D-Fructose")] <- 0
diet_corn_grower_cecum$met_mM[diet_corn_grower_cecum$met_name %in%  c("Galactose", "Melitose","Sucrose", "Maltose","Stachyose")] <- 
  0.02 * diet_corn_grower$met_mM[diet_corn_grower$met_name %in% c("Galactose", "Melitose","Sucrose", "Maltose","Stachyose")]


###  2) Proteins (amino acids) (ref - https://doi.org/10.1093/jn/117.8.1459, https://doi.org/10.1093/jn/107.10.1775)
# By the time diet reaches ileum, 80% of ingested amino acids are absorbed
# Assuming 15 more % absorbed in ileum, retain 5% of the dietary AAs:
diet_corn_grower_cecum$met_mM[diet_corn_grower_cecum$category == "Proteins"] <-
  0.05 * diet_corn_grower$met_mM[diet_corn_grower$category == "Proteins"]

###  3) Fats (ref -  https://doi.org/10.3382/ps.2013-03344,  https://doi.org/10.3382/ps/pey458)
# Retain 15% of all fats:
diet_corn_grower_cecum$met_mM[diet_corn_grower_cecum$category == "Fats"] <-
  0.15 * diet_corn_grower$met_mM[diet_corn_grower$category == "Fats"]

### 4) Fibre (cellulose and celloheptaose)
# retain 25% of celloheptaose and cellulose:
diet_corn_grower_cecum$met_mM[diet_corn_grower_cecum$met_name == "cellulose"] <- 
  0.125 * diet_corn_grower$met_mM[diet_corn_grower$met_name == "cellulose"]

cellulose_conc <- diet_corn_grower_cecum$met_mM[diet_corn_grower_cecum$met_name == "cellulose"]
diet_corn_grower_cecum <- rbind(diet_corn_grower_cecum,
                                data.frame(met_name = "Celloheptaose", unit = "g_per100g",
                                           met_id = "EX_cpd03718_e0", met_mM = cellulose_conc,
                                           category = "Carbohydrates")
)

# Based on shadow prices, increase concentration of Riboflavin and Folate by 3x
diet_corn_grower_cecum$met_mM[diet_corn_grower_cecum$met_name %in% c("Riboflavin", "Folate")] <- 
  3 * diet_corn_grower_cecum$met_mM[diet_corn_grower_cecum$met_name %in% c("Riboflavin", "Folate")]

# Save tables for input:
write.table(diet_corn_grower_cecum, "ModelGeneration_Files/corn_absorbed_2fibre_diet_grower_aggregated_2023.txt",
            sep = '\t',quote=F)


########################################################################################################
#################        Generate diet file for gapfilling step (gapseq)      ##########################
########################################################################################################
## Generate media file for gapseq gap-filling steps:
diet_file <- diet_corn_grower_cecum[, c("met_id", "met_name", "met_mM")]
colnames(diet_file) <- c("compounds","name","maxFlux")

# add mucins and urea:
essential_metabolites <- c("EX_cpd00122_e0", "EX_cpd00232_e0","EX_cpd00832_e0",
                           "EX_cpd02992_e0", "EX_cpd11842_e0", "EX_cpd21520_e0","EX_cpd00073_e0")
essential_conc <- 1e-02
mucins_urea_to_add <- cbind("compounds" = essential_metabolites, 
                            "name" = c("acgam","acnam","acgal","core4","core3","Tn_antigen","urea"),
                            "maxFlux" = essential_conc)

diet_file <- rbind(diet_file,mucins_urea_to_add)

diet_file$compounds <- gsub("EX_|_e0","", diet_file$compounds)
diet_file$name[diet_file$compounds == "cpd90003"] <- "starch_n19"
diet_file$name[diet_file$compounds == "cpd90004"] <- "starch_n27"

write.csv(diet_file, paste0( "ModelGeneration_Files/",
                             diet, "_forGapfill_v2.csv"), quote=F,row.names=F)


