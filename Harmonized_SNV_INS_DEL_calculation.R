# Calculating the harmonized SNV, INS and DEL values that are used as features in the ERCC2-panel based model.
# This is an example code showing how to calculate these features using your own datasets. 

# Load a genomic information file (BED file) for the desired panel(s)
bed_file <- read.table("bed_file.txt", header = TRUE, sep = "\t")

# Calcualte the panel coverage in megabases (Mb)
panel_size <- (sum(unique(bed_file)$End_Position - unique(bed_file)$Start_Position)) / 1e6 

# Calculate the harmonized features from the maf (or maf-like style) dataframe. 
# This is done for each panel version separately and can be merged in one dataframe later on.
# maf_file should contain mutation information only about the samples of interest (after all processing steps as shown in Processing_maf_file.R).
# The resulting dataframe will contain raw number of SNV, INS, DEL, as as well as normalized and standardized to z-score values.

library(bestNormalize)

panel_harmonized_features <- maf_file %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(
    SNV = sum(Variant_Type == "SNP"),
    INS = sum(Variant_Type == "INS"),
    DEL = sum(Variant_Type == "DEL"), 
    TMB_SNV = SNV / panel_size,
    TMB_INS = INS / panel_size,
    TMB_DEL = DEL / panel_size
  ) %>%
  ungroup() %>%  
  mutate(
    TMB_SNV_jittered = TMB_SNV + runif(nrow(.), min = -1e-6, max = 1e-6),
    TMB_INS_jittered = TMB_INS + runif(nrow(.), min = -1e-6, max = 1e-6),
    TMB_DEL_jittered = TMB_DEL + runif(nrow(.), min = -1e-6, max = 1e-6),
    TMB_SNV_transformed = orderNorm(TMB_SNV_jittered)$x.t,
    TMB_INS_transformed = orderNorm(TMB_INS_jittered)$x.t,
    TMB_DEL_transformed = orderNorm(TMB_DEL_jittered)$x.t
  ) %>%
   mutate(
    TMB_SNV_z = (TMB_SNV_transformed - mean(TMB_SNV_transformed, na.rm = TRUE)) / sd(TMB_SNV_transformed, na.rm = TRUE),
    TMB_INS_z = (TMB_INS_transformed - mean(TMB_INS_transformed, na.rm = TRUE)) / sd(TMB_INS_transformed, na.rm = TRUE),
    TMB_DEL_z = (TMB_DEL_transformed - mean(TMB_DEL_transformed, na.rm = TRUE)) / sd(TMB_DEL_transformed, na.rm = TRUE)
  )

# Join the hamronized values dataframes for all panels in one dataframe
panel_mutation_data <- rbind(panel_harmonized_features, panel_harmonized_features_2, panel_harmonized_features_3)

# Add a total_num_mutations column for each sample, calculating the number of SNV+INS+DEL for each sample
panel_mutation_data <- panel_mutation_data %>% 
  group_by(Tumor_Sample_Barcode) %>%
  mutate(total_num_mutations = sum(SNV, INS, DEL)) %>% 
  ungroup()

# Filter out samples with <5 mutations
panel_mutation_data <- panel_mutation_data %>% filter(total_num_mutations >= 5)
