# This function enables the generation of a dataframe with all of the possible features (SNV_TMB_z, INS_TMB_z, DEL_TMB_z, and SBS and ID cosine simialrties with signatures found in bladder cancer)
# This dataframe was used for trainng purposes, specifically for the feature selection process, to identify the features with the highest weight that should be used in the model
# For testing the model, only the features included in the final model can be subseted (SNV_TMB_z, DEL_TMB_z, cosine_sbs5, cosine_sbs2, cosine_sbs13, cosine_id8)
# This code uses dataframes generated using the following scripts: panel_mutation_data (generated in Harmonized_SNV_INS_DEL_calculation.R), sbs_cosine_similarity and id_cosine_similarity (generated in Mutational_matrices_and_cosine_similarity.R)

# Make sure the panel_mutation_data (generated in Harmonized_SNV_INS_DEL_calculation.R) has the sample names are rownames
rownames(panel_mutation_data) <- panel_mutation_data$Tumor_Sample_Barcode 

# Create a function that combines all possible features in one dataframe
# For samples that are missing specific features (ex. samples without deletions or indels would have NA values for all ID cosine similarity values), input NA in the dataframe

calculate_features <- function(sample_idx, 
                     sbs_cosine_similarity_matrix, 
                     id_cosine_similarity_matrix, 
                     tmb_z_score_df,
                     all_samples) {
  
  # Check if sample exists in SBS cosine similarity matrix
  if (sample_idx %in% rownames(sbs_cosine_similarity_matrix)) {
    cosine_sbs5 <- sbs_cosine_similarity_matrix[sample_idx, "SBS5"]  
    cosine_sbs2 <- sbs_cosine_similarity_matrix[sample_idx, "SBS2"]
    cosine_sbs1 <- sbs_cosine_similarity_matrix[sample_idx, "SBS1"]  
    cosine_sbs8 <- sbs_cosine_similarity_matrix[sample_idx, "SBS8"]
    cosine_sbs13 <- sbs_cosine_similarity_matrix[sample_idx, "SBS13"]  
    cosine_sbs29 <- sbs_cosine_similarity_matrix[sample_idx, "SBS29"]
    cosine_sbs40 <- sbs_cosine_similarity_matrix[sample_idx, "SBS40"]  
    total_sbs_count <- subset(tmb_z_score_df, Tumor_Sample_Barcode == sample_idx)$TMB_SNV_z
  } else {
    cosine_sbs5 <- NA
    cosine_sbs2 <- NA
    cosine_sbs1 <- NA
    cosine_sbs8 <- NA
    cosine_sbs13 <- NA
    cosine_sbs29 <- NA
    cosine_sbs40 <- NA
    total_sbs_count <- NA
  }
  
  # Check if sample exists in indel cosine similarity matrix
  if (sample_idx %in% rownames(id_cosine_similarity_matrix)) {
    cosine_id8 <- id_cosine_similarity_matrix[sample_idx, "ID8"] 
    cosine_id4 <- id_cosine_similarity_matrix[sample_idx, "ID4"] 
    cosine_id5 <- id_cosine_similarity_matrix[sample_idx, "ID5"]
    cosine_id9 <- id_cosine_similarity_matrix[sample_idx, "ID9"] 
    cosine_id10 <- id_cosine_similarity_matrix[sample_idx, "ID10"]
    cosine_id1 <- id_cosine_similarity_matrix[sample_idx, "ID1"] 
    cosine_id2 <- id_cosine_similarity_matrix[sample_idx, "ID2"]
    cosine_id3 <- id_cosine_similarity_matrix[sample_idx, "ID3"] 
  } else {
    cosine_id8 <- NA
    cosine_id4 <- NA
    cosine_id5 <- NA
    cosine_id9 <- NA
    cosine_id10 <- NA
    cosine_id1 <- NA
    cosine_id2 <- NA
    cosine_id3 <- NA
  }
  
  if (sample_idx %in% all_samples) {
    total_ins_count <- subset(tmb_z_score_df, Tumor_Sample_Barcode == sample_idx)$TMB_INS_z
    total_del_count <- subset(tmb_z_score_df, Tumor_Sample_Barcode == sample_idx)$TMB_DEL_z
  } else { 
    total_ins_count <- NA
    total_del_count <- NA
  }

  # Combine features into a data frame
  panel_features <- data.frame(
    cosine_sbs5 = cosine_sbs5,  
    cosine_sbs2 = cosine_sbs2,
    cosine_sbs1 = cosine_sbs1,
    cosine_sbs8 = cosine_sbs8,
    cosine_sbs13 = cosine_sbs13,
    cosine_sbs29 = cosine_sbs29,
    cosine_sbs40 = cosine_sbs40,
    cosine_id8 = cosine_id8,
    cosine_id4 = cosine_id4,
    cosine_id5 = cosine_id5,
    cosine_id9 = cosine_id9,
    cosine_id10 = cosine_id10,
    cosine_id1 = cosine_id1,
    cosine_id2 = cosine_id2,
    cosine_id3 = cosine_id3,
    SNV_TMB_z = as.numeric(total_sbs_count),
    INS_TMB_z = as.numeric(total_ins_count),
    DEL_TMB_z = as.numeric(total_del_count),
    stringsAsFactors = FALSE
  )
  
  return(panel_features)
}

all_panel_features_more_sigs <- list()  

for (sample_idx in panel_mutation_data$Tumor_Sample_Barcode) {
  sample_features <- calculate_features_for_sample_panel_more_sigs(
    sample_idx = sample_idx, 
    sbs_cosine_similarity_matrix = sbs_cosine_similarity,
    id_cosine_similarity_matrix = id_cosine_similarity,
    tmb_z_score_df = panel_mutation_data,
    all_samples = panel_mutation_data$Tumor_Sample_Barcode
  )
  all_panel_features_more_sigs[[sample_idx]] <- sample_features
}


# Combine all features into one data frame
panel_features_df_more_sigs <- as.data.frame(do.call(rbind, all_panel_features_more_sigs))

# Add a sample name column
panel_features_df_more_sigs$sample_names <- panel_mutation_data$Tumor_Sample_Barcode

# Add a panel column witht the name of each panel, in order to be able to separately test the model on each panel
panel_features_df_more_sigs <- panel_features_df_more_sigs %>%
  mutate(panel = case_when(
    # Group the samples from different version of the same panel under one panel name
    sample_names %in% c(panel1_samples, panel2_samples, panel2_samples) ~ "TARGET_SEQ_1",
    sample_names %in% c(targeted_panel_v1, targeted_panel_v2, tagreted_panel_v3) ~ "TARGET_SEQ_2",
    TRUE ~ "Other"
  ))

# Make separate dataframes for each panel
panel_features_df_more_sigs <- panel_features_df_more_sigs %>% filter(panel == "TARGET_SEQ_1")

targeted_panel_features_df_more_sigs <- panel_features_df_more_sigs %>% filter(panel == "TARGET_SEQ_2")
