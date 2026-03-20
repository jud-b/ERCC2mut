# This code shows how to test the panel-based ERCC2 model in your own dataset. 
# The model can be tested with or without knowing the ERCC2 mutation status of each sample (MUT, WT, WT-like, etc), 
# however in order to be able to evaluate the performace in your own dataset and to find cohort-specific cutpoint,
# the model needs to be tested only on the MUT and WT classes and the true labels need to be known.
# Below you can find the implpementation for testing and evaluating the performance, if not needed the Cohort-spceific cutpoint and Evaluate the model in own dataset steps can be skipped

# Load necessary libraries
library(caret)
library(gbm3)
library(gbm)
library(cutpointr)
library(pROC)

###################################################################### Testing and evaluating the model in own dataset ##################################################################
# Assuming the ERCC2 status of the samples is known, subset only for the true ERCC2-MUT and ERCC2-WT ssamples
ERCC2_mut_wt_panel_more_sigs <- panel_features_df_more_sigs %>% filter(label %in% c("ERCC2-MUT", "ERCC2-WT")) 

# Select the 6 features for the model
ERCC2_mut_wt_panel_more_sigs_6f <- ERCC2_mut_wt_panel_more_sigs %>% dplyr::select(cosine_sbs5, cosine_sbs13, cosine_id8, SNV_TMB_z, cosine_sbs2, DEL_TMB_z,label)

# Add a binary label for the two classes (MUT=1, WT=0) and input 0's for all NA values
ERCC2_mut_wt_panel_more_sigs_6f$label <- ifelse(ERCC2_mut_wt_panel_more_sigs_6f$label == "ERCC2-MUT", 1, 0)
ERCC2_mut_wt_panel_more_sigs_6f$label <- as.factor(ERCC2_mut_wt_panel_more_sigs_6f$label)
ERCC2_mut_wt_panel_more_sigs_6f[is.na(ERCC2_mut_wt_panel_more_sigs_6f)] <- 0

# To avoid retraining the model, download the already trained model "ERCC2mut.rds" that can be found in the repository

# Load pretrained model
model_obj <- readRDS("ERCC2mut.rds")

gbm_model <- model_obj$model
best_params <- model_obj$best_params
feature_names <- model_obj$feature_names

set.seed(123)
ERCC2mut_scores_panel <- predict(
  gbm_model, 
  newdata = ERCC2_mut_wt_panel_more_sigs_6f %>% dplyr::select(-label), 
  n.trees = best_params$n.trees, 
  type = "response"
)

################################################################################ Cohort-spceific cutpoint ##############################################################################
# This step is optional and can be replaced with the cutpoint that we have established from the training set (0.04821027)

# Make a dataframe with sample labels and the scroes from ERCC2mut model
scores_labels_panel <- data.frame(scores = ERCC2mut_scores_panel, labels = ERCC2_mut_wt_panel_more_sigs_6f$label)

set.seed(123)
cutpoint_result_panel <- cutpointr(scores_labels_panel, x = scores, class = labels, pos_class = 1)

# Extract best threshold
set.seed(123)
optimal_cutpoint_panel <- cutpoint_result_panel$optimal_cutpoint

################################################################################ Classify the samples as MUT or WT based ##############################################################################
# Apply the cutpoit to classify the samples as either ERCC2-WT or ERCC2-MUT
# Samples with scores >= cutpoint are classified as MUT, otherwise as WT

set.seed(123)
predicted_class_panel <- ifelse(ERCC2mut_scores_panel >= optimal_cutpoint_panel, 1, 0)

################################################################################ Evaluate the model in own dataset  ##############################################################################
# Generate a confusion matrix
evaluate_ERCC2mut_panel <- confusionMatrix(as.factor(predicted_class_panel), as.factor(ERCC2_mut_wt_panel_more_sigs_6f$label), positive = "1") 

# Generate AUC ROC plot
roc_obj <- pROC::roc(
  response = ERCC2_mut_wt_panel_more_sigs_6f$label,  # true labels (binary)
  predictor = ERCC2mut_scores_panel,                 # predicted scores from the model
  levels = c("0", "1"),                              # specify order: WT, MUT
  direction = "<"                                    # "<" means higher score = MUT
)

panel_roc_df <- data.frame(sensitivity = roc_obj$sensitivities, specificity = roc_obj$specificities)
panel_roc_df$fpr <- 1-panel_roc_df$specificity

# Calculate AUC
panel_auc <- AUC(x = panel_roc_df$fpr, y = panel_roc_df$sensitivity, method = "trapezoid")

# Make the ROC AUC plot
ggplot(panel_roc_df, aes(x = fpr, y = sensitivity)) +
  geom_line(color = "red", size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  xlab("1 - specificity") +
  ylab("sensitivity") +
  theme_classic()+
  geom_label(
    aes(x = 0.9, y = 0.15, label = paste0("AUC = ", round(panel_auc, 2))),
    fill = "white",
    color = "black",
    label.size = 0.5
  )
