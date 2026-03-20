# ERCC2mut
Composite mutational signature of ERCC2 deficiency (and more broadly Nucleotide Excision Repair (NER) deficiency) derived from whole exome sequencing (WES) or targeted panel sequencing data.

# Implementing ERCC2mut
To successfully implement the ERCC2mut classifier to another dataset, the following workflow is recommended

### Step 1
Loading the MAF-style dataframe and processing it following the instructions in [Processing_maf_file.txt](Processing_maf_file.txt). The MAF-style input dataframe should contain somatic mutation information for all samples where each row represents one mutation. An example input is shown in [ex_MAF_file.tsv](example_inputs/ex_MAF_file.tsv). A file with patient clinical information is necessary, because it is possible that there are multiple samples per single patient, so we recommend performing this filtering step in order to have only one sample per patient. An example patient clinical information file can is shown in [ex_clinical_information_file.tsv](example_inputs/ex_clinical_information_file.tsv).
      
### Step 2
Creating harmonized SNV, INS and DEL TMB values that are used as features for testing the ERCC2mut model. The script [Harmonized_SNV_INS_DEL_calcularion.R](scripts/Harmonized_SNV_INS_DEL_calculation.R) provides the code needed to perform this step. To sucessfully harmonize the SNV, INS and DEL values, the size of the panel would need to be calculated. This can be done by using a BED file which contains genomic information regarding the panel(s) in question (and the genes included in the panel). An example BED file is shown in [ex_bed_file.tsv](example_inputs/ex_bed_file.tsv). The script shows how to first calculate the panel size and then how to harmonize the features. The input (MAF file) necessary for the feature harmonization step is the filtered MAF dataframe from step 1. The output of this step (step 2) is a dataframe with all samples (from all panels; 1 row per sample) which contains information about the raw, normalized and hamronized TMB values. An example output from this step is shown in [ex_step2_output.tsv](example_output_from_scripts/ex_step2_output.tsv).

### Step 3
Creating mutational matrices (SBS and ID) and afterwards calculating the cosine similarity values (some of which are used as features in the ERCC2mut model). [Mutatational_matrices_and_cosine_similarity.R](scripts/Mutational_matrices_and_cosine_similarity.R) script shows how to properly make the mutational matrices and get the cosine similarity values. The input for this step should be the filtered MAF file from step 1, whereas the final outputs are 2 matrices (SBS and ID) containing cosine similarity values with all mutational signatures found previously in bladder cancer. 
     
### Step 4
Collecting all features for the ERCC2mut model into one dataframe. A function to correctly do this can be found in [Create_features_dataframe.R](scripts/Create_features_dataframe.R) script. The inputs for this script are the outputs from the scripts in steps 2 and 3 (specifically the dataframe containing harmonized values and cosine similarity matrix). An example output from this step is shown in [ex_features_df.tsv](example_inputs/ex_features_df.tsv) (NOTE: the script has one additional step to create mutliple dataframes, one per panel type, so that the model can be tested on every panel individually). 
     
### Step 5
Download the trained ERCC2mut classifier. This can be done by downloading the [ERCC2mut.rds](ERCC2mut.rds) file.
     
### Step 6
Load the ERCC2mut classifier and test on your own data. The script first selects only the 6 features necessary for running the ERCC2mut classifier (the output from step 5 can be used as an input in here). Further details on how to correctly load the classifier, get prediction scores, calculate cohort-specific cutpoint and plot the AUC ROC curve are provided in the [Test_ERCC2mut_model.R](scripts/Test_ERCC2mut_model.R) script.
