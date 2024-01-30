# Calculation of Horvath DNAm Clock, DNAmCorticalClock and neuronal cell type 
# proportions from EPIC array script for the paper: "Cell-type-specific 
# ageing effects in the human OFC and implications for psychiatric disease" by 
# Fröhlich et al.
# ==============================================================================

# Author: Natan Yusupov, natan_yusupov@psych.mpg.de
# Supervision: Dr. Darina Czamara

# The Horvath clock was calculated using the methylclock package as described by 
# Pelegí-Sisó et al. Bioinformatics 2021 (DOI: 10.1093/bioinformatics/btaa825). 
# See also details in https://github.com/isglobal-brge/methylclock/blob/main/README.md
# The CorticalClock was calculated according to Shireby et al., Brain 2020 
# (DOI:10.1093/brain/awaa334). See also details in https://github.com/gemmashireby/CorticalClock
# Neuronal cell type proportions were calculated from DNAm data according to 
# Guintivano et al. Epigenetcis 2013 (DOI: 10.4161/epi.23924)

################################################################################

# Load libraries
library(tidyverse)
library(methylclock)
library(impute)
library(Rcpp)
library(tibble)
library(ggpubr)
library(readr)
library(cowplot)
library(Biobase)
library(ggpmisc)
library(GEOquery)
library(Hmisc)
library(htmlTable)
library(broom)
library(DataExplorer)
library(corrplot)
library(AnnotationDbi)
library(GeneAnswers)
library(performance)
library(FlowSorted.DLPFC.450k)
library(sjPlot)
getwd()

# Load data
load("BMIQ.quantileN_combated.Rdata")
load("betas_normalized_filtered_combated.Rdata")
load("pd_qual.Rdata")
load("RGSet_cleaned_final.Rdata")

dataset_norm1 <- BMIQ.quantileN_combated # for Horvath clock
dataset_norm2 <- betas_normalized_filtered_combated # for Cortical clock
pheno <- pd_qual %>%
  rownames_to_column(var = "Sample_ID")

# Calculate DNAm age estimation and acceleration
# Comment: since freezer storage and brain pH showed batch effects on DNAm PCs 
# in the preprocessing procedure, we calculate DNAm age acceleration (AgeAccel)
# by regressing these variables along age as well

# Calculate Horvath DNAm age estimation and DNAm age acceleration
Horvath_age_norm1 <- DNAmAge(dataset_norm1, clocks = c("Horvath"),
                             age = pheno$Age, normalize = FALSE, cell.count = FALSE)
Final_df <- Horvath_age_norm1 %>%
  dplyr::select(Sample_ID = "id", DNAmAgeHorvath_norm1 = "Horvath") %>%
  left_join(pheno, by = "Sample_ID")

regressors <- c("Age", "Freezer_storage", "Brain_ph")
f_Horvath_norm1 <- as.formula(paste0("DNAmAgeHorvath_norm1", "~", paste(regressors, collapse = "+")))

Final_df <- Final_df %>%
  mutate(AgeAccelHorvath_norm1_corrected = augment(lm(formula = f_Horvath_norm1, 
                                                      data = Final_df)) %>% pull(.resid))

# Calculate DNAmClockCortical and DNAm age acceleration
source("addFiles/PredCorticalAge/CorticalClock.r") # source 'CorticalClock.r' script
# Comment: WE use the Github code of CorticalClock function from 
# https://github.com/gemmashireby/CorticalClock with minimal changes: remove autoplot
# and include two new arguments: 'output' (location to save) and 'suffix' (add suffix to the name)
CorticalClock(betas = dataset_norm2, pheno = pheno, dir = "addFiles/PredCorticalAge/", 
              IDcol = "Sample_ID", Agecol = "Age", output = "addFiles/PredCorticalAge/", 
              suffix = "norm2")
# Output:
# "some probes are missing, we will need to impute values here - the final predictions will be less accurate"
#                 Cortical Clock
# Correlation (r)           0.97
# RMSE (years)              6.05
# MAD (years)               4.56
# How many & which CpGs are missing in Norm2 for CorticalClock
CorticalClock_coef <- read.delim("addFiles/PredCorticalAge/CorticalClockCoefs.txt", sep = "")
length(setdiff(CorticalClock$probe, rownames(dataset_norm2)))
setdiff(CorticalClock$probe, rownames(dataset_norm2))

CorticalPred_norm2 <- read_delim("addFiles/CorticalPred_norm2.csv", delim = ",")

Final_df <- Final_df %>%
   left_join(CorticalPred_norm2 %>% dplyr::select(Sample_ID = "ID", DNAmAgeCorticalClock_norm2 = "brainpred"),
            by = "Sample_ID")

f_CorticalClock_norm2 <- as.formula(paste0("DNAmAgeCorticalClock_norm2", "~", paste(regressors, collapse = "+")))

Final_df <- Final_df %>%
  mutate(AgeAccelCorticalClock_norm2_corrected = augment(lm(formula = f_CorticalClock_norm2, 
                                                            data = Final_df)) %>% pull(.resid))

# Calculate neuronal cell type proportions
cell_counts <- minfi::estimateCellCounts(RGSet_cleaned_final,
                                         compositeCellType = "DLPFC",
                                         probeSelect = "auto",
                                         processMethod = "preprocessQuantile",
                                         cellTypes = c("NeuN_pos","NeuN_neg"),
                                         referencePlatform = "IlluminaHumanMethylation450k")
cell_proportions_df <- as.data.frame(cell_counts)
cell_proportions_df <- cell_proportions_df %>%
  left_join(pheno %>% dplyr::select(Sample_ID, Sample_Name, propNeuron=NeuN_pos))

Final_df <- Final_df %>%
  left_join(cell_proportions_df %>% dplyr::select(Sample_ID, propNeuron), by = "Sample_ID")

Final_df_no_duplicates_selected <- Final_df %>% filter(Technical_duplicate == "no") %>%
  dplyr::select(Sample_ID, Sample_Name, DNAmAgeHorvath_norm1, 
                AgeAccelHorvath_norm1_corrected, DNAmAgeCorticalClock_norm2, 
                AgeAccelCorticalClock_norm2_corrected, propNeuron, Sex = "Gender", Age,
                Smoking_status, Disease_status)
write.xlsx(Final_df_no_duplicates_selected, "Final_df_no_duplicates_selected.xlsx")

# End of Script

