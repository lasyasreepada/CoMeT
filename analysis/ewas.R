library(dplyr)
library(stringr)
library(sesame)
library(lubridate)
library(SummarizedExperiment)
library(tibble)
library(devtools)
library(RSQLite)
library(maxprobes)
library(neuroCombat)

sesameDataCache()

# READ DATA

# Set WD
setwd('/Users/lasyasreepada/Projects/CoMeT/')

# Read in betas matrix
se <- readRDS("data/ADNI/methylation/methylation.rds")
comet <- read.csv('data/ADNI/MRI_DNA_COG_XSH_COMET_Grouped.csv')

# Get phenotype data
coldata <- as.data.frame(colData(se))
rownames(coldata) <- NULL

comet <- comet %>%
  select(-c(Sample,PlateNumber))

coldata <- coldata %>%
  select(c(RID, Sample, Age, PlateNumber))

# Merge
adni <- comet %>%
  left_join(coldata, by=c('RID'), relationship = 'many-to-many') %>%
  mutate(Diff = abs(AGE_SAMPLE_nearest - Age)) %>%
  group_by(RID) %>%
  slice_min(Diff, with_ties = FALSE) %>%
  ungroup() %>%
  select(-c(Diff, X, X.1))

# PREPARE DATASET

# Remove missing features
adni <- adni %>%
  filter(if_all(c(Age, Sample, CoMeTv4b), ~ !is.na(.x)))

# Final Data
adni <- adni[adni$Case==1,]

# SE 

# Select samples
samples <- adni$Sample
se <- se[,samples]

# CELL TYPE DECONVOLUTION
# install_github("zhou-lab/knowYourCG")
betas <- assays(se)[[1]]

# Cell type proportions
library(EpiDISH)
library(FlowSorted.Blood.EPIC)
data(cent12CT.m)

# Set reference matrix
ref <- cent12CT.m

# Intersect
common_cpgs <- intersect(rownames(betas), rownames(ref))
betas <- betas[common_cpgs,]
ref <- ref[common_cpgs,]

# Run deconvolution
deconv <- epidish(beta.m = betas, 
                  ref.m = ref, 
                  method = "RPC")

# Extract cell type proportions
cell_proportions <- as.data.frame(deconv$estF)
cell_proportions <- rownames_to_column(cell_proportions,var = "Sample")

# Merge
adni <- left_join(adni,cell_proportions,by="Sample")

# PREPARE FINAL SE
betas <- assays(se)[[1]]
coldata <- as.data.frame(adni)
rownames(coldata) <- coldata$Sample

# Create SE
se <- SummarizedExperiment(assays = betas, colData=coldata)

# Limit to CpGs in common between EPIC v1 and EPIC v2
common_cpgs <- read.csv('/Users/lasyasreepada/Library/CloudStorage/Box-Box/INDD_Methylation_Processed/cpgs_v1_v2_shared.csv')
common_cpgs <- common_cpgs$CPG
se <- se[common_cpgs,]

# Filter cross-reactive probes
library(maxprobes)
xloci <- maxprobes::xreactive_probes(array_type = "EPIC")

cg_ok <- !rownames(se) %in% xloci
se <- se[cg_ok,]

# MISSING VALUE HANDLING

# Function to remove rows and features with too much missingness
cleanMatrixForClusterSE <- function(se, f_row = 0, f_col = 0.4) {
  mtx = assays(se)[[1]]
  cat(sprintf("Filter rows with >%1.2f missingness and columns with >%1.2f missingness.\n",
              f_row, f_col))
  cat("Before: ", nrow(mtx), "rows and ", ncol(mtx),"columns.\n")
  namtx = is.na(mtx)
  good_row = rowSums(namtx) <= ncol(mtx) * f_row
  good_col = colSums(namtx) <= nrow(mtx) * f_col
  cat("After: ", sum(good_row), "rows and ", sum(good_col),"columns.\n")
  return(se[good_row, good_col])
}

# Filter
se <- cleanMatrixForClusterSE(se)

# Categorical level checks
colData(se)$Sex <- relevel(factor(colData(se)$Sex), "Female")
cg_ok <- (checkLevels(assay(se), colData(se)$Sex))

# Filter
se <- se[cg_ok,]

# SAVE
saveRDS(se,'~/Projects/CoMeT/data/ADNI/20250504_COMET_SE.rds')

##############################################################

# READ DATA
se <- readRDS('~/Projects/CoMeT/data/ADNI/20250504_COMET_SE.rds')

# RUN EWAS

# CoMeT
res <- DML(se,~ AGE_SAMPLE_nearest+PTGENDER_M+CoMeTv4b, BPPARAM=BiocParallel::MulticoreParam(6))

# Test
smry <- summaryExtractTest(res)
smry$FDR <- p.adjust(smry$Pval_CoMeTv4b, method="fdr")
smry %>% arrange(Pval_CoMeTv4b) %>% dplyr::select(Probe_ID,Est_CoMeTv4b,Pval_CoMeTv4b,FDR) %>%
  filter(Pval_CoMeTv4b < 1e-5)

# SAVE 
saveRDS(smry, '~/Projects/CoMeT/data/ADNI/20250504_COMET_SMRY_NOAGE.rds')
