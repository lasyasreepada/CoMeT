library(dplyr)
library(stats)

# Set WD
setwd('/Users/lasyasreepada/Projects/CoMeT/')

# Read in the latest version of the longitudinal data
adni <- read.csv('data/ADNI/MRI_DNA_COG_LONG.csv')

# Define cases and controls based on DX and amyloid status
adni <- adni %>% mutate(Case = case_when(
  (AMYLOID_SUMMARY==0 & DX_nearest_1.0=='CN') ~ 0, # Control
  (AMYLOID_SUMMARY==1 & (DX_nearest_1.0=='Dementia' | DX_nearest_1.0=='MCI')) ~ 1, # Case
)) %>% filter(if_all(c(Case), ~ !is.na(.x)))

# Select features
names <- colnames(adni)
rois <- grep('^ST\\d+TA$', names(adni), value=TRUE) # All FreeSurfer ROIs
rois.missing <- c("ST123TA", "ST22TA", "ST64TA", "ST81TA") # Exclude ROIs with large amounts of missing data
rois <- rois[!(rois %in% rois.missing)]
volumes <- c('ST29SV','ST88SV','ST12SV','ST71SV','ST10CV')
rois <- c(rois,volumes)

# Complete Cases
features <- c('RID','AGE','Case','DX_nearest_1.0','PTGENDER','FLDSTRENG',rois)
adniNoNA <- adni[complete.cases(adni[ ,features]),] 

# Select cross sectional data

# DNA only
adniDNAmX <- adniNoNA %>%
  group_by(RID) %>%
  dplyr::slice(which.min(Years_to_DNAm)) %>% # For subjects with DNAm, select row with min Years_to_DNAm
  ungroup()

# MRI only
adniMRIX <- adniNoNA %>%
  group_by(RID) %>%
  filter(is.na(Years_to_DNAm)) %>% # Subjects with no DNAm, only MRI
  dplyr::slice(which.min(AGE)) %>% # For subjects without DNAm, select row with min AGE
  ungroup()

# Combine
adniX <- rbind(adniDNAmX,adniMRIX) # Concatenate the dataframes by rows

# WRITE
write.csv(adniX,'data/ADNI/MRI_DNA_COG_XS.csv',row.names = FALSE)

