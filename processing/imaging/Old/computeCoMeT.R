library(dplyr)
library(stats)
library(lme4)

# Set WD
setwd('/Users/lasyasreepada/Projects/CoMeT/')

# Read in the latest version of the longitudinal harmonized data
adni <- read.csv('data/ADNI/MRI_DNA_COG_LONG_H.csv')
adni <- subset(adni, select = -c(X.1,X))

# Change type of some variables
adni$DX_nearest_1.0_ordinal <- as.factor(adni$DX_nearest_1.0_ordinal) # DX
adni$PTGENDER_ordinal <- as.factor(adni$PTGENDER_ordinal) # Gender

# Define ROIs
rois.combat <- grep('^ST\\d+TA.combat$', names(adni), value=TRUE) # All harmonized FreeSurfer ROIs

# Remove rows without harmonized imaging
adniMRI <- adni[complete.cases(adni[,rois.combat]),]

# DATA WRANGLING

# Define cognitive composite 
adniMRI <- adniMRI %>% 
  rowwise() %>% 
  mutate(M_Cortical_Cog = mean(c(PHC_EXF_nearest,PHC_VSP_nearest,PHC_LAN_nearest)))

adniMRI <- adniMRI %>% 
  group_by(RID) %>%
  mutate(Visit = n()) %>%
  ungroup()
 
# Compute mean bilateral thickness measures
adniMRI <- adniMRI %>%
  rowwise() %>%
  mutate(M_Pericalcarine_Thk = mean(c(ST48TA.combat,ST107TA.combat)),
         M_Insula_Thk = mean(c(ST129TA.combat,ST130TA.combat)),
         M_Cuneus_Thk = mean(c(ST23TA.combat,ST82TA.combat)),
         M_Fusiform_Thk = mean(c(ST26TA.combat,ST85TA.combat)),
         M_LateralOccipital_Thk = mean(c(ST35TA.combat,ST94TA.combat)),
         M_Precentral_Thk = mean(c(ST51TA.combat,ST110TA.combat)),
         M_InferiorFrontal_Thk = mean(c(ST45TA.combat,ST104TA.combat,ST47TA.combat,ST106TA.combat)),
         M_InferiorTemporal_Thk = mean(c(ST32TA.combat,ST91TA.combat)),
         M_TemporalPole_Thk = mean(c(ST60TA.combat,ST119TA.combat)),
         M_SuperiorParietal_Thk = mean(c(ST57TA.combat,ST116TA.combat)),
         M_Precuneus_Thk = mean(c(ST52TA.combat,ST111TA.combat)),
         M_Supramarginal_Thk = mean(c(ST59TA.combat,ST118TA.combat)),
         M_SuperiorFrontal_Thk = mean(c(ST56TA.combat,ST115TA.combat)),
         M_MidFrontal_Thk = mean(c(ST15TA.combat,ST74TA.combat,ST55TA.combat,ST114TA.combat)), # Caudal LR, Rostral LR
         M_InferiorParietal_Thk = mean(c(ST31TA.combat,ST90TA.combat)),
         M_PosteriorCingulate_Thk = mean(c(ST50TA.combat,ST109TA.combat)),
         M_SuperiorTemporal_Thk = mean(c(ST58TA.combat,ST117TA.combat)),
         M_Parahippocampal_Thk = mean(c(ST44TA.combat,ST103TA.combat)),
         M_Entorhinal_Thk = mean(c(ST24TA.combat,ST83TA.combat)),
         M_SuperiorTemporal_Thk = mean(c(ST58TA.combat,ST117TA.combat)),
         M_MiddleTemporal_Thk = mean(c(ST40TA.combat,ST99TA.combat)))

# Compute signatures
adniMRI <- adniMRI %>%
  rowwise() %>%
  mutate(CorticalLOAD = mean(c(M_InferiorTemporal_Thk, 
                               M_TemporalPole_Thk, 
                               M_SuperiorParietal_Thk, 
                               M_Precuneus_Thk, 
                               M_Supramarginal_Thk, 
                               M_SuperiorFrontal_Thk,
                               M_MidFrontal_Thk)),
         CorticalEOAD = mean(c(M_Fusiform_Thk,
                               M_SuperiorParietal_Thk,
                               M_Precuneus_Thk,
                               M_SuperiorFrontal_Thk,
                               M_MidFrontal_Thk,
                               M_InferiorParietal_Thk,
                               M_SuperiorTemporal_Thk,
                               M_PosteriorCingulate_Thk,
                               M_MiddleTemporal_Thk)),
         MTLComposite = mean(c(M_Entorhinal_Thk,
                               M_Parahippocampal_Thk)))

adniMRI <- adniMRI %>%
  rowwise() %>%
  mutate(L_CorticalLOAD = sum(ST32TA.combat, 
                              ST60TA.combat, 
                              ST57TA.combat, 
                              ST52TA.combat, 
                              ST59TA.combat, 
                              ST56TA.combat,
                              ST15TA.combat,ST55TA.combat),
         R_CorticalLOAD = sum(ST91TA.combat, 
                              ST119TA.combat, 
                              ST116TA.combat, 
                              ST111TA.combat, 
                              ST118TA.combat, 
                              ST117TA.combat,
                              ST74TA.combat,ST114TA.combat),
         CorticalLOAD_AI = abs(L_CorticalLOAD - R_CorticalLOAD) / (L_CorticalLOAD + R_CorticalLOAD))

adniMRI <- adniMRI %>%
  rowwise() %>%
  mutate(L_CorticalEOAD = sum(ST26TA.combat,
                              ST57TA.combat,
                              ST52TA.combat,
                              ST56TA.combat,
                              ST15TA.combat,ST55TA.combat,
                              ST31TA.combat,
                              ST58TA.combat,
                              ST50TA.combat,
                              ST40TA.combat),
         R_CorticalEOAD = sum(ST85TA.combat,
                              ST116TA.combat,
                              ST111TA.combat,
                              ST115TA.combat,
                              ST74TA.combat,ST114TA.combat,
                              ST90TA.combat,
                              ST117TA.combat,
                              ST109TA.combat,
                              ST99TA.combat),
         CorticalEOAD_AI = abs(L_CorticalEOAD - R_CorticalEOAD) / (L_CorticalEOAD + R_CorticalEOAD))

# CROSS SECTIONAL DATA EXTRACTION
adniDNAmX <- adniMRI %>%
  group_by(RID) %>%
  dplyr::slice(which.min(Years_to_DNAm)) %>% # For subjects with DNAm, select row with min Years_to_DNAm
  ungroup()

adniMRIX <- adniMRI %>%
  group_by(RID) %>%
  filter(is.na(Years_to_DNAm)) %>% # Subjects with no DNAm, only MRI
  dplyr::slice(which.min(AGE)) %>% # For subjects without DNAm, select row with min AGE
  ungroup()

adniX <- rbind(adniDNAmX,adniMRIX) # Concatenate the dataframes by rows

# ADJUSTED Z-SCORES
# Measures for z-scoring
rois.bilateral <- grep("^M_.+_Thk$", names(adniMRI), value=TRUE)
signatures <- c('CorticalLOAD','CorticalEOAD','MTLComposite','CorticalLOAD_AI','CorticalEOAD_AI')
cog <- c('M_Cortical_Cog','PHC_EXF_nearest','PHC_LAN_nearest','PHC_VSP_nearest','PHC_MEM_nearest','AVRECTOT','AVDELERR2','AVDELTOT')
 
# Loop
for (roi in c(rois.bilateral,signatures)) {

  # Fit regression to controls
  formula <- as.formula(paste0(roi,'~ AGE + PTGENDER_ordinal'))
  lmem <- lm(formula, data = subset(adniX,Case==0))

  # Predict for all subjects
  prediction <- predict(lmem,newdata = adniX)

  # Adjust values
  residual <- adniX[roi] - prediction
  # adj <- residual + coef(lm)[1]
  name_residual <- paste0(roi,'_Residual')
  adniX[name_residual] <- as.numeric(unlist(residual))

  # Standardize
  zscore <- (residual - mean(unlist(na.omit(adniX[adniX$Case==0,name_residual])))) / sd(unlist(na.omit(adniX[adniX$Case==0,name_residual])))
  name_z <- paste0(roi,'_Zscore')
  adniX[name_z] <- as.numeric(unlist(zscore))

}

for (roi in cog) {
  
  # Fit regression to controls
  formula <- as.formula(paste0(roi,'~ AGE + PTGENDER_ordinal + PTEDUCAT'))
  lmem <- lm(formula, data = subset(adniX,Case==0))
  
  # Predict for all subjects
  prediction <- predict(lmem,newdata = adniX)
  
  # Adjust values
  residual <- adniX[roi] - prediction
  # adj <- residual + coef(lm)[1]
  name_residual <- paste0(roi,'_Residual')
  adniX[name_residual] <- as.numeric(unlist(residual))
  
  # Standardize
  zscore <- (residual - mean(unlist(na.omit(adniX[adniX$Case==0,name_residual])))) / sd(unlist(na.omit(adniX[adniX$Case==0,name_residual])))
  name_z <- paste0(roi,'_Zscore')
  adniX[name_z] <- as.numeric(unlist(zscore))
  
}

# COMPUTE TYPICALITY INDEX
rois.z <- grep("^M_.+_Thk_Zscore$", names(adniX), value=TRUE)
signatures.z <- grep("OAD_Zscore$", names(adniX), value=TRUE)
cog.z <- c('M_Cortical_Cog_Zscore','PHC_EXF_nearest_Zscore','PHC_LAN_nearest_Zscore','PHC_VSP_nearest_Zscore')

# Imaging TI
for (roi in c(rois.z,signatures.z)) {
  ti <- adniX[roi] - adniX['MTLComposite_Zscore']
  name_ti <- paste0(gsub('_Zscore','',roi),'_TI')
  adniX[name_ti] <- ti
}

# Cognitive TI
for (score in cog.z) {
  ti <- adniX[score] - adniX['PHC_MEM_nearest_Zscore']
  name_ti <- paste0(gsub('_Zscore','',score),'_TI')
  adniX[name_ti] <- ti
}

# Save file 
write.csv(adniX,'data/ADNI/MRI_DNA_COG_X.csv')
