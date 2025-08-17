library(dplyr)
library(stats)
library(lme4)

# Set WD
setwd('/Users/lasyasreepada/Projects/CoMeT/')

# Read in the latest version of the longitudinal harmonized data
adni <- read.csv('data/ADNI/MRI_DNA_COG_XSH.csv')
adni <- subset(adni, select = -c(X.1,X))

# Define ROIs
rois_H <- grep("_H", colnames(adni), value = TRUE)

# Remove rows without harmonized imaging
adniMRI <- adni[complete.cases(adni[,rois_H]),]

# DATA WRANGLING

# Compute mean bilateral thickness measures
adniMRI <- adniMRI %>%
  rowwise() %>%
  mutate(M_Pericalcarine_Thk = mean(c(ST48TA_H,ST107TA_H)),
         M_Insula_Thk = mean(c(ST129TA_H,ST130TA_H)),
         M_Cuneus_Thk = mean(c(ST23TA_H,ST82TA_H)),
         M_Fusiform_Thk = mean(c(ST26TA_H,ST85TA_H)),
         M_LateralOccipital_Thk = mean(c(ST35TA_H,ST94TA_H)),
         M_Precentral_Thk = mean(c(ST51TA_H,ST110TA_H)),
         M_InferiorFrontal_Thk = mean(c(ST45TA_H,ST104TA_H,ST47TA_H,ST106TA_H)),
         M_InferiorTemporal_Thk = mean(c(ST32TA_H,ST91TA_H)),
         M_TemporalPole_Thk = mean(c(ST60TA_H,ST119TA_H)),
         M_SuperiorParietal_Thk = mean(c(ST57TA_H,ST116TA_H)),
         M_Precuneus_Thk = mean(c(ST52TA_H,ST111TA_H)),
         M_Supramarginal_Thk = mean(c(ST59TA_H,ST118TA_H)),
         M_SuperiorFrontal_Thk = mean(c(ST56TA_H,ST115TA_H)),
         M_MidFrontal_Thk = mean(c(ST15TA_H,ST74TA_H,ST55TA_H,ST114TA_H)), # Caudal LR, Rostral LR
         M_InferiorParietal_Thk = mean(c(ST31TA_H,ST90TA_H)),
         M_PosteriorCingulate_Thk = mean(c(ST50TA_H,ST109TA_H)),
         M_SuperiorTemporal_Thk = mean(c(ST58TA_H,ST117TA_H)),
         M_MiddleTemporal_Thk = mean(c(ST40TA_H,ST99TA_H)),
         M_MedialOrbitoFrontal_Thk = mean(c(ST39TA_H,ST98TA_H)),
         M_IsthmusCingulate_Thk = mean(c(ST34TA_H,ST93TA_H)),
         M_Lingual_Thk = mean(c(ST38TA_H,ST97TA_H)),
         M_Parahippocampal_Thk = mean(c(ST44TA_H,ST103TA_H)),
         M_Entorhinal_Thk = mean(c(ST24TA_H,ST83TA_H)),
         M_Hippocampus_Vol = mean(c(ST29SV_H,ST88SV_H)),
         M_Amygdala_Vol = mean(c(ST12SV_H,ST71SV_H)))

thk <- c('M_Pericalcarine_Thk','M_Insula_Thk','M_Cuneus_Thk','M_Fusiform_Thk',
         'M_LateralOccipital_Thk','M_Precentral_Thk','M_InferiorFrontal_Thk','M_InferiorTemporal_Thk',
         'M_TemporalPole_Thk','M_SuperiorParietal_Thk','M_Precuneus_Thk','M_Supramarginal_Thk',
         'M_SuperiorFrontal_Thk','M_MidFrontal_Thk','M_InferiorParietal_Thk','M_PosteriorCingulate_Thk',
         'M_SuperiorTemporal_Thk','M_MiddleTemporal_Thk','M_MedialOrbitoFrontal_Thk','M_IsthmusCingulate_Thk',
         'M_Lingual_Thk','M_Parahippocampal_Thk','M_Entorhinal_Thk')

vol <- c('M_Hippocampus_Vol','M_Amygdala_Vol')
  
# Loop for Thickness
for (roi in thk) {
  
  # Fit regression to controls
  formula <- as.formula(paste0(roi,'~ AGE + PTGENDER_M'))
  lmem <- lm(formula, data = subset(adniMRI,Case==0))
  
  # Predict for all subjects
  prediction <- predict(lmem,newdata = adniMRI)
  
  # Adjust values
  residual <- adniMRI[roi] - prediction
  # adj <- residual + coef(lm)[1]
  name_residual <- paste0(roi,'_res')
  adniMRI[name_residual] <- as.numeric(unlist(residual))
  
  # Standardize
  zscore <- (residual - mean(unlist(na.omit(adniMRI[adniMRI$Case==0,name_residual])))) / sd(unlist(na.omit(adniMRI[adniMRI$Case==0,name_residual])))
  name_z <- paste0(roi,'_wscore')
  adniMRI[name_z] <- as.numeric(unlist(zscore))
  
}

# Loop for Volume
for (roi in vol) {
  
  # Fit regression to controls
  formula <- as.formula(paste0(roi,'~ AGE + PTGENDER_M + ST10CV_H'))
  lmem <- lm(formula, data = subset(adniMRI,Case==0))
  
  # Predict for all subjects
  prediction <- predict(lmem,newdata = adniMRI)
  
  # Adjust values
  residual <- adniMRI[roi] - prediction
  # adj <- residual + coef(lm)[1]
  name_residual <- paste0(roi,'_res')
  adniMRI[name_residual] <- as.numeric(unlist(residual))
  
  # Standardize
  zscore <- (residual - mean(unlist(na.omit(adniMRI[adniMRI$Case==0,name_residual])))) / sd(unlist(na.omit(adniMRI[adniMRI$Case==0,name_residual])))
  name_z <- paste0(roi,'_wscore')
  adniMRI[name_z] <- as.numeric(unlist(zscore))
  
}

# Compute signatures
adniMRI <- adniMRI %>%
  rowwise() %>%
  mutate(Corticalv1 = mean(c(M_InferiorTemporal_Thk_wscore, 
                             M_TemporalPole_Thk_wscore, 
                             M_SuperiorParietal_Thk_wscore, 
                             M_Precuneus_Thk_wscore, 
                             M_Supramarginal_Thk_wscore, 
                             M_SuperiorFrontal_Thk_wscore,
                             M_MidFrontal_Thk_wscore)),
         Corticalv2 = mean(c(M_InferiorTemporal_Thk_wscore, 
                             M_TemporalPole_Thk_wscore, 
                             M_SuperiorParietal_Thk_wscore, 
                             M_Precuneus_Thk_wscore, 
                             M_Supramarginal_Thk_wscore, 
                             M_SuperiorFrontal_Thk_wscore,
                             M_MidFrontal_Thk_wscore,
                             M_Fusiform_Thk_wscore,
                             M_InferiorParietal_Thk_wscore,
                             M_SuperiorTemporal_Thk_wscore,
                             M_PosteriorCingulate_Thk_wscore,
                             M_MiddleTemporal_Thk_wscore)),
         Corticalv3 = mean(c(M_InferiorTemporal_Thk_wscore, 
                             M_TemporalPole_Thk_wscore, 
                             M_SuperiorParietal_Thk_wscore, 
                             M_Precuneus_Thk_wscore, 
                             M_Supramarginal_Thk_wscore, 
                             M_SuperiorFrontal_Thk_wscore,
                             M_MidFrontal_Thk_wscore,
                             M_Fusiform_Thk_wscore,
                             M_InferiorParietal_Thk_wscore,
                             M_SuperiorTemporal_Thk_wscore,
                             M_PosteriorCingulate_Thk_wscore,
                             M_MiddleTemporal_Thk_wscore,
                             M_MedialOrbitoFrontal_Thk_wscore)),
         Corticalv4 = mean(c(M_InferiorTemporal_Thk_wscore, 
                             M_TemporalPole_Thk_wscore, 
                             M_SuperiorParietal_Thk_wscore, 
                             M_Precuneus_Thk_wscore, 
                             M_Supramarginal_Thk_wscore, 
                             M_SuperiorFrontal_Thk_wscore,
                             M_MidFrontal_Thk_wscore,
                             M_Fusiform_Thk_wscore,
                             M_InferiorParietal_Thk_wscore,
                             M_SuperiorTemporal_Thk_wscore,
                             M_PosteriorCingulate_Thk_wscore,
                             M_IsthmusCingulate_Thk_wscore,
                             M_LateralOccipital_Thk_wscore,
                             M_MiddleTemporal_Thk_wscore,
                             M_MedialOrbitoFrontal_Thk_wscore,
                             M_Cuneus_Thk_wscore,
                             M_Lingual_Thk_wscore,
                             M_Pericalcarine_Thk_wscore)),
         MTLa = mean(c(M_Entorhinal_Thk_wscore,M_Parahippocampal_Thk_wscore)),
         MTLb = mean(c(M_Entorhinal_Thk_wscore,M_Parahippocampal_Thk_wscore,M_Hippocampus_Vol_wscore,M_Amygdala_Vol_wscore)))

# COMPUTE CoMeTs
# 1a
comet_v1a <- adniMRI['Corticalv1'] - adniMRI['MTLa']
name <- 'CoMeTv1a'
adniMRI[name] <- comet_v1a
# 1b
comet_v1b <- adniMRI['Corticalv1'] - adniMRI['MTLb']
name <- 'CoMeTv1b'
adniMRI[name] <- comet_v1b
# 2a
comet_v2a <- adniMRI['Corticalv2'] - adniMRI['MTLa']
name <- 'CoMeTv2a'
adniMRI[name] <- comet_v2a
# 2b
comet_v2b <- adniMRI['Corticalv2'] - adniMRI['MTLb']
name <- 'CoMeTv2b'
adniMRI[name] <- comet_v2b
# 3a
comet_v3a <- adniMRI['Corticalv3'] - adniMRI['MTLa']
name <- 'CoMeTv3a'
adniMRI[name] <- comet_v3a
# 3b
comet_v3b <- adniMRI['Corticalv3'] - adniMRI['MTLb']
name <- 'CoMeTv3b'
adniMRI[name] <- comet_v3b
# 4a
comet_v4a <- adniMRI['Corticalv4'] - adniMRI['MTLa']
name <- 'CoMeTv4a'
adniMRI[name] <- comet_v4a
# 4b
comet_v4b <- adniMRI['Corticalv4'] - adniMRI['MTLb']
name <- 'CoMeTv4b'
adniMRI[name] <- comet_v4b

# COGNITIVE
adniMRI <- adniMRI %>%
  rowwise() %>%
  mutate(MIN = min(c(PHC_EXF_nearest, PHC_LAN_nearest, PHC_VSP_nearest), na.rm = TRUE))

cog <- c('PHC_EXF_nearest','PHC_LAN_nearest','PHC_VSP_nearest','PHC_MEM_nearest','AVRECTOT_nearest','MIN')

# Loop for cognition
for (roi in cog) {
  
  # Fit regression to controls
  formula <- as.formula(paste0(roi,'~ AGE + PTGENDER_M'))
  lmem <- lm(formula, data = subset(adniMRI,Case==0))
  
  # Predict for all subjects
  prediction <- predict(lmem,newdata = adniMRI)
  
  # Adjust values
  residual <- adniMRI[roi] - prediction
  # adj <- residual + coef(lm)[1]
  name_residual <- paste0(roi,'_res')
  adniMRI[name_residual] <- as.numeric(unlist(residual))
  
  # Standardize
  zscore <- (residual - mean(unlist(na.omit(adniMRI[adniMRI$Case==0,name_residual])))) / sd(unlist(na.omit(adniMRI[adniMRI$Case==0,name_residual])))
  name_z <- paste0(roi,'_wscore')
  adniMRI[name_z] <- as.numeric(unlist(zscore))
  
}

# COMPUTE CoMeTs

# EXF
comet_exf <- adniMRI['PHC_EXF_nearest_wscore'] - adniMRI['AVRECTOT_nearest_wscore']
name <- 'CoMeT_EXF'
adniMRI[name] <- comet_exf

# LAN
comet_lan <- adniMRI['PHC_LAN_nearest_wscore'] - adniMRI['AVRECTOT_nearest_wscore']
name <- 'CoMeT_LAN'
adniMRI[name] <- comet_lan

# VSP
comet_vsp <- adniMRI['PHC_VSP_nearest_wscore'] - adniMRI['AVRECTOT_nearest_wscore']
name <- 'CoMeT_VSP'
adniMRI[name] <- comet_vsp

# MIN
comet_min <- adniMRI['MIN_wscore'] - adniMRI['AVRECTOT_nearest_wscore']
name <- 'CoMeT_MIN'
adniMRI[name] <- comet_min

# Save file 
write.csv(adniMRI,'data/ADNI/MRI_DNA_COG_XSH_COMET.csv')

# adniMRI <- adniMRI %>%
#   rowwise() %>%
#   mutate(L_CorticalLOAD = sum(ST32TA_H, 
#                               ST60TA_H, 
#                               ST57TA_H, 
#                               ST52TA_H, 
#                               ST59TA_H, 
#                               ST56TA_H,
#                               ST15TA_H,ST55TA_H),
#          R_CorticalLOAD = sum(ST91TA_H, 
#                               ST119TA_H, 
#                               ST116TA_H, 
#                               ST111TA_H, 
#                               ST118TA_H, 
#                               ST117TA_H,
#                               ST74TA_H,ST114TA_H),
#          CorticalLOAD_AI = abs(L_CorticalLOAD - R_CorticalLOAD) / (L_CorticalLOAD + R_CorticalLOAD))
