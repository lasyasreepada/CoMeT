library(longCombat)
library(invgamma)
library(lme4)
library(dplyr)

# Set WD
setwd('/Users/lasyasreepada/Projects/CoMeT/')

# Data IO
adniLongFile <- "data/ADNI/MRI_DNA_COG_LONG.csv"
adniLong <- read.csv(adniLongFile)

# Define time variable for harmonization
adniLong$Time <- adniLong$AGE - adniLong$AGE_BL

# Select features
names <- colnames(adniLong)
rois <- grep('^ST\\d+TA$', names(adniLong), value=TRUE) # All FreeSurfer ROIs
rois.missing <- c("ST123TA", "ST22TA", "ST64TA", "ST81TA") # Exclude ROIs with large amounts of missing data
rois <- rois[!(rois %in% rois.missing)]

# Encode categorical variables
encode_ordinal <- function(x, order = unique(x)) {
  x <- as.numeric(factor(x, levels = order, exclude = NULL))
  x
}

adniLong$DX_nearest_1.0_ordinal <- encode_ordinal(adniLong[["DX_nearest_1.0"]], order = c("CN", "MCI", "Dementia"))
adniLong$PTGENDER_ordinal <- encode_ordinal(adniLong[["PTGENDER"]], order = c("Male", "Female"))

# Define cases and controls based on DX and amyloid status
adniLong <- adniLong %>% mutate(Case = case_when(
  (AMYLOID_SUMMARY==0 & DX_nearest_1.0_ordinal==1) ~ 0, # Control
  (AMYLOID_SUMMARY==1 & DX_nearest_1.0_ordinal!=1) ~ 1, # Case
))

adniLong$Case_ordinal <- encode_ordinal(adniLong[["Case"]], order = c("0", "1"))

featureNames <- c('RID','AGE','Case_ordinal','DX_BL','DX_nearest_1.0_ordinal','Time','PTGENDER_ordinal','FLDSTRENG',rois)
adniLongNoNA <- adniLong[complete.cases(adniLong[ ,featureNames]),] 

# Plots

# Visualize change in batch over time
batchTimeViz(batchvar='FLDSTRENG',timevar='Time',data=adniLongNoNA)

# Test for additive scanner effects
addTestTable <- addTest(idvar="RID",
                        batchvar = 'FLDSTRENG',
                        features = rois,
                        formula = 'AGE + PTGENDER_ordinal + DX_BL + DX_nearest_1.0_ordinal*Time',
                        ranef = '(Time|RID)',
                        data = adniLongNoNA)

# Test for multiplicative scanner effects
multTestTable <- multTest(idvar="RID",
                          batchvar = 'FLDSTRENG',
                          features = rois,
                          formula = 'AGE + PTGENDER_ordinal + DX_BL + DX_nearest_1.0_ordinal*Time',
                          ranef = '(Time|RID)',
                          data = adniLongNoNA)

# # Another example - highest multiplicative effect var
# batchBoxplot(idvar='RID',batchvar='FLDSTRENG',feature = 'ST51TA',
#              formula='AGE + PTGENDER_ordinal + DX_BL + DX_nearest_1.0_ordinal*Time',ranef='(1 + Time|RID)',
#              data=adniLongNoNA)
# 
# batchBoxplot(idvar='RID',batchvar='FLDSTRENG',feature = 'ST51TA',
#              formula='AGE + PTGENDER_ordinal + DX_BL + DX_nearest_1.0_ordinal*Time',ranef='(1 + Time|RID)',
#              data=adniLongNoNA, adjustBatch = TRUE, orderby='var')

# Visualize trajectories
# trajPlot(idvar='RID',timevar='Time',batchvar='FLDSTRENG',feature = 'ST51TA',
#          data=adniLongNoNA, 
#          line.col=adniLongNoNA$DX_nearest_1.0_ordinal[!duplicated(adniLongNoNA$RID)]+1) # everyone

# Fit combat
adniLongCombat <- longCombat(idvar='RID',
                             timevar = 'Time',
                             batchvar = 'FLDSTRENG',
                             features = rois,
                             formula = 'AGE + PTGENDER_ordinal + DX_nearest_1.0_ordinal*Time',
                             ranef='(Time|RID)',
                             data = subset(adniLongNoNA,select=featureNames))

# Get the harmonized data
rois.combat <- paste(rois, ".combat", sep="")
adniLongHarmonized <- adniLongCombat$data_combat

par(mfrow=c(2,2))
# Visualize the change before and after adjusting for batch - highest additive effect var
batchBoxplot(idvar='RID',batchvar='FLDSTRENG',feature = 'ST111TA',
             formula='AGE + PTGENDER_ordinal + DX_BL + DX_nearest_1.0_ordinal*Time',ranef='(Time|RID)',
             data=adniLongNoNA, adjustBatch = FALSE, orderby = 'var', xlabel = 'Field Strength')

batchBoxplot(idvar='RID',batchvar='FLDSTRENG',feature = 'ST111TA',
             formula='AGE + PTGENDER_ordinal + DX_BL + DX_nearest_1.0_ordinal*Time',ranef='(Time|RID)',
             data=adniLongNoNA, adjustBatch = TRUE, orderby = 'var', xlabel = 'Field Strength')

trajPlot(idvar='RID',timevar='Time',batchvar='FLDSTRENG',feature = 'ST111TA',
         data=adniLongNoNA,
         line.col=adniLongNoNA$DX_nearest_1.0_ordinal[!duplicated(adniLongNoNA$RID)]+1) # everyone

trajPlot(idvar='RID',timevar='Time',batchvar='FLDSTRENG',feature = 'ST111TA.combat',
         data=na.omit(adniLongHarmonized),
         line.col=adniLongNoNA$DX_nearest_1.0_ordinal[!duplicated(adniLongNoNA$RID)]+1) # everyone

adniLongHarmonized <- subset(adniLongHarmonized, select=c('RID','Time',rois.combat))

# Merge into full dataframe
library(dplyr)
adniLongHarmonizedFull <- left_join(adniLong,adniLongHarmonized, by=c('RID','Time'),relationship='many-to-many')

# Save CSV of harmonized output 
write.csv(adniLongHarmonizedFull,'data/ADNI/MRI_DNA_COG_LONG_H.csv')

