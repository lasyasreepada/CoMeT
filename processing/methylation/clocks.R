# Libraries
library(dnaMethyAge)
library(SummarizedExperiment)

# Set WD
setwd('/Users/lasyasreepada/Projects/CoMeT/')

# Read in betas matrix
adni <- readRDS("data/ADNI/methylation/methylation.rds")

# Extract betas
betas <- assays(adni)[[1]]

# Extract colData as its own dataframe
coldata <- as.data.frame(colData(adni))

# Compute clock ages and age accelerations
HorvathS2013 <- methyAge(betas, clock='HorvathS2013', age_info=coldata,simple_mode = TRUE)
HorvathS2018 <- methyAge(betas, clock='HorvathS2018', age_info=coldata,simple_mode = TRUE)
HannumG2013 <- methyAge(betas, clock="HannumG2013", age_info=coldata)
LevineM2018 <- methyAge(betas, clock="LevineM2018", age_info=coldata)
ShirebyG2020 <- methyAge(betas, clock="ShirebyG2020", age_info=coldata)
DunedinPACE <- methyAge(betas, clock="DunedinPACE", age_info=coldata)
ZhangQ2019 <- methyAge(betas, clock='ZhangQ2019', age_info=coldata)

# Add clock ages and residuals to dataframe
# Horvath 1
coldata$HorvathS2013 <- HorvathS2013$mAge
coldata$HorvathS2013_res <- HorvathS2013$Age_Acceleration

# Horvath 2
coldata$HorvathS2018 <- HorvathS2018$mAge
coldata$HorvathS2018_res <- HorvathS2018$Age_Acceleration

# Hannum
coldata$HannumG2013 <- HannumG2013$mAge
coldata$HannumG2013_res <- HannumG2013$Age_Acceleration

# Levine (Pheno)
coldata$LevineM2018 <- LevineM2018$mAge
coldata$LevineM2018_res <- LevineM2018$Age_Acceleration

# Shireby (Cortical)
coldata$ShirebyG2020 <- ShirebyG2020$mAge
coldata$ShirebyG2020_res <- ShirebyG2020$Age_Acceleration

# DunedinPACE
coldata$DunedinPACE <- DunedinPACE$mAge

# Lu (Telomere)
coldata$ZhangQ2019 <- ZhangQ2019$mAge
coldata$ZhangQ2019_res <- ZhangQ2019$Age_Acceleration


# Write to CSV
write.csv(coldata,'data/ADNI/methylation/clocks.csv')
