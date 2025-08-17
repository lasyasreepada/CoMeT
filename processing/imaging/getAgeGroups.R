library(dplyr)

# Set WD
setwd('/Users/lasyasreepada/Projects/CoMeT/')

# Read in data
adni <- read.csv('data/ADNI/MRI_DNA_COG_XSH_COMET.csv')

# Categorize subjects by chronological age groups
adni <- adni %>% mutate(GroupChronAge = case_when(AGE >= 80 ~ ">= 80",
                                                  AGE <= 65 ~ "<= 65",
                                                  TRUE ~ "Average"))
# SHIREBY CORTICAL
# Compute mean and SD EAA
mean.shireby2020 <- mean(adni$ShirebyG2020_res_nearest,na.rm=T)
sd.shireby2020 <- sd(adni$ShirebyG2020_res_nearest,na.rm=T)

# Compute confidence intervals
upper.shireby2020 <- mean.shireby2020 + 0.5*sd.shireby2020
lower.shireby2020 <- mean.shireby2020 - 0.5*sd.shireby2020

# Categorize subjects by epigenetic age group
adni <- adni %>% mutate(GroupShireby2020 = case_when(ShirebyG2020_res_nearest >= upper.shireby2020 ~ "Accelerated",
                                                     ShirebyG2020_res_nearest <= lower.shireby2020 ~ "Decelerated",
                                                     (ShirebyG2020_res_nearest > lower.shireby2020 & ShirebyG2020_res_nearest < upper.shireby2020) ~ "Neutral"))
# HORVATH 2013
# Compute mean and SD EAA
mean.horvath2013 <- mean(adni$HorvathS2013_res_nearest,na.rm=T)
sd.horvath2013 <- sd(adni$HorvathS2013_res_nearest,na.rm=T)

# Compute confidence intervals
upper.horvath2013 <- mean.horvath2013 + 0.5*sd.horvath2013
lower.horvath2013 <- mean.horvath2013 - 0.5*sd.horvath2013

# Categorize subjects by epigenetic age group
adni <- adni %>% mutate(GroupHorvath2013 = case_when(HorvathS2013_res_nearest >= upper.horvath2013 ~ "Accelerated",
                                                     HorvathS2013_res_nearest <= lower.horvath2013 ~ "Decelerated",
                                                     (HorvathS2013_res_nearest > lower.horvath2013 & HorvathS2013_res_nearest < upper.horvath2013) ~ "Neutral"))
# HORVATH 2018
# Compute mean and SD EAA
mean.horvath2018 <- mean(adni$HorvathS2018_res_nearest,na.rm=T)
sd.horvath2018 <- sd(adni$HorvathS2018_res_nearest,na.rm=T)

# Compute confidence intervals
upper.horvath2018 <- mean.horvath2018 + 0.5*sd.horvath2018
lower.horvath2018 <- mean.horvath2018 - 0.5*sd.horvath2018

# Categorize subjects by epigenetic age group
adni <- adni %>% mutate(GroupHorvath2018 = case_when(HorvathS2018_res_nearest >= upper.horvath2018 ~ "Accelerated",
                                                     HorvathS2018_res_nearest <= lower.horvath2018 ~ "Decelerated",
                                                     (HorvathS2018_res_nearest > lower.horvath2018 & HorvathS2018_res_nearest < upper.horvath2018) ~ "Neutral"))
# HANNUM 2013
# Compute mean and SD EAA
mean.hannum2013 <- mean(adni$HannumG2013_res_nearest,na.rm=T)
sd.hannum2013 <- sd(adni$HannumG2013_res_nearest,na.rm=T)

# Compute confidence intervals
upper.hannum2013 <- mean.hannum2013 + 0.5*sd.hannum2013
lower.hannum2013 <- mean.hannum2013 - 0.5*sd.hannum2013

# Categorize subjects by epigenetic age group
adni <- adni %>% mutate(GroupHannum2013 = case_when(HannumG2013_res_nearest >= upper.hannum2013 ~ "Accelerated",
                                                    HannumG2013_res_nearest <= lower.hannum2013 ~ "Decelerated",
                                                    (HannumG2013_res_nearest > lower.hannum2013 & HannumG2013_res_nearest < upper.hannum2013) ~ "Neutral"))
# DUNEDIN
# Compute mean and SD EAA
mean.dunedin <- mean(adni$DunedinPACE_nearest,na.rm=T)
sd.dunedin <- sd(adni$DunedinPACE_nearest,na.rm=T)

# Compute confidence intervals
upper.dunedin <- mean.dunedin + 0.5*sd.dunedin
lower.dunedin <- mean.dunedin - 0.5*sd.dunedin

# Categorize subjects by epigenetic age group
adni <- adni %>% mutate(GroupDunedinPACE = case_when(DunedinPACE_nearest > upper.dunedin ~ "Accelerated",
                                                     DunedinPACE_nearest < lower.dunedin ~ "Decelerated",
                                                     (DunedinPACE_nearest < upper.dunedin & DunedinPACE_nearest > lower.dunedin) ~ "Neutral"))
# AVG CLOCK
adni <- adni %>%
  rowwise() %>%
  mutate(AgeClockCombo = mean(c(HorvathS2013_nearest,HannumG2013_nearest,ShirebyG2020_nearest)), 
         AgeAccClockCombo = mean(c(HorvathS2013_res_nearest,HannumG2013_res_nearest,ShirebyG2020_res_nearest)))

# model <- lm(AccClockCombo ~ PTGENDER_M, data = adni)
# adni$AgeAccClockCombo <- resid(model)

mean.combo <- mean(adni$AgeAccClockCombo,na.rm=T)
sd.combo <- sd(adni$AgeAccClockCombo,na.rm=T)

upper.combo <- mean.combo + 0.5*sd.combo
lower.combo <- mean.combo - 0.5*sd.combo

adni <- adni %>% mutate(GroupClockCombo = case_when(AgeAccClockCombo > upper.combo ~ "Accelerated",
                                                    AgeAccClockCombo < lower.combo ~ "Decelerated",
                                                    (AgeAccClockCombo < upper.combo & AgeAccClockCombo > lower.combo) ~ "Neutral"))

write.csv(adni,'data/ADNI/MRI_DNA_COG_XSH_COMET_Grouped.csv')
