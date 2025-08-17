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
tertiles <- quantile(adni$ShirebyG2020_res_nearest, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)

# Assign labels
adni$GroupShireby2020 <- cut(adni$ShirebyG2020_res_nearest, 
                       breaks = tertiles, 
                       labels = c("Decelerated", "Neutral", "Accelerated"), 
                       include.lowest = TRUE)

# HORVATH 2013
tertiles <- quantile(adni$HorvathS2013_res_nearest, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)

# Assign labels
adni$GroupHorvath2013 <- cut(adni$HorvathS2013_res_nearest, 
                             breaks = tertiles, 
                             labels = c("Decelerated", "Neutral", "Accelerated"), 
                             include.lowest = TRUE)

# HORVATH 2018
tertiles <- quantile(adni$HorvathS2018_res_nearest, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)

# Assign labels
adni$GroupHorvath2018 <- cut(adni$HorvathS2018_res_nearest, 
                             breaks = tertiles, 
                             labels = c("Decelerated", "Neutral", "Accelerated"), 
                             include.lowest = TRUE)

# HANNUM 2013
tertiles <- quantile(adni$HannumG2013_res_nearest, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)

# Assign labels
adni$GroupHannum2013 <- cut(adni$HannumG2013_res_nearest, 
                             breaks = tertiles, 
                             labels = c("Decelerated", "Neutral", "Accelerated"), 
                             include.lowest = TRUE)

# DUNEDIN
tertiles <- quantile(adni$DunedinPACE_nearest, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)

# Assign labels
adni$GroupDunedinPACE <- cut(adni$DunedinPACE_nearest, 
                            breaks = tertiles, 
                            labels = c("Decelerated", "Neutral", "Accelerated"), 
                            include.lowest = TRUE)

# AVG CLOCK
adni <- adni %>%
  rowwise() %>%
  mutate(AgeClockCombo = mean(c(HorvathS2013_nearest,HannumG2013_nearest,ShirebyG2020_nearest)),
         AgeAccClockCombo = mean(c(HorvathS2013_res_nearest,HannumG2013_res_nearest,ShirebyG2020_res_nearest)))

tertiles <- quantile(adni$AgeAccClockCombo, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)

# Assign labels
adni$GroupClockCombo <- cut(adni$AgeAccClockCombo, 
                            breaks = tertiles, 
                            labels = c("Decelerated", "Neutral", "Accelerated"), 
                            include.lowest = TRUE)

write.csv(adni,'data/ADNI/MRI_DNA_COG_XSH_COMET_Grouped.csv')
