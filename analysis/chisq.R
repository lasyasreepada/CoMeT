library(dplyr)
library(ggplot2)
library(ggpubr)
library(EnvStats)
library(tableone)

setwd('/Users/lasyasreepada/Projects/CoMeT/')
adni <- read.csv('data/ADNI/MRI_DNA_COG_XSH_COMET_Grouped.csv')

# Factor some variables
adni$Case <- as.factor(ifelse(adni$Case == 0, "Cognitively Unimpaired", "Symptomatic")) 
adni$APOE4 <- as.factor(adni$APOE4)
adni$GroupChronAge <- factor(adni$GroupChronAge, levels=c("<= 65","Average",">= 80"))
adni$GroupClockCombo <- factor(adni$GroupClockCombo, levels=c("Decelerated","Neutral","Accelerated"))
adni$PTEDUCAT_binned <- cut(adni$PTEDUCAT, breaks=c(0, 12, 16, max(adni$PTEDUCAT)), right = TRUE,labels=c('hs','college','grad'))

# Extract subsets
cn <- adni[(!is.na(adni$GroupClockCombo) & adni$Case=='Cognitively Unimpaired'),]
cases <- adni[(!is.na(adni$GroupClockCombo) & adni$Case=='Symptomatic'),]

# X SQ TESTS

# Control
apoe.cn <- table(cn$GroupClockCombo,cn$APOE4) # APOE4
edu.cn <- table(cn$GroupClockCombo,cn$PTEDUCAT_binned) # Education

chisq.apoe.cn <- chisq.test(apoe.cn)
chisq.edu.cn <- chisq.test(edu.cn)

# Case
apoe.case <- table(cases$GroupClockCombo,cases$APOE4) # APOE4
edu.case <- table(cases$GroupClockCombo,cases$PTEDUCAT_binned) # Education

chisq.apoe.case <- chisq.test(apoe.case)
chisq.edu.case <- chisq.test(edu.case)
