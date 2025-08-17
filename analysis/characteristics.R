library(dplyr)
library(ggplot2)
library(ggpubr)
library(EnvStats)
library(tableone)

setwd('/Users/lasyasreepada/Projects/CoMeT/')
adni <- read.csv('data/ADNI/MRI_DNA_COG_XSH_COMET_Grouped.csv')

adni <- as.data.frame(adni)

# Factor some variables
adni$Case <- as.factor(ifelse(adni$Case == 0, "Cognitively Unimpaired", "Symptomatic")) 
adni$GroupChronAge <- factor(adni$GroupChronAge, levels=c("<= 65","Average",">= 80"))
adni$GroupClockCombo <- factor(adni$GroupClockCombo, levels=c("Decelerated","Neutral","Accelerated"))

# Extract AD only
mci <- adni[adni$DX_nearest_1.0=="MCI",]
ad <- adni[adni$DX_nearest_1.0=="Dementia",]

# TABLES

# By Groups
cu <- adni[(adni$Case=='Cognitively Unimpaired'),]
sym <- adni[(adni$Case=='Symptomatic'),]

# Methylation
adni$DNAM <- !is.na(adni$Years_to_DNAm)
dnam <- adni[adni$DNAM,]
dnamCU <- adni[(adni$Case=='Cognitively Unimpaired' & adni$DNAM),]
dnamSym <- adni[(adni$Case=='Symptomatic' & adni$DNAM),]
dnamMCI <- dnamSym[dnamSym$DX_nearest_1.0=='MCI']
dnamAD <- dnamSym[dnamSym$DX_nearest_1.0=='Dementia',]

# Variables to summarize
myVars <- c('AGE','PTGENDER','DX_nearest_1.0','APOE4','PTEDUCAT','CDRSB','Years_to_DNAm',"GroupClockCombo",'AgeAccClockCombo','Corticalv4','MTLb','CoMeTv4b')

# Median variables
medianVars <- c('CDRSB','PTEDUCAT')

# Categorical variables that need transformation
catVars <- c('PTGENDER','DX_nearest_1.0','APOE4','GroupClockCombo')

# TableOne object

# Overall 
tab1 <- CreateTableOne(vars = myVars, data = adni, factorVars = catVars)
tab1

# CU
tab1a <- CreateTableOne(vars = myVars, data = cu, factorVars = catVars)
tab1a

# Symptomatic 
tab1b <- CreateTableOne(vars = myVars, data = sym, factorVars = catVars)
tab1b 

# Overall DNA methylation 
# tab4_full <- CreateTableOne(vars = myVars, data = dnam, factorVars = catVars)
# tab5_full <- CreateTableOne(vars = myVars, data = dnamCase, factorVars = catVars)

# Symptomatic DNA methylation
tab5_groups <- CreateTableOne(vars = myVars, data = dnamCase, strata = 'DX_nearest_1.0_ordinal',factorVars = catVars)
tab5_groups_cn <- CreateTableOne(vars = myVars, data = dnamCN,factorVars = catVars)
