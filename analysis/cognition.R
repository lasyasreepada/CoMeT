library(dplyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(EnvStats)
library(tableone)
library(stringr)

setwd('/Users/lasyasreepada/Projects/CoMeT/')
adni <- read.csv('data/ADNI/MRI_DNA_COG_XSH_COMET_Grouped.csv')

cu <- adni[adni$Case==0,]
mci <- adni[adni$DX_nearest_1.0=='MCI',]
ad <- adni[adni$DX_nearest_1.0=='Dementia',]
sym <- adni[adni$Case==1,]

# Factor some variables
adni$Case <- as.factor(ifelse(adni$Case == 0, "Cognitively Unimpaired", "Symptomatic")) 
adni$APOE4 <- as.factor(adni$APOE4)
adni$GroupChronAge <- factor(adni$GroupChronAge, levels=c("<= 65","Average",">= 80"))
adni$GroupClockCombo <- factor(adni$GroupClockCombo, levels=c("Decelerated","Neutral","Accelerated"))

# Create tables by age groups
# Variables to summarize
myVars <- c('PHC_EXF_nearest_wscore', 'PHC_LAN_nearest_wscore', 'PHC_VSP_nearest_wscore', 'AVRECTOT_nearest_wscore')

# CU vs Symptomatic
tab <- CreateTableOne(vars = myVars, strata = 'Case', data = adni)

# Chronological Age Groups in Symptomatic
tab_sym <- CreateTableOne(vars = myVars, strata = 'GroupChronAge', data = sym)

# Biological Age Groups in Symptomatic
tab_sym_bio <- CreateTableOne(vars = myVars, strata = 'GroupDunedinPACE', data = sym)

# CoMeT
myVars <- c('CoMeT_EXFb','CoMeT_LANb','CoMeT_VSPb')

# Chronological Age Groups in Symptomatic
tab_sym_comet <- CreateTableOne(vars = myVars, strata = 'GroupChronAge', data = sym)

# Biological Age Groups in Symptomatic
tab_sym_comet_bio <- CreateTableOne(vars = myVars, strata = 'GroupClockCombo', data = sym)

# Compute mean cognitive score
adni <- adni %>%
  rowwise() %>%
  mutate(PHC_mean = mean(c(PHC_EXF_nearest_Zscore,PHC_MEM_nearest_Zscore,PHC_LAN_nearest_Zscore,PHC_VSP_nearest_Zscore)))

# Compute delta scores per domain
adni <- adni %>%
  rowwise() %>%
  mutate(EXF_delta = abs(PHC_EXF_nearest_Zscore - PHC_mean),
         MEM_delta = abs(PHC_MEM_nearest_Zscore - PHC_mean),
         LAN_delta = PHC_LAN_nearest_Zscore - PHC_mean,
         VSP_delta = PHC_VSP_nearest_Zscore - PHC_mean)

# Deltas by Chronological Age
ggplot(adni[(adni$GroupChronAge!='Average' & !is.na(adni$Case)),],aes(as.factor(GroupChronAge),EXF_delta)) +
  geom_violin(outlier.shape=NA,show.legend = FALSE) + stat_summary(fun='mean',geom='crossbar',width=0.5,color='blue3') +facet_grid(cols = vars(Case)) + geom_jitter(aes(color=as.factor(Case)),size=0.7, alpha=0.6, width = 0.25,show.legend = FALSE) +
  stat_n_text() + labs(title="Executive Function Delta by Chronological Age Groups",
                       x ="Chronological Age Group", y = "Executive Function Delta", fill = "Diagnosis") + theme_pubclean()

ggplot(adni[(adni$GroupChronAge!='Average' & !is.na(adni$Case)),],aes(as.factor(GroupChronAge),LAN_delta)) +
  geom_violin(outlier.shape=NA,show.legend = FALSE) + stat_summary(fun='mean',geom='crossbar',width=0.5,color='blue3') +facet_grid(cols = vars(Case)) + geom_jitter(aes(color=as.factor(Case)),size=0.7, alpha=0.6, width = 0.25,show.legend = FALSE) +
  stat_n_text() + labs(title="Language Delta by Chronological Age Groups",
                       x ="Chronological Age Group", y = "Language Delta", fill = "Diagnosis") + theme_pubclean()

ggplot(adni[(adni$GroupChronAge!='Average' & !is.na(adni$Case)),],aes(as.factor(GroupChronAge),VSP_delta)) +
  geom_violin(outlier.shape=NA,show.legend = FALSE) + stat_summary(fun='mean',geom='crossbar',width=0.5,color='blue3') +facet_grid(cols = vars(Case)) + geom_jitter(aes(color=as.factor(Case)),size=0.7, alpha=0.6, width = 0.25,show.legend = FALSE) +
  stat_n_text() + labs(title="Visuospatial Delta by Chronological Age Groups",
                       x ="Chronological Age Group", y = "Visuospatial Delta", fill = "Diagnosis") + theme_pubclean()

ggplot(adni[(adni$GroupChronAge!='Average' & !is.na(adni$Case)),],aes(as.factor(GroupChronAge),MEM_delta)) +
  geom_violin(outlier.shape=NA,show.legend = FALSE) + stat_summary(fun='mean',geom='crossbar',width=0.5,color='blue3') +facet_grid(cols = vars(Case)) + geom_jitter(aes(color=as.factor(Case)),size=0.7, alpha=0.6, width = 0.25,show.legend = FALSE) +
  stat_n_text() + labs(title="Memory Delta by Chronological Age Groups",
                       x ="Chronological Age Group", y = "Memory Delta", fill = "Diagnosis") + theme_pubclean()

# Deltas by Biological Age
ggplot(adni[((!is.na(adni$GroupClockCombo) & adni$Case=='Symptomatic')),],aes(as.factor(GroupClockCombo),EXF_delta)) +
  geom_boxplot(outlier.shape=NA,show.legend = FALSE) + geom_jitter(aes(color=as.factor(Case)),size=0.7, alpha=0.6, width = 0.25,show.legend = FALSE) +
  stat_n_text() + labs(title="Executive Function Delta by Biological Age Groups",
                       x ="Biological Age Group", y = "Executive Function Delta", fill = "Diagnosis") 

ggplot(adni[((!is.na(adni$GroupClockCombo) & adni$Case=='Symptomatic')),],aes(as.factor(GroupClockCombo),LAN_delta)) +
  geom_boxplot(outlier.shape=NA,show.legend = FALSE) + geom_jitter(aes(color=as.factor(Case)),size=0.7, alpha=0.6, width = 0.25,show.legend = FALSE) +
  stat_n_text() + labs(title="Language Delta by Biological Age Groups",
                       x ="Biological Age Group", y = "Language Delta", fill = "Diagnosis") 

ggplot(adni[((!is.na(adni$GroupClockCombo) & adni$Case=='Symptomatic')),],aes(as.factor(GroupClockCombo),VSP_delta)) +
  geom_boxplot(outlier.shape=NA,show.legend = FALSE) + geom_jitter(aes(color=as.factor(Case)),size=0.7, alpha=0.6, width = 0.25,show.legend = FALSE) +
  stat_n_text() + labs(title="Visuospatial Delta by Biological Age Groups",
                       x ="Biological Age Group", y = "Visuospatial Delta", fill = "Diagnosis") 

ggplot(adni[((!is.na(adni$GroupClockCombo) & adni$Case=='Symptomatic')),],aes(as.factor(GroupClockCombo),MEM_delta)) +
  geom_boxplot(outlier.shape=NA,show.legend = FALSE) + geom_jitter(aes(color=as.factor(Case)),size=0.7, alpha=0.6, width = 0.25,show.legend = FALSE) +
  stat_n_text() + labs(title="Memory Delta by Biological Age Groups",
                       x ="Biological Age Group", y = "Memory Delta", fill = "Diagnosis") 

# STATS

# Groups for Comparison
exfN <- na.omit(as.matrix(adni[(adni$GroupClockCombo=='Negatively Accelerated' & adni$Case=='Symptomatic'),'EXF_delta']))
exfP <- na.omit(as.matrix(adni[(adni$GroupClockCombo=='Positively Accelerated' & adni$Case=='Symptomatic'),'EXF_delta']))

lanN <- na.omit(as.matrix(adni[(adni$GroupClockCombo=='Negatively Accelerated' & adni$Case=='Symptomatic'),'LAN_delta']))
lanP <- na.omit(as.matrix(adni[(adni$GroupClockCombo=='Positively Accelerated' & adni$Case=='Symptomatic'),'LAN_delta']))

vspN <- na.omit(as.matrix(adni[(adni$GroupClockCombo=='Negatively Accelerated' & adni$Case=='Symptomatic'),'VSP_delta']))
vspP <- na.omit(as.matrix(adni[(adni$GroupClockCombo=='Positively Accelerated' & adni$Case=='Symptomatic'),'VSP_delta']))

memN <- na.omit(as.matrix(adni[(adni$GroupClockCombo=='Negatively Accelerated' & adni$Case=='Symptomatic'),'MEM_delta']))
memP <- na.omit(as.matrix(adni[(adni$GroupClockCombo=='Positively Accelerated' & adni$Case=='Symptomatic'),'MEM_delta']))

# Wilcox signed rank tests (one-tail)
res.exf.combo <- wilcox.test(exfN, exfP, paired=FALSE, exact = FALSE)
res.lan.combo <- wilcox.test(lanN, lanP, paired=FALSE, exact = FALSE)
res.vsp.combo <- wilcox.test(vspN, vspP, paired=FALSE, exact = FALSE)
res.mem.combo <- wilcox.test(memN, memP, paired=FALSE, exact = FALSE)

# Frequency Tables
cases <- adni[(adni$Case=='Patient') & !(adni$GroupClockCombo == 'Neutral'),]

cases$mask.vs <- cases$PHC_VSP_Zscore < -1
table(cases$GroupClockCombo,cases$mask.vs)

cases$mask.mem <- cases$PHC_MEM_Zscore < -1
table(cases$GroupClockCombo,cases$mask.mem)

cases$mask.exf <- cases$PHC_EXF_Zscore < -1
table(cases$GroupClockCombo,cases$mask.exf)

cases$mask.lan <- cases$PHC_LAN_Zscore < -1
table(cases$GroupClockCombo,cases$mask.lan)

cases <- adni[(adni$Case==1) & !(adni$GroupDunedinPACE == 'Neutral'),]

cases$mask.vs <- cases$PHC_VSP_Zscore < -0.5
table(cases$GroupDunedinPACE,cases$mask.vs)

cases$mask.mem <- cases$AVRECTOT_Zscore < -0.5
table(cases$GroupDunedinPACE,cases$mask.mem)

cases$mask.exf <- cases$PHC_EXF_Zscore < -0.5
table(cases$GroupDunedinPACE,cases$mask.exf)

cases$mask.lan <- cases$PHC_LAN_Zscore < -0.5
table(cases$GroupDunedinPACE,cases$mask.lan)

cases <- adni[(adni$Case==1) & !(adni$GroupChronAge == 'Average'),]

cases$mask.vs <- cases$PHC_VSP_Zscore < -1
table(cases$GroupChronAge,cases$mask.vs)

cases$mask.mem <- cases$AVRECTOT_Zscore < -1
table(cases$GroupChronAge,cases$mask.mem)

cases$mask.exf <- cases$PHC_EXF_Zscore < -1
table(cases$GroupChronAge,cases$mask.exf)

cases$mask.lan <- cases$PHC_LAN_Zscore < -1
table(cases$GroupChronAge,cases$mask.lan)

