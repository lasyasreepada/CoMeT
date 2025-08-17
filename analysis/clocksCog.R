library(dplyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(EnvStats)
library(tableone)
library(stringr)

setwd('/Users/lasyasreepada/Projects/CoMeT/')
adni <- read.csv('data/ADNI/MRI_DNA_COG_XSH_COMET_Grouped.csv')

# Factor some variables
adni$Case <- as.factor(ifelse(adni$Case == 0, "Cognitively Unimpaired", "Symptomatic")) 
adni$APOE4 <- as.factor(adni$APOE4)
adni$GroupChronAge <- factor(adni$GroupChronAge, levels=c("<= 65","Average",">= 80"))
adni$GroupClockCombo <- factor(adni$GroupClockCombo, levels=c("Decelerated","Neutral","Accelerated"))
adni$GroupDunedinPACE <- factor(adni$GroupDunedinPACE, levels=c("Decelerated","Neutral","Accelerated"))

# Extract clinical dx groups
cn <- adni[adni$Case=='Cognitively Unimpaired',]
case <- adni[adni$Case=='Symptomatic',]
mci <- adni[adni$DX_nearest_1.0=='MCI',]
ad <- adni[adni$DX_nearest_1.0=='Dementia',]

mycolor <- c("#011F5B", "#A8A8A8", "#990000") 

# 4A
# 
# ggplot(adni,aes(factor(GroupChronAge),PHC_EXF_nearest)) +
#   geom_violin(show.legend = FALSE) + geom_jitter(aes(col=Case),size=1, alpha=0.7, width = 0.25,show.legend = FALSE) +
#   stat_summary(fun='mean',geom='crossbar',width=0.5,color='black') + 
#   facet_grid(cols = vars(Case)) +
#   stat_n_text() +
#   stat_compare_means(comparisons = list(c("<= 65",">= 80")), label = "p.signif") +
#   scale_color_manual(values = c("#A8A8A8","#011F5B")) + 
#   labs(title="Executive Function by Chronological Age",x ="Chronological Age Group", y = "Executive Function Score") + 
#   theme_pubclean() 

ggplot(adni,aes(factor(GroupChronAge),PHC_EXF_nearest_wscore)) +
  geom_violin(show.legend = FALSE) + geom_jitter(aes(col=Case),size=1, alpha=0.7, width = 0.25,show.legend = FALSE) +
  stat_summary(fun='mean',geom='crossbar',width=0.5,color='black') + 
  facet_grid(cols = vars(Case)) +
  stat_n_text() +
  stat_compare_means(comparisons = list(c("<= 65",">= 80")), label = "p.signif") +
  scale_color_manual(values = c("#A8A8A8","#011F5B")) + 
  labs(title="Executive Function by Chronological Age",x ="Chronological Age Group", y = "Executive Function w-score") + 
  theme_pubclean() 

# ggplot(adni,aes(factor(GroupChronAge),PHC_LAN_nearest)) +
#   geom_violin(show.legend = FALSE) + geom_jitter(aes(col=Case),size=1, alpha=0.7, width = 0.25,show.legend = FALSE) +
#   stat_summary(fun='mean',geom='crossbar',width=0.5,color='black') + 
#   facet_grid(cols = vars(Case)) +
#   stat_n_text() +
#   stat_compare_means(comparisons = list(c("<= 65",">= 80")), label = "p.signif") +
#   scale_color_manual(values = c("#A8A8A8","#011F5B")) + 
#   labs(title="Language by Chronological Age",x ="Chronological Age Group", y = "Language Score") + 
#   theme_pubclean() 

ggplot(adni,aes(factor(GroupChronAge),PHC_LAN_nearest_wscore)) +
  geom_violin(show.legend = FALSE) + geom_jitter(aes(col=Case),size=1, alpha=0.7, width = 0.25,show.legend = FALSE) +
  stat_summary(fun='mean',geom='crossbar',width=0.5,color='black') + 
  facet_grid(cols = vars(Case)) +
  stat_n_text() +
  stat_compare_means(comparisons = list(c("<= 65",">= 80")), label = "p.signif") +
  scale_color_manual(values = c("#A8A8A8","#011F5B")) + 
  labs(title="Language by Chronological Age",x ="Chronological Age Group", y = "Language w-score") + 
  theme_pubclean() 

# ggplot(adni,aes(factor(GroupChronAge),PHC_VSP_nearest)) +
#   geom_violin(show.legend = FALSE) + geom_jitter(aes(col=Case),size=1, alpha=0.7, width = 0.25,show.legend = FALSE) +
#   stat_summary(fun='mean',geom='crossbar',width=0.5,color='black') + 
#   facet_grid(cols = vars(Case)) +
#   stat_n_text() +
#   stat_compare_means(comparisons = list(c("<= 65",">= 80")), label = "p.signif") +
#   scale_color_manual(values = c("#A8A8A8","#011F5B")) + 
#   labs(title="Visuospatial by Chronological Age",x ="Chronological Age Group", y = "Visuospatial Score") + 
#   theme_pubclean() 

ggplot(adni,aes(factor(GroupChronAge),PHC_VSP_nearest_wscore)) +
  geom_violin(show.legend = FALSE) + geom_jitter(aes(col=Case),size=1, alpha=0.7, width = 0.25,show.legend = FALSE) +
  stat_summary(fun='mean',geom='crossbar',width=0.5,color='black') + 
  facet_grid(cols = vars(Case)) +
  stat_n_text() +
  stat_compare_means(comparisons = list(c("<= 65",">= 80")), label = "p.signif") +
  scale_color_manual(values = c("#A8A8A8","#011F5B")) + 
  labs(title="Visuospatial by Chronological Age",x ="Chronological Age Group", y = "Visuospatial w-score") + 
  theme_pubclean() 

# ggplot(adni,aes(factor(GroupChronAge),PHC_MEM_nearest)) +
#   geom_violin(show.legend = FALSE) + geom_jitter(aes(col=Case),size=1, alpha=0.7, width = 0.25,show.legend = FALSE) +
#   stat_summary(fun='mean',geom='crossbar',width=0.5,color='black') + 
#   facet_grid(cols = vars(Case)) +
#   stat_n_text() +
#   stat_compare_means(comparisons = list(c("<= 65",">= 80")), label = "p.signif") +
#   scale_color_manual(values = c("#A8A8A8","#990000")) + 
#   labs(title="Memory by Chronological Age",x ="Chronological Age Group", y = "Memory Score") + 
#   theme_pubclean() 

# ggplot(adni,aes(factor(GroupChronAge),PHC_MEM_nearest_wscore)) +
#   geom_violin(show.legend = FALSE) + geom_jitter(aes(col=Case),size=1, alpha=0.7, width = 0.25,show.legend = FALSE) +
#   stat_summary(fun='mean',geom='crossbar',width=0.5,color='black') + 
#   facet_grid(cols = vars(Case)) +
#   stat_n_text() +
#   stat_compare_means(comparisons = list(c("<= 65",">= 80")), label = "p.signif") +
#   scale_color_manual(values = c("#A8A8A8","#990000")) + 
#   labs(title="Memory by Chronological Age",x ="Chronological Age Group", y = "Memory w-score") + 
#   theme_pubclean() 

ggplot(adni,aes(factor(GroupChronAge),AVRECTOT_nearest_wscore)) +
  geom_violin(show.legend = FALSE) + geom_jitter(aes(col=Case),size=1, alpha=0.7, width = 0.25,show.legend = FALSE) +
  stat_summary(fun='mean',geom='crossbar',width=0.5,color='black') + 
  facet_grid(cols = vars(Case)) +
  stat_n_text() +
  stat_compare_means(comparisons = list(c("<= 65",">= 80")), label = "p.signif") +
  scale_color_manual(values = c("#A8A8A8","#990000")) + 
  labs(title="Recall Memory by Chronological Age",x ="Chronological Age Group", y = "Recall Memory w-score") + 
  theme_pubclean()

ggplot(adni,aes(factor(GroupChronAge),MIN_wscore)) +
  geom_violin(show.legend = FALSE) + geom_jitter(aes(col=Case),size=1, alpha=0.7, width = 0.25,show.legend = FALSE) +
  stat_summary(fun='mean',geom='crossbar',width=0.5,color='black') + 
  facet_grid(cols = vars(Case)) +
  stat_n_text() +
  stat_compare_means(comparisons = list(c("<= 65",">= 80")), label = "p.signif") +
  scale_color_manual(values = c("#A8A8A8","#990000")) + 
  labs(title="MIN Cortical by Chronological Age",x ="Chronological Age Group", y = "MIN Cortical w-score") + 
  theme_pubclean()

ggplot(adni,aes(factor(GroupChronAge),CoMeT_EXF)) +
  geom_violin(show.legend = FALSE) + geom_jitter(aes(col=Case),size=1, alpha=0.7, width = 0.25,show.legend = FALSE) +
  stat_summary(fun='mean',geom='crossbar',width=0.5,color='black') + 
  facet_grid(cols = vars(Case)) +
  stat_n_text() +
  stat_compare_means(comparisons = list(c("<= 65",">= 80")), label = "p.signif") +
  scale_color_manual(values = c("#A8A8A8","#660099")) + 
  labs(title="CoMeT EXF by Chronological Age",x ="Chronological Age Group", y = "CoMeT EXF") + 
  theme_pubclean() 

ggplot(adni,aes(factor(GroupChronAge),CoMeT_LAN)) +
  geom_violin(show.legend = FALSE) + geom_jitter(aes(col=Case),size=1, alpha=0.7, width = 0.25,show.legend = FALSE) +
  stat_summary(fun='mean',geom='crossbar',width=0.5,color='black') + 
  facet_grid(cols = vars(Case)) +
  stat_n_text() +
  stat_compare_means(comparisons = list(c("<= 65",">= 80")), label = "p.signif") +
  scale_color_manual(values = c("#A8A8A8","#660099")) + 
  labs(title="CoMeT LAN by Chronological Age",x ="Chronological Age Group", y = "CoMeT LAN") + 
  theme_pubclean() 

ggplot(adni,aes(factor(GroupChronAge),CoMeT_VSP)) +
  geom_violin(show.legend = FALSE) + geom_jitter(aes(col=Case),size=1, alpha=0.7, width = 0.25,show.legend = FALSE) +
  stat_summary(fun='mean',geom='crossbar',width=0.5,color='black') + 
  facet_grid(cols = vars(Case)) +
  stat_n_text() +
  stat_compare_means(comparisons = list(c("<= 65",">= 80")), label = "p.signif") +
  scale_color_manual(values = c("#A8A8A8","#660099")) + 
  labs(title="CoMeT VSP by Chronological Age",x ="Chronological Age Group", y = "CoMeT VSP") + 
  theme_pubclean() 

ggplot(adni,aes(factor(GroupChronAge),CoMeT_MIN)) +
  geom_violin(show.legend = FALSE) + geom_jitter(aes(col=Case),size=1, alpha=0.7, width = 0.25,show.legend = FALSE) +
  stat_summary(fun='mean',geom='crossbar',width=0.5,color='black') + 
  facet_grid(cols = vars(Case)) +
  stat_n_text() +
  stat_compare_means(comparisons = list(c("<= 65",">= 80")), label = "p.signif") +
  scale_color_manual(values = c("#A8A8A8","#660099")) + 
  labs(title="CoMeT MIN Cortical by Chronological Age",x ="Chronological Age Group", y = "CoMeT MIN") + 
  theme_pubclean() 

# 4B
ggplot(adni[adni$Case=='Symptomatic',],aes(AGE,CoMeT_EXF)) +
  geom_point() + 
  geom_smooth(method = 'lm') +
  stat_cor(method='pearson') +
  labs(title="CoMeT vs Chronological Age - Symptomatic",x ="Chronological Age (years)", y = "CoMeT") +
  theme_pubclean()

ggplot(adni[adni$Case=='Symptomatic',],aes(AGE,CoMeT_LAN)) +
  geom_point() + 
  geom_smooth(method = 'lm') +
  stat_cor(method='pearson') +
  labs(title="CoMeT vs Chronological Age - Symptomatic",x ="Chronological Age (years)", y = "CoMeT") +
  theme_pubclean()

ggplot(adni[adni$Case=='Symptomatic',],aes(AGE,CoMeT_VSP)) +
  geom_point() + 
  geom_smooth(method = 'lm') +
  stat_cor(method='pearson') +
  labs(title="CoMeT vs Chronological Age - Symptomatic",x ="Chronological Age (years)", y = "CoMeT") +
  theme_pubclean()

ggplot(adni[adni$Case=='Symptomatic',],aes(AGE,CoMeT_MIN)) +
  geom_point() + 
  geom_smooth(method = 'lm') +
  stat_cor(method='pearson') +
  labs(title="CoMeT vs Chronological Age - Symptomatic",x ="Chronological Age (years)", y = "CoMeT") +
  theme_pubclean()

ggplot(ad[ad$Case=='Symptomatic',],aes(AGE,CoMeT_EXF)) +
  geom_point() + 
  geom_smooth(method = 'lm') +
  stat_cor(method='pearson') +
  labs(title="CoMeT vs Chronological Age - Dementia",x ="Chronological Age (years)", y = "CoMeT") +
  theme_pubclean()

ggplot(ad[ad$Case=='Symptomatic',],aes(AGE,CoMeT_LAN)) +
  geom_point() + 
  geom_smooth(method = 'lm') +
  stat_cor(method='pearson') +
  labs(title="CoMeT vs Chronological Age - Dementia",x ="Chronological Age (years)", y = "CoMeT") +
  theme_pubclean()

ggplot(ad[ad$Case=='Symptomatic',],aes(AGE,CoMeT_VSP)) +
  geom_point() + 
  geom_smooth(method = 'lm') +
  stat_cor(method='pearson') +
  labs(title="CoMeT vs Chronological Age - Dementia",x ="Chronological Age (years)", y = "CoMeT") +
  theme_pubclean()

# 5A
p1 <- ggplot(adni[((!is.na(adni$GroupClockCombo) & adni$Case=='Symptomatic')),],aes(GroupClockCombo,CoMeT_EXF)) +
  geom_violin(aes(color=GroupClockCombo),show.legend = FALSE,trim = FALSE) +
  stat_summary(fun='mean',geom='crossbar',width=0.6,color='black') +
  geom_jitter(aes(color=GroupClockCombo),size=1, alpha=0.7, width=0.25, height=0, show.legend = FALSE) +
  stat_n_text() + 
  stat_compare_means(comparisons = list(c("Decelerated","Accelerated")), label = "p.signif") +
  labs(title="CoMeT Executive Function - Symptomatic",x ="Biological Age Group", y = "CoMeT Executive Function-Memory") + 
  scale_color_manual(values = mycolor) + 
  theme_pubclean()

p2 <- ggplot(adni[((!is.na(adni$GroupClockCombo) & adni$Case=='Symptomatic')),],aes(GroupClockCombo,CoMeT_LAN)) +
  geom_violin(aes(color=GroupClockCombo),show.legend = FALSE,trim = FALSE) +
  stat_summary(fun='mean',geom='crossbar',width=0.6,color='black') +
  geom_jitter(aes(color=GroupClockCombo),size=1, alpha=0.7, width=0.25, height=0, show.legend = FALSE) +
  stat_n_text() + 
  stat_compare_means(comparisons = list(c("Decelerated","Accelerated")), label = "p.signif") +
  labs(title="CoMeT Language - Symptomatic",x ="Biological Age Group", y = "CoMeT Language-Memory") + 
  scale_color_manual(values = mycolor) + 
  theme_pubclean()

p3 <- ggplot(adni[((!is.na(adni$GroupClockCombo) & adni$Case=='Symptomatic')),],aes(GroupClockCombo,CoMeT_VSP)) +
  geom_violin(aes(color=GroupClockCombo),show.legend = FALSE,trim = FALSE) +
  stat_summary(fun='mean',geom='crossbar',width=0.6,color='black') +
  geom_jitter(aes(color=GroupClockCombo),size=1, alpha=0.7, width=0.25, height=0, show.legend = FALSE) +
  stat_n_text() + 
  stat_compare_means(comparisons = list(c("Decelerated","Accelerated")), label = "p.signif") +
  labs(title="CoMeT Visuospatial - Symptomatic",x ="Biological Age Group", y = "CoMeT Visuospatial-Memory") + 
  scale_color_manual(values = mycolor) + 
  theme_pubclean()

ggplot(adni[((!is.na(adni$GroupClockCombo) & adni$Case=='Symptomatic')),],aes(GroupClockCombo,CoMeT_MIN)) +
  geom_violin(aes(color=GroupClockCombo),show.legend = FALSE,trim = FALSE) +
  stat_summary(fun='mean',geom='crossbar',width=0.6,color='black') +
  geom_jitter(aes(color=GroupClockCombo),size=1, alpha=0.7, width=0.25, height=0, show.legend = FALSE) +
  stat_n_text() + 
  stat_compare_means(comparisons = list(c("Decelerated","Accelerated")), label = "p.signif") +
  labs(title="CoMeT MIN - Symptomatic",x ="Biological Age Group", y = "CoMeT MIN-Memory") + 
  scale_color_manual(values = mycolor) + 
  theme_pubclean()

grid.arrange(p1, p2, p3, nrow = 1)

ggplot(ad[((!is.na(ad$GroupClockCombo) & ad$Case=='Symptomatic')),],aes(GroupClockCombo,CoMeT_EXF)) +
  geom_violin(aes(color=GroupClockCombo),show.legend = FALSE,trim = FALSE) +
  stat_summary(fun='mean',geom='crossbar',width=0.5,color='black') +
  geom_jitter(aes(color=GroupClockCombo),size=1.1, alpha=0.8, width=0.25, height=0, show.legend = FALSE) +
  stat_n_text() +
  stat_compare_means(comparisons = list(c("Decelerated","Accelerated")), label = "p.signif") +
  labs(title="CoMeT EXF by Biological Age Groups - Dementia",x ="Biological Age Group", y = "CoMeT") +
  scale_color_manual(values = mycolor) +
  theme_pubclean()

ggplot(ad[((!is.na(ad$GroupClockCombo) & ad$Case=='Symptomatic')),],aes(GroupClockCombo,CoMeT_LAN)) +
  geom_violin(aes(color=GroupClockCombo),show.legend = FALSE,trim = FALSE) +
  stat_summary(fun='mean',geom='crossbar',width=0.5,color='black') +
  geom_jitter(aes(color=GroupClockCombo),size=1.1, alpha=0.8, width=0.25, height=0, show.legend = FALSE) +
  stat_n_text() +
  stat_compare_means(comparisons = list(c("Decelerated","Accelerated")), label = "p.signif") +
  labs(title="CoMeT LAN by Biological Age Groups - Dementia",x ="Biological Age Group", y = "CoMeT") +
  scale_color_manual(values = mycolor) +
  theme_pubclean()

ggplot(ad[((!is.na(ad$GroupClockCombo) & ad$Case=='Symptomatic')),],aes(GroupClockCombo,CoMeT_VSP)) +
  geom_violin(aes(color=GroupClockCombo),show.legend = FALSE,trim = FALSE) +
  stat_summary(fun='mean',geom='crossbar',width=0.5,color='black') +
  geom_jitter(aes(color=GroupClockCombo),size=1.1, alpha=0.8, width=0.25, height=0, show.legend = FALSE) +
  stat_n_text() +
  stat_compare_means(comparisons = list(c("Decelerated","Accelerated")), label = "p.signif") +
  labs(title="CoMeT VSP by Biological Age Groups - Dementia",x ="Biological Age Group", y = "CoMeT") +
  scale_color_manual(values = mycolor) +
  theme_pubclean()


# 5B
ggplot(adni[adni$Case=='Symptomatic',],aes(AgeAccClockCombo,CoMeT_EXF)) +
  geom_point() + 
  geom_smooth(method='lm') + 
  stat_cor(method='pearson') + 
  labs(title="CoMeT vs BAG - Symptomatic",x ="BAG (years)", y = "CoMeT") +
  theme_pubclean()

ggplot(adni[adni$Case=='Symptomatic',],aes(AgeAccClockCombo,CoMeT_LAN)) +
  geom_point() + 
  geom_smooth(method='lm') + 
  stat_cor(method='pearson') + 
  labs(title="CoMeT vs BAG - Symptomatic",x ="BAG (years)", y = "CoMeT") +
  theme_pubclean()

ggplot(adni[adni$Case=='Symptomatic',],aes(AgeAccClockCombo,CoMeT_VSP)) +
  geom_point() + 
  geom_smooth(method='lm') + 
  stat_cor(method='pearson') + 
  labs(title="CoMeT vs BAG - Symptomatic",x ="BAG (years)", y = "CoMeT") +
  theme_pubclean()

ggplot(ad,aes(AgeAccClockCombo,CoMeT_EXF)) +
  geom_point() + 
  geom_smooth(method='lm') +
  stat_cor(method='pearson') +
  labs(title="CoMeT vs BAG - Dementia",x ="BAG (years)", y = "CoMeT") +
  theme_pubclean()

ggplot(ad,aes(AgeAccClockCombo,CoMeT_LAN)) +
  geom_point() + 
  geom_smooth(method='lm') +
  stat_cor(method='pearson') +
  labs(title="CoMeT vs BAG - Dementia",x ="BAG (years)", y = "CoMeT") +
  theme_pubclean()

ggplot(ad,aes(AgeAccClockCombo,CoMeT_VSP)) +
  geom_point() + 
  geom_smooth(method='lm') +
  stat_cor(method='pearson') +
  labs(title="CoMeT vs BAG - Dementia",x ="BAG (years)", y = "CoMeT") +
  theme_pubclean()
