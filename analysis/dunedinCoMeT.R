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
adni$GroupDunedinPACE <- factor(adni$GroupDunedinPACE, levels=c("Decelerated","Neutral","Accelerated"))
adni$GroupClockCombo <- factor(adni$GroupClockCombo, levels=c("Decelerated","Neutral","Accelerated"))

# Extract clinical dx groups
cn <- adni[adni$Case=='Cognitively Unimpaired',]
case <- adni[adni$Case=='Symptomatic',]
mci <- adni[adni$DX_nearest_1.0=='MCI',]
ad <- adni[adni$DX_nearest_1.0=='Dementia',]

mycolor <- c("#011F5B", "#A8A8A8", "#990000") 

# 3
mean.combo <- mean(adni$DunedinPACE_nearest,na.rm=T)
sd.combo <- sd(adni$DunedinPACE_nearest,na.rm=T)

upper.combo <- mean.combo + 0.5*sd.combo
lower.combo <- mean.combo - 0.5*sd.combo

p3 <- ggplot(adni,aes(AGE,DunedinPACE_nearest)) +
  geom_point(aes(color=GroupDunedinPACE),na.rm = TRUE) + 
  geom_hline(yintercept = lower.combo,linetype = "dashed") +
  geom_hline(yintercept = upper.combo,linetype = "dashed") +
  labs(title="BAG vs Chronological Age", x ="Chronological Age (years)", y = "BAG (years)", fill = "Biological Age Group") +
  scale_color_manual(values = mycolor, na.translate = F) + 
  theme_pubclean() 

grid.arrange(p3, nrow = 1)

# 4A

p4ai <- ggplot(adni,aes(factor(GroupChronAge),Corticalv4)) +
  geom_violin(show.legend = FALSE) + geom_jitter(aes(col=Case),size=1, alpha=0.7, width = 0.25,show.legend = FALSE) +
  stat_summary(fun='mean',geom='crossbar',width=0.75,color='black') + 
  facet_grid(cols = vars(Case)) +
  stat_n_text() +
  stat_compare_means(comparisons = list(c("<= 65",">= 80")), label = "p.signif") +
  scale_color_manual(values = c("#A8A8A8","#011F5B")) + 
  labs(title="Cortical Composite by Chronological Age Groups",x ="Chronological Age Group", y = "Cortical Composite thickness w-score") + 
  theme_pubclean() 

p4aii <- ggplot(adni,aes(factor(GroupChronAge),MTLb)) +
  geom_violin(show.legend = FALSE) + geom_jitter(aes(col=Case),size=1, alpha=0.7, width = 0.25,show.legend = FALSE) +
  stat_summary(fun='mean',geom='crossbar',width=0.75,color='black') + 
  facet_grid(cols = vars(Case)) +
  stat_n_text() +
  stat_compare_means(comparisons = list(c("<= 65",">= 80")), label = "p.signif") +
  scale_color_manual(values = c("#A8A8A8","#990000")) + 
  labs(title="MTL Composite by Chronological Age Groups",x ="Chronological Age Group", y = "MTL Composite thickness w-score") + 
  theme_pubclean() 

p4aiii <- ggplot(adni,aes(factor(GroupChronAge),CoMeTv4b)) +
  geom_violin(show.legend = FALSE) + geom_jitter(aes(col=Case),size=1, alpha=0.7, width = 0.25,show.legend = FALSE) +
  stat_summary(fun='mean',geom='crossbar',width=0.75,color='black') + 
  facet_grid(cols = vars(Case)) +
  stat_n_text() +
  stat_compare_means(comparisons = list(c("<= 65",">= 80")), label = "p.signif") +
  scale_color_manual(values = c("#A8A8A8","#660099")) + 
  labs(title="CoMeT by Chronological Age Groups", x ="Chronological Age Group", y = "CoMeT") + 
  theme_pubclean() 

# 4B
p4bi <- ggplot(adni[adni$Case=='Symptomatic',],aes(AGE,CoMeTv4b)) +
  geom_point() + 
  geom_smooth(method = 'lm') +
  stat_cor(method='pearson') +
  labs(title="CoMeT vs Chronological Age - Symptomatic",x ="Chronological Age (years)", y = "CoMeT") +
  theme_pubclean()

p4bii <- ggplot(ad[ad$Case=='Symptomatic',],aes(AGE,CoMeTv4b)) +
  geom_point() + 
  geom_smooth(method = 'lm') +
  stat_cor(method='pearson') +
  labs(title="CoMeT vs Chronological Age - Dementia",x ="Chronological Age (years)", y = "CoMeT") +
  theme_pubclean()

# 4
grid.arrange(p4ai, p4aii, p4aiii, p4bi, p4bii, layout_matrix = matrix(c(1, 4, 2, 4, 3, 5), nrow = 2))

# 5A
p5ai <- ggplot(adni[((!is.na(adni$GroupDunedinPACE) & adni$Case=='Symptomatic')),],aes(GroupDunedinPACE,CoMeTv4b)) +
  geom_violin(aes(color=GroupDunedinPACE),show.legend = FALSE,trim = FALSE) +
  stat_summary(fun='mean',geom='crossbar',width=0.6,color='black') +
  geom_jitter(aes(color=GroupDunedinPACE),size=1, alpha=0.7, width=0.25, height=0, show.legend = FALSE) +
  stat_n_text() + 
  stat_compare_means(comparisons = list(c("Decelerated","Accelerated")), label = "p.signif") +
  labs(title="CoMeT by Biological Age Groups - Symptomatic",x ="Biological Age Group", y = "CoMeT") + 
  scale_color_manual(values = mycolor) + 
  theme_pubclean()

p5aii <- ggplot(ad[((!is.na(ad$GroupDunedinPACE) & ad$Case=='Symptomatic')),],aes(GroupDunedinPACE,CoMeTv4b)) +
  geom_violin(aes(color=GroupDunedinPACE),show.legend = FALSE,trim = FALSE) +
  stat_summary(fun='mean',geom='crossbar',width=0.75,color='black') +
  geom_jitter(aes(color=GroupDunedinPACE),size=1.1, alpha=0.8, width=0.25, height=0, show.legend = FALSE) +
  stat_n_text() +
  stat_compare_means(comparisons = list(c("Decelerated","Accelerated")), label = "p.signif") +
  labs(title="CoMeT by Biological Age Groups - Dementia",x ="Biological Age Group", y = "CoMeT") +
  scale_color_manual(values = mycolor) +
  theme_pubclean()

# 5B
p5bi <- ggplot(adni[adni$Case=='Symptomatic',],aes(AgeAccClockCombo,CoMeTv4b)) +
  geom_point() + 
  geom_smooth(method='lm') + 
  stat_cor(method='pearson') + 
  labs(title="CoMeT vs BAG - Symptomatic",x ="BAG (years)", y = "CoMeT") +
  theme_pubclean()

p5bii <- ggplot(ad,aes(AgeAccClockCombo,CoMeTv4b)) +
  geom_point() + 
  geom_smooth(method='lm') +
  stat_cor(method='pearson') +
  labs(title="CoMeT vs BAG - Dementia",x ="BAG (years)", y = "CoMeT") +
  theme_pubclean()

# 5
grid.arrange(p5ai, p5aii, p5bi, p5bii, layout_matrix = matrix(c(1,2,3,4), nrow = 2))

# S4
p9 <- ggplot(cn[cn$Case=='Cognitively Unimpaired',],aes(AGE,CoMeTv4b)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_cor(method='pearson') +
  labs(title="CoMeT vs Chronological Age - CU",x ="Chronological Age (years)", y = "CoMeT") +
  theme_pubclean()

p10 <- ggplot(cn[cn$Case=='Cognitively Unimpaired',],aes(AgeAccClockCombo,CoMeTv4b)) +
  geom_point() + 
  geom_smooth(method = 'lm') +
  stat_cor(method='pearson') +
  labs(title="CoMeT vs BAG - CU",x ="BAG (years)", y = "CoMeT") +
  theme_pubclean()

# p11 <- ggplot(cn[!is.na(cn$GroupDunedinPACE),],aes(GroupDunedinPACE,CoMeTv4b)) +
#   geom_violin(aes(color=GroupDunedinPACE),show.legend = FALSE) +
#   stat_summary(fun='mean',geom='crossbar',width=0.6,color='black') +
#   geom_jitter(aes(color=factor(GroupDunedinPACE)),size=1.2, alpha=0.8, width=0.25, height=0, show.legend = FALSE) +
#   stat_n_text() +
#   stat_compare_means(comparisons = list(c("Decelerated","Accelerated")), label = "p.signif") +
#   labs(title="CoMeT by Biological Age Groups - CU",x ="Biological Age Group", y = "CoMeT") +
#   scale_color_manual(values = c("royalblue1","grey45","firebrick1")) +
#   theme_pubclean()

grid.arrange(p9, p10, nrow = 2)

# S5
p11 <- ggplot(mci,aes(AGE,CoMeTv4b)) +
  geom_point() + 
  geom_smooth(method='lm') +
  stat_cor(method='pearson') +
  labs(title="CoMeT vs Chronological Age - MCI",x ="Chronological Age (years)", y = "CoMeT") +
  theme_pubclean()

p12 <- ggplot(mci,aes(AgeAccClockCombo,CoMeTv4b)) +
  geom_point() + 
  geom_smooth(method='lm') +
  stat_cor(method='pearson') +
  labs(title="CoMeT vs BAG - MCI",x ="BAG (years)", y = "CoMeT") +
  theme_pubclean()

# p14 <- ggplot(mci[!is.na(mci$GroupDunedinPACE),],aes(GroupDunedinPACE,CoMeTv4b)) +
#   geom_violin(aes(color=GroupDunedinPACE),show.legend = FALSE) +
#   stat_summary(fun='mean',geom='crossbar',width=0.6,color='black') +
#   geom_jitter(aes(color=factor(GroupDunedinPACE)),size=1.2, alpha=0.8, width=0.25, height=0, show.legend = FALSE) +
#   stat_n_text() +
#   stat_compare_means(comparisons = list(c("Decelerated","Accelerated")), label = "p.signif") +
#   labs(title="CoMeT by Biological Age Groups - MCI",x ="Biological Age Group", y = "CoMeT") +
#   scale_color_manual(values = c("royalblue1","grey45","firebrick1")) +
#   theme_pubclean()

grid.arrange(p11, p12, nrow = 2)

# STATS

# Two-sample t-tests
age.cortical <- wilcox.test(adni[(adni$GroupChronAge=='<= 65' & adni$Case=='Symptomatic'),'Corticalv4'], adni[(adni$GroupChronAge=='>= 80' & adni$Case=='Symptomatic'),'Corticalv4'], paired = FALSE)
age.comet <- wilcox.test(adni[(adni$GroupChronAge=='<= 65' & adni$Case=='Symptomatic'),'CoMeTv4b'], adni[(adni$GroupChronAge=='>= 80' & adni$Case=='Symptomatic'),'CoMeTv4b'], paired = FALSE)
bag.comet <- wilcox.test(adni[(adni$GroupDunedinPACE=='Decelerated' & adni$Case=='Symptomatic'),'CoMeTv4b'], adni[(adni$GroupDunedinPACE=='Accelerated' & adni$Case=='Symptomatic'),'CoMeTv4b'], paired = FALSE)
# bag.comet.cu <- wilcox.test(cn[(cn$GroupDunedinPACE=='Decelerated'),'CoMeTv4b'], cn[(cn$GroupDunedinPACE=='Accelerated'),'CoMeTv4b'], paired = FALSE)
# bag.comet.mci <- wilcox.test(mci[(mci$GroupDunedinPACE=='Decelerated' & mci$Case=='Symptomatic'),'CoMeTv4b'], mci[(mci$GroupDunedinPACE=='Accelerated' & mci$Case=='Symptomatic'),'CoMeTv4b'], paired = FALSE)
bag.comet.ad <- wilcox.test(ad[(ad$GroupDunedinPACE=='Decelerated' & ad$Case=='Symptomatic'),'CoMeTv4b'], ad[(ad$GroupDunedinPACE=='Accelerated' & ad$Case=='Symptomatic'),'CoMeTv4b'], paired = FALSE)

# Effect Sizes

# Age Cortical
U <- age.cortical$statistic # Extract U statistic
n1 <- length(adni[(adni$GroupChronAge=='<= 65' & adni$Case=='Symptomatic'),'Corticalv4']) # Sample Sizes
n2 <- length(adni[(adni$GroupChronAge=='>= 80' & adni$Case=='Symptomatic'),'Corticalv4']) 
r_rb_cortical <- (2 * U) / (n1 * n2) - 1 # Compute rank-biserial correlation

# Age CoMeT
U <- age.comet$statistic # Extract U statistic
n1 <- length(adni[(adni$GroupChronAge=='<= 65' & adni$Case=='Symptomatic'),'CoMeTv4b']) # Sample Sizes
n2 <- length(adni[(adni$GroupChronAge=='>= 80' & adni$Case=='Symptomatic'),'CoMeTv4b']) 
r_rb_comet <- (2 * U) / (n1 * n2) - 1 # Compute rank-biserial correlation

# BAG CoMeT
U <- bag.comet$statistic # Extract U statistic
n1 <- length(adni[(adni$GroupDunedinPACE=='Decelerated' & adni$Case=='Symptomatic'),'CoMeTv4b']) # Sample Sizes
n2 <- length(adni[(adni$GroupDunedinPACE=='Accelerated' & adni$Case=='Symptomatic'),'CoMeTv4b']) 
r_rb_bag <- (2 * U) / (n1 * n2) - 1 # Compute rank-biserial correlation

# BAG CoMeT AD
U <- bag.comet.ad$statistic # Extract U statistic
n1 <- length(ad[(ad$GroupDunedinPACE=='Decelerated' & ad$Case=='Symptomatic'),'CoMeTv4b']) # Sample Sizes
n2 <- length(ad[(ad$GroupDunedinPACE=='Accelerated' & ad$Case=='Symptomatic'),'CoMeTv4b']) 
r_rb_bag_ad <- (2 * U) / (n1 * n2) - 1 # Compute rank-biserial correlation

# BAG CoMeT CU
# U <- bag.comet.cu$statistic # Extract U statistic
# n1 <- length(cn[(cn$GroupDunedinPACE=='Decelerated'),'CoMeTv4b']) # Sample Sizes
# n2 <- length(cn[(cn$GroupDunedinPACE=='Accelerated'),'CoMeTv4b'])
# r_rb_bag_cu <- (2 * U) / (n1 * n2) - 1 # Compute rank-biserial correlation

# BAG CoMeT MCI
# U <- bag.comet.mci$statistic # Extract U statistic
# n1 <- length(mci[(mci$GroupDunedinPACE=='Decelerated' & mci$Case=='Symptomatic'),'CoMeTv4b']) # Sample Sizes
# n2 <- length(mci[(mci$GroupDunedinPACE=='Accelerated' & mci$Case=='Symptomatic'),'CoMeTv4b']) 
# r_rb_bag_mci <- (2 * U) / (n1 * n2) - 1 # Compute rank-biserial correlation
