library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(EnvStats)

setwd('/Users/lasyasreepada/Projects/CoMeT/')

# READ DATA
adni <- read.csv('data/ADNI/MRI_DNA_COG_XSH.csv')
adni <- subset(adni, select = -c(X.1,X))
adni <- adni[adni$Case==0,]

# PIVOT VOLUME
regions <- c('ST29SV','ST29SV_H')
data <- adni[,c('RID','FLDSTRENG',regions)]
data <- data %>% pivot_longer(cols=regions,
                                 names_to='name',
                                 values_to='atrophy')

# HARMONIZATION
data$harmonized <- ""
data$harmonized[grep("_H", data$name)] <- "Harmonized"
data[data$harmonized=="","harmonized"] <- 'Raw'

# RENAME COLS
data$name <- sapply(strsplit(data$name, "\\."), `[`, 1)

# STATS
res.raw <- t.test(data[(data$FLDSTRENG==1.5 & data$harmonized=='Raw'),'atrophy'], data[(data$FLDSTRENG==3.0 & data$harmonized=='Raw'),'atrophy'], paired = FALSE)
res.harmonized <- t.test(data[(data$FLDSTRENG==1.5 & data$harmonized=='Harmonized'),'atrophy'], data[(data$FLDSTRENG==3.0 & data$harmonized=='Harmonized'),'atrophy'], paired = FALSE)

# PLOT
p1 <- ggplot(data[data$harmonized=='Raw',],aes(factor(FLDSTRENG),atrophy)) +
  geom_violin(show.legend = FALSE) + 
  stat_summary(fun='mean',geom='crossbar',width=0.5,color='blue') + 
  stat_n_text() + stat_compare_means(method = 't') +
  labs(title = 'Left Hippocampus Volume - Raw', x = 'Field Strength', y = 'Volume (mm3)') +
  geom_jitter(size=1.1, alpha=0.8, width=0.3, height=0, show.legend = FALSE) +
  theme_pubclean()

p2 <- ggplot(data[data$harmonized=='Harmonized',],aes(factor(FLDSTRENG),atrophy)) +
  geom_violin(show.legend = FALSE) + 
  stat_summary(fun='mean',geom='crossbar',width=0.5,color='blue') + 
  stat_n_text() + stat_compare_means(method = 't') +
  labs(title = 'Left Hippocampus Volume - Harmonized', x = 'Field Strength', y = 'Volume (mm3)') +
  geom_jitter(size=1.1, alpha=0.8, width=0.3, height=0, show.legend = FALSE) +
  theme_pubclean()


# PIVOT THICKNESS
regions <- c('ST52TA','ST52TA_H')
data <- adni[,c('RID','FLDSTRENG',regions)]
data <- data %>% pivot_longer(cols=regions,
                              names_to='name',
                              values_to='atrophy')

# HARMONIZATION
data$harmonized <- ""
data$harmonized[grep("_H", data$name)] <- "Harmonized"
data[data$harmonized=="","harmonized"] <- 'Raw'

# RENAME COLS
data$name <- sapply(strsplit(data$name, "\\."), `[`, 1)

# STATS
res.raw <- t.test(data[(data$FLDSTRENG==1.5 & data$harmonized=='Raw'),'atrophy'], data[(data$FLDSTRENG==3.0 & data$harmonized=='Raw'),'atrophy'], paired = FALSE)
res.harmonized <- t.test(data[(data$FLDSTRENG==1.5 & data$harmonized=='Harmonized'),'atrophy'], data[(data$FLDSTRENG==3.0 & data$harmonized=='Harmonized'),'atrophy'], paired = FALSE)


p3 <- ggplot(data[data$harmonized=='Raw',],aes(factor(FLDSTRENG),atrophy)) +
  geom_violin(show.legend = FALSE) + 
  stat_summary(fun='mean',geom='crossbar',width=0.5,color='blue') + 
  stat_n_text() + stat_compare_means(method = 't') +
  labs(title = 'Left Precuneus Thickness - Raw', x = 'Field Strength', y = 'Thickness (mm)') +
  geom_jitter(size=1.1, alpha=0.8, width=0.3, height=0, show.legend = FALSE) +
  theme_pubclean()

p4 <- ggplot(data[data$harmonized=='Harmonized',],aes(factor(FLDSTRENG),atrophy)) +
  geom_violin(show.legend = FALSE) + 
  stat_summary(fun='mean',geom='crossbar',width=0.5,color='blue') + 
  stat_n_text() + stat_compare_means(method = 't') +
  labs(title = 'Left Precuneus Thickness - Harmonized', x = 'Field Strength', y = 'Thickness (mm)') +
  geom_jitter(size=1.1, alpha=0.8, width=0.3, height=0, show.legend = FALSE) +
  theme_pubclean()

grid.arrange(p1, p2, p3, p4, nrow=2)



