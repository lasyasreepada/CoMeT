library(dplyr)
library(ggplot2)
library(ggpubr)
library(EnvStats)
library(GGally)
library(gridExtra)
library(patchwork)

setwd('/Users/lasyasreepada/Projects/CoMeT/')

# Get ADNI Cases
adni <- read.csv('data/ADNI/MRI_DNA_COG_XSH_COMET_Grouped.csv')

clockCols <- c('AGE_SAMPLE_nearest',
               'HorvathS2013_nearest',
               'HannumG2013_nearest',
               'ShirebyG2020_nearest',
               'AgeClockCombo',
               'HorvathS2013_res_nearest',
               'HannumG2013_res_nearest',
               'ShirebyG2020_res_nearest',
               'AgeAccClockCombo')

adniClocks <- subset(adni,select=clockCols)

my_func <- function(data, mapping, ...) {
  ggplot(data, mapping) + 
    geom_point(size = 0.5) +
    geom_smooth(formula = y~x, color = "blue") + 
    theme_pubclean()
}

my_fn <- function(data, mapping, method="p", use="pairwise", ...){
  
  # grab data
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  # calculate correlation
  corr <- cor(x, y, method=method, use=use)
  
  # calculate colour based on correlation value
  # Here I have set a correlation of minus one to blue, 
  # zero to white, and one to red 
  # Change this to suit: possibly extend to add as an argument of `my_fn`
  colFn <- colorRampPalette(c("blue", "white", "red"), interpolate ='spline')
  fill <- colFn(100)[findInterval(corr, seq(-1, 1, length=100))]
  
  ggally_cor(data = data, mapping = mapping, ...) + 
    theme_pubr() +
    theme(panel.background = element_rect(fill=fill))
}

p1 <- ggpairs(
  adniClocks,
  columns = c(1,6:9),
  columnLabels = c('Age','Horvath BAG','Hannum BAG','Shireby BAG','BAG'),
  lower = list(continuous = my_func,displayGrid = FALSE),
  upper = list(continuous = my_fn,displayGrid = FALSE),
  title = "Correlation Matrix of Biological Age Gaps and Chronological Age",
) 

p2 <- ggpairs(
  adniClocks,
  columns = c(1:5),
  columnLabels = c('Age','Horvath 2013','Hannum 2013','Shireby 2020','Clock Average'),
  lower = list(continuous = my_func,displayGrid = FALSE),
  upper = list(continuous = my_fn,displayGrid = FALSE),
  title = "Correlation Matrix of Epigenetic Clock Ages and Chronological Age",
)

p1
p2
