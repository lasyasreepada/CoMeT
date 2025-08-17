library(dplyr)
library(sesame)
library(knowYourCG)
library(SummarizedExperiment)
library(GenomicRanges)

sesameDataCache()
setwd('/Users/lasyasreepada/Projects/CoMeT/')

smry <- readRDS('~/Projects/CoMeT/data/ADNI/20250504_COMET_SMRY_NOAGE.rds')

qry_age <- smry %>% arrange(Pval_CoMeTv4b) %>% dplyr::select(Probe_ID,Est_CoMeTv4b,Pval_CoMeTv4b,FDR) %>%
  filter(Pval_CoMeTv4b < 1e-5)

qry_noage <- smry %>% arrange(Pval_CoMeTv4b) %>% dplyr::select(Probe_ID,Est_CoMeTv4b,Pval_CoMeTv4b,FDR) %>%
  filter(Pval_CoMeTv4b < 1e-5)

source("analysis/linkgenes.R")

cg <- qry_age$Probe_ID
bp <- 5000

cg_gene_link <- linkGenes(
  qry=cg,
  platform="EPIC",
  distance=bp
)

library(enrichR)

websiteLive <- getOption("enrichR.live")

if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes
}

dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", 
         "GO_Biological_Process_2023")

input <- cg_gene_link$Gene

if (websiteLive) {
  enriched <- enrichr(input, dbs)
}

biol <- enriched[["GO_Biological_Process_2023"]]
mol <- enriched[["GO_Molecular_Function_2023"]]
cell <- enriched[["GO_Cellular_Component_2023"]]

if (websiteLive) {
  plotEnrich(biol, showTerms = 30, numChar = 80, 
             y = "Count", orderBy = "P.value")
}

if (websiteLive) {
  plotEnrich(mol, showTerms = 15, numChar = 80, 
             y = "Count", orderBy = "P.value")
}

if (websiteLive) {
  plotEnrich(cell, showTerms = 3, numChar = 60, 
             y = "Count", orderBy = "P.value")
}

# CSV
enriched_df <- rbind(biol,mol, cell)

