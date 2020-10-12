# the code to get all the 3 cases
# evaluated visually through
# a line plot in DukeNUS data
# 25 Jun 2020: Samanwoy
# 02 Jul 2020: Modified by SKM

rm(list=ls())

# Load the libraries
#---------------------
library(Biobase)
library(gplots)
library(org.Hs.eg.db)

# Get Duke-NUS data
source("Rcode/02_get_DukeNUS_data.R")
source("Rcode/funcDefs/drawGeneLinePlot.R")

# gene symbols for NETosis genes
netSyms = c("CR1","ITGAM","CCL5","CEACAM8","ITGB1", "CFH","CFB","C5","C5AR1","C3","CFP","MPO","ELANE","CTSG","HMGB1","AGER","TLR2","TLR4","H4C1","TF","TFPI","F2","FGB","PLG")
egSyms = links(org.Hs.egSYMBOL2EG[netSyms])
egs = egSyms[,"gene_id"]
egs = intersect(egs, featureNames(esetDN))
for(eg in egs) {
  drawGeneLinePlot(eg);scan()}
