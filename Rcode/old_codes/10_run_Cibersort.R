# run CIBERSORT
# Preapare the gene expression
# matrix which acts as
# an input file to CIBERSORTx
# author: Samanwoy
# "Mon Jul 06 11:31:33 2020"
############################
# Location: SKM_Lab/COVID-19
rm(list = ls())
graphics.off()
gc()
options(digits=2)

# Load the libraries
#---------------------
library(Biobase)
library(gplots)
library(org.Hs.eg.db)

# function definitions
#source("Rcode/funcDefs/plotFuncs.R")

# Load data and metadata
source("Rcode/getData.R")

# start with duke cohort:
# then wuhan PBMC cohort:
# then wuhan BALF cohort:
# then Ny lung Epithelium:

# get gene expression matrix #UN-Logged
MIT.ethmoid.sinus.ncov_gexp = 2^exprs(esetMIT)
# create a directory for CIBERSORT analysis
dir.create("cibersort_COVID-19")

# Change the entrez gene ids to gene symbol
MIT.ethmoid.sinus.gsyms =as.character(unlist(mget(rownames(MIT.ethmoid.sinus.ncov_gexp), org.Hs.egSYMBOL, ifnotfound = NA)))

#sel = which(is.na(wuBALF.gsyms))
#wuhan_BALF_ncov_gexp.new = wuhan_BALF_ncov_gexp[-c(sel),]

# change the rowname of gene expression matrix
rownames(MIT.ethmoid.sinus.ncov_gexp) = as.character(unlist(mget(rownames(MIT.ethmoid.sinus.ncov_gexp), org.Hs.egSYMBOL)))

# write the data into a tsv file
write.table(MIT.ethmoid.sinus.ncov_gexp, "cibersort_COVID-19/MIT.ethmoid.sinus.ncov_gexp_unlogged.tsv", sep= "\t", row.names = T, col.names = NA)


