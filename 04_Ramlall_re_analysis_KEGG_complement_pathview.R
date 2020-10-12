# code to reanalysis Ramllal et. al. Nature Medicine
# Link: https://www.nature.com/articles/s41591-020-1021-2
# Supplementary data used
#dat = read.xlsx("C:/Users/Samanwoy/Documents/Ramllal_et_al/SuppTable_3_4_DEGandPathways_preprint.xlsx", sheet = 4)

# Samanwoy: "Sat Aug 29 05:15:21 2020"
#------------------------------------
# Supplementary Figure 1: Pathview of
# KEGG complement Coagulation pathway
# from Ramllal et. al.2020
#------------------------------------
library(openxlsx)
library(Biobase)
library(org.Hs.eg.db)

# read the data
#---------------
load("Data/Ramlal.lfcdata.rda")

# check KEGG complement coagulation pathway genes
#------------------------------------------------
load("metadata/gbyp.rda")
comp.gns = genes.by.pathway[["hsa04610"]]
comp.gns.gsyms = as.character(unlist(mget(comp.gns,
                                          org.Hs.egSYMBOL)))

comp.lfc = dat$log2FoldChange[which(dat$Gene%in%comp.gns.gsyms)]
names(comp.lfc) = dat$Gene[which(dat$Gene%in%comp.gns.gsyms)]

# Changing the gene symbols to ENTREZID
names(comp.lfc) = as.character(unlist(mget(names(comp.lfc), org.Hs.egSYMBOL2EG)))

if(file.exists("Results/Ramlall_et_al_results/Ramllal_et_al_KEGG_hsa04610.pathview.png")){
  cat("No need to run Pathview Image file exists")
}else{
  cat("Running Pathview")
  library(pathview)
  pathview(comp.lfc, pathway.id = "hsa04610")
}

# Change the downregulated genes value to 0 for visualisation with gray in pathview
#-------------------------------------------------------------------
comp.lfc.up = comp.lfc
comp.lfc.up[comp.lfc.up<0] = 0
library(pathview)
pathview(comp.lfc.up, pathway.id = "hsa04610", out.suffix = "up_regulated_gns_only")

