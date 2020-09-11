# Get data
# 03 Jul 2020

# Location: SKM_Lab/COVID-19

# Get Duke-NUS data
source("Rcode/02_get_DukeNUS_data.R")

# Get Wuhan_PBMC data
source("Rcode/03_get_CRA002390_PBMC_data.R")

# Get Wuhan_BALF data
load(file = "Data/Wuhan_BALF.rda")

# Get Blanco_melo lung epitheliam infected with
# covid-19  data
load(file = "Data/esetNYlung_epithelium.rda") # GSE147507/
load(file = "Data/esetNYlung_biopsy.rda") # GSE147507/

# Shorten the colnames for plotting
# Lung epithelium
colnames(esetNYlung_epithelium) = gsub("Series1_NHBE_", "",colnames(esetNYlung_epithelium))
colnames(esetNYlung_epithelium) = paste0(colnames(esetNYlung_epithelium), "_lung_epithelium")
colnames(esetNYlung_epithelium) = c("Mock_1_lung_epithelium", "Mock_2_lung_epithelium", "Mock_3_lung_epithelium", "nCOV_1_lung_epithelium", "nCOV_2_lung_epithelium", "nCOV_3_lung_epithelium")
#--------------------
# Lung biopsy
colnames(esetNYlung_biopsy) = gsub("Series15_", "",colnames(esetNYlung_biopsy))
colnames(esetNYlung_biopsy) = paste0(colnames(esetNYlung_biopsy), "_biopsy")

# Wuhan PBMC
colnames(esetWuhan)  = paste0(colnames(esetWuhan), "Wu_PBMC")

# Wuhan BALF
colnames(eset.Wuhan.Balf) = paste0(colnames(eset.Wuhan.Balf), "Wu_BALF")

# DUKE Blood
colnames(esetDN) = paste0(colnames(esetDN), "_Whole.Blood")


# Merge all the data to a large list
#-----------------------------------
ncov = list()
ncov$DN.Blood = esetDN
ncov$Wuhan.pbmc = esetWuhan
ncov$Wuhan.balf = eset.Wuhan.Balf
ncov$NY.Lung.epith = esetNYlung_epithelium
ncov$NY.Lung.biops = esetNYlung_biopsy

# dat.list.cov$rv_asthama = eset
rm("dds", "eset.Wuhan.Balf", "esetDN", "esetNYlung_epithelium", "esetNYlung_biopsy", "esetWuhan")

# Get pathway gene lists
load("metadata/gbyp.rda")
names(pathways.list) = as.character(sapply(names(pathways.list), function(mystr)
  strsplit(mystr, split=":")[[1]][2]))
pathways.list = sapply(pathways.list, function(mystr)
  strsplit(mystr, split=" - ")[[1]][1])

# Get sepsis-associated module genelist
net.dat <- read.table("../nibmg2020/Manuscript.adult.sepsis/Analysis/metadata/network_input.csv", sep="\t", header=T)
net.list = with(net.dat, split(x=Egid, f=Process))
rm(net.dat)
