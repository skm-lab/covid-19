# This code compares genes common between two COVID datasets
# Then draws comparitive boxplots for theose common genes
# samanwoy "Wed Jun 10 14:43:03 2020"
####################################
rm(list = ls())
graphics.off()
gc()
# Load the three biological module genes
#----------------------------------------
net.dat <- read.table("C:/Users/Samanwoy/Dropbox/SKM_LAB/nibmg2020/Manuscript.adult.sepsis/Analysis/metadata/network_input.csv", sep="\t", header=T)
# net.dat <- read.table("C:/Users/Saroj/Dropbox/Work/SKM_LAB/nibmg2020/Manuscript.adult.sepsis/Analysis/metadata/network_input.csv", sep="\t", header=T)
net.list = with(net.dat, split(x=Egid, f=Process))

# pathid = "coagulation"
# pathid = "immunosuppression"
 pathid = "inflammation"
# Load the libraries
#---------------------
library(Biobase)
library(gplots)
#library(ggplot2)
#library(ggpubr)
library(org.Hs.eg.db)

# Get Duke-NUS data
source("Rcode/02_get_DukeNUS_data.R")

# Get Wuhan_PBMC data
source("Rcode/03_get_CRA002390_PBMC_data.R")

# Get Wuhan_BALF data
load(file = "Wuhan_BALF/Data/Wuhan_BALF.rda")

# Get Blanco_melo lung epitheliam infected with
# covid-19  data
load(file = "GSE147507/Data/esetNYlung_epithelium.rda")
load(file = "GSE147507/Data/esetNYlung_biopsy.rda")


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


# Get RV infections in asthmatics data
# load(file = "GSE149273/GSE149273.rda")

# Merge all the data to a large list
#-----------------------------------
dat.list.cov = list()
dat.list.cov$ncov_DN = esetDN
dat.list.cov$ncov_Wuhan_pbmc = esetWuhan
dat.list.cov$ncov_Wuhan_balf = eset.Wuhan.Balf
dat.list.cov$ncov_Lung_epithelium = esetNYlung_epithelium
dat.list.cov$ncov_Lung_biopsy = esetNYlung_biopsy

# dat.list.cov$rv_asthama = eset

rm("dds", "eset.Wuhan.Balf", "esetDN", "esetNYlung_epithelium", "esetNYlung_biopsy", "esetWuhan", "net.dat")

# Make the data set changed such a
# way only common genes exist
#-----------------------------------
# Common genes between two datsets
cmn.Cov.Gns.1 = intersect(featureNames(dat.list.cov$ncov_DN), featureNames(dat.list.cov$ncov_Wuhan_pbmc))

cmn.Cov.Gns.2 = intersect(featureNames(dat.list.cov$ncov_Wuhan_balf), featureNames(dat.list.cov$ncov_Lung_epithelium))

cmn.Cov.Gns.3 = intersect(cmn.Cov.Gns.1, cmn.Cov.Gns.2)

cmn.Cov.Gns = intersect(cmn.Cov.Gns.3,
                        featureNames(dat.list.cov$ncov_Lung_biopsy))

rm(cmn.Cov.Gns.1, cmn.Cov.Gns.2, cmn.Cov.Gns.3)

# get the genes common in both datasets and Module
#--------------------------------------------------
path.cmn.gns = intersect(cmn.Cov.Gns,
                         as.character(net.list[[pathid]]))
dat.list.cov.cmngns = sapply(dat.list.cov,
                             function(eset){eset[featureNames(eset)%in%path.cmn.gns,]})
#--------------------------------------
# add suffix for cohort on each samples
#--------------------------------------
#-----------------------------------
# create a function
# Given an expression set onject
# it subtracts median control gene
# expression for each gene
# and subtracts that value from all cases
# returns a matrix of only genes by cases
#------------------------------------------
getCtrlNormalised = function(eset =eset){
  ctrl.mdn = rowMedians(exprs(eset[,eset$Group=="Zcontrol"]))
  names(ctrl.mdn) = path.cmn.gns
  case.mat = exprs(eset[,eset$Group!="Zcontrol"]) - ctrl.mdn
  # return  control normalised case matrix
  return(case.mat)
}

plot.dat = do.call(cbind, sapply(dat.list.cov.cmngns,
                                 function(eset){getCtrlNormalised(eset)}))

# Draw Sample level boxplot of all COVID data
###############################################
ttl.str = paste0("Module: ", pathid)

colvec =   c(rep("pink1",9),
             rep("orchid1",9),rep("pink1",4), rep("red1", 3), rep("deepskyblue", 2), rep("purple", 3), rep("cyan3", 2))
x11()
par(mar=c(12,4,1,1),
    cex.main = 1.2)
boxplot(plot.dat,
        outline = F,
        las=2,
        main = ttl.str,
        col = colvec,
        ylab = "Log fold change"
)
# add Legend
legend("bottomleft", c("DN_Case-1", "DN_Case-2", "DN_Case-3", "Wuhan_PBMC", "Wuhan_BALF", "NY_Lung_Epithelium", "NY_Lung_Biopsy"), fill = c("pink1", "orchid1", "pink1", "red1", "deepskyblue", "purple", "cyan3"), bty = "n")
abline(h=0, lty = 2, lwd = 2 , col= "grey40")

#-----------------------
# Draw pathway heatmap
#-----------------------
egs = path.cmn.gns
gnames = AnnotationDbi::select(org.Hs.eg.db,
                               keys=egs,
                               column=c("SYMBOL","GENENAME"))
labRowStr = with(gnames,
                 paste0(GENENAME, " | ", ENTREZID))

par(cex.main=0.65)
x11()
heatmap.2(plot.dat,
          Rowv=,
          Colv=FALSE,
          dendrogram = "row",
          scale="row",
          trace="none",
          density.info="none",
          keysize = 0.91,
          key.title = " ",
          key.xlab = "Log fold change",
          key.ylab = NULL,
          col = bluered(256),
          margins = c(12,25),
          main = ttl.str,
          labRow = labRowStr,
          ColSideColors = colvec,
          sepwidth = c(0.05,0.05),
          colsep = c(9, 18, 22, 25, 27),
          sepcolor = "white",
          rowsep = 1:nrow(plot.dat))
legend(y=1.17, x=.75, xpd=TRUE,
       legend = c("DN_Case-1", "DN_Case-2", "DN_Case-3", "Wuhan_PBMC", "Wuhan_BALF", "NY_Lung_Epithelium", "NY_Lung_Biopsy"),
       fill = c("pink1", "orchid1", "pink1", "red1", "deepskyblue", "purple", "cyan3"),
       #lty= 1,
       #lwd = 5,
       cex=.87, bty = "n"
)
