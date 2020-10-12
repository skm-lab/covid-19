# This code compares genes common between two COVID datasets
# Then draws comparitive boxplots for theose common genes
# samanwoy "Wed Jun 10 14:43:03 2020"
####################################
rm(list = ls())
graphics.off()
gc()
# Load the three biological module genes
#----------------------------------------
net.dat <- read.table("D:/Dropbox/SKM_LAB/nibmg2020/Manuscript.adult.sepsis/Analysis/metadata/network_input.csv", sep="\t", header=T)
#net.dat <- read.table("C:/Users/Saroj/Dropbox/Work/SKM_LAB/nibmg2020/Manuscript.adult.sepsis/Analysis/metadata/network_input.csv", sep="\t", header=T)
net.list = with(net.dat, split(x=Egid, f=Process))

# pathid = "coagulation"
# pathid = "immunosuppression"
# pathid = "inflammation"
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

# Get Blanco_melo lung epetheliam infected with
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

# DUKE PBMC
colnames(esetDN) = paste0(colnames(esetDN), "_PBMC")


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

#-----------------------------------
# create a function
# Given an expression set onject
# it subtracts median control gene
# expression for each gene
# and subtracts that value from all cases
# returns a matrix of only genes by cases
#------------------------------------------
getCtrlNormalised = function(eset = esetDN, pathid = pathid){
  path.cmn.gns = intersect(featureNames(eset),
                           net.list[[pathid]])
  eset.nw = eset[which(featureNames(eset)%in%path.cmn.gns),]
  ctrl.mdn = rowMedians(exprs(eset.nw[,eset.nw$Group=="Zcontrol"]))
  names(ctrl.mdn) = path.cmn.gns
  case.mat = exprs(eset.nw[,eset.nw$Group!="Zcontrol"]) - ctrl.mdn
  # return  control normalised case matrix
  return(case.mat)
}
# apply the function to return
# list of LFC matrices
#----------------------
plot.dat.list = sapply(dat.list.cov,
                                 function(eset){getCtrlNormalised(eset = eset, pathid = pathid)})
#---------------------------------------------
# Draw Sample level boxplot of all COVID data
#---------------------------------------------
drawModuleBoxplot = function(study.name = "ncov_DN" ){
  plot.dat = plot.dat.list[[study.name]] # createthe plot data
  # Create the title string
  ttl.str = paste0("Module: ", pathid, ": ", names(plot.dat.list[study.name]))
  # Create color vector
  #------------------------
  library(RColorBrewer)
  if(names(plot.dat.list[study.name]) == "ncov_DN"){
    colvec =   c(rep("#9E0142",9),
                 rep("#FDAE61",9),rep("#D53E4F",4))
  }else{colvec = brewer.pal(n = ncol(plot.dat), name = "Spectral")
  }
  par(mar=c(12,4,1,1),
      cex.main = 1.2)
  boxplot(plot.dat,
          outline = F,
          las=2,
          main = ttl.str,
          col = colvec,
          ylab = "Log fold change")
  abline(h=0, lty = 2, lwd = 2 , col= "grey40")
  # add Legend
  #-----------
  if(names(plot.dat.list[study.name]) == "ncov_DN"){
    legend("bottomleft",
           legend = c("Case_1", "Case_2", "Case_3"),
           fill = c("#9E0142", "#FDAE61", "#D53E4F"),
           bty = "n")
  }else{legend("bottomleft",
               legend = colnames(plot.dat),
               fill = colvec,
               #lty= 1,
               #lwd = 5,
               cex=.7, bty = "n")}
  }
# Save the result in a pdf file
#-------------------------------
# pdf(file = paste0("Module_", pathid,"_Boxplot_Covid-19.pdf"), onefile = T)
# sapply(names(plot.dat.list), function(x){
#   drawModuleBoxplot(study.name = x)
# })
# dev.off()
rm("drawModuleBoxplot", "getCtrlNormalised", "net.list", "dat.list.cov")



#-----------------------
# Draw pathway heatmap
#-----------------------
drawModuleHeatmap = function(study.name = "ncov_DN"){
  plot.dat = plot.dat.list[[study.name]] # createthe plot data
  # Create the title string
  ttl.str = paste0("Module: ", pathid, names(plot.dat.list[study.name]))
  # Create color vector
  #------------------------
  library(RColorBrewer)
  if(names(plot.dat.list[study.name]) == "ncov_DN"){
    colvec =   c(rep("#9E0142",9),
                 rep("#FDAE61",9),rep("#D53E4F",4))
  }else(colvec = rainbow(ncol(plot.dat))
  )
  egs = rownames(plot.dat)
  gnames = AnnotationDbi::select(org.Hs.eg.db,
                                 keys=egs,
                                 column=c("SYMBOL","GENENAME"))
  labRowStr = with(gnames,
                   paste0(GENENAME, " | ", ENTREZID))

  par(cex.main=0.85)
  #x11()
  heatmap.2(plot.dat,
            Rowv=FALSE,
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
            margins = c(9,25),
            main = ttl.str,
            labRow = labRowStr,
            cexCol = 0.7,
            ColSideColors = colvec,
            sepwidth = c(0.05,0.05),
            #colsep = c(9, 18, 22, 25, 27),
            sepcolor = "white",
            rowsep = 1:nrow(plot.dat))
  if(names(plot.dat.list[study.name]) == "ncov_DN"){
    legend(y=1.17, x=.75, xpd=TRUE,
           legend = c("Case_1", "Case_2", "Case_3"),
           fill = c("#9E0142", "#FDAE61", "#D53E4F"),
           bty = "n")
  }else{legend(y=1.1, x=.65, xpd=TRUE,
               legend = colnames(plot.dat),
               fill = colvec,
               cex=.87,
               bty = "n")
  }
  }
# Save the result in a pdf file
#-------------------------------
pdf(file = paste0("Module_", pathid,"_Heatmap_Covid-19.pdf"), onefile = T, width=9, height= 7)
sapply(names(plot.dat.list), function(x){
  drawModuleHeatmap(study.name = x)
})
dev.off()
