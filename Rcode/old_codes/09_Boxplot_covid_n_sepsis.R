# This code tests the global data
# for assessment of significance
# of modules in global data
# for both Sepsis and COVID-19
# Samanwoy
# Thu Jul 02 17:51:15 2020
########################################################################
# Boxplot
########################################################################
rm(list = ls())
graphics.off()
# get the data for validation
library(ggplot2)
library(cowplot)
library(ssgeosurv)
data(ss.list)

# sorce the function definition file
source("Rcode/prelim.R")
# select the studies with Control samples
# "GSE26378" "GSE26440" "GSE4607"
# "GSE9692"  "GSE48080" "GSE54514"
# "GSE95233"
# sel.studies = c(1:4,6)
# Load the three biological module genes
#----------------------------------------
net.dat <- read.table("C:/Users/Samanwoy/Dropbox/SKM_LAB/nibmg2020/Manuscript.adult.sepsis/Analysis/metadata/network_input.csv", sep="\t", header=T)
net.list = with(net.dat, split(x=Egid, f=Process))
# Drawing monotonic boxplots
comb.plot.lst = list()
mod.nms = c("Immunosuppression", "Coagulation", "Inflammation")
for (j in 1:length(ss.list)) {
  eset.s = ss.list[[j]]
  #eset.s = ss.list[[j]][,eset.s$Outcome%in%c("survivor", "nonsurvivor")]


  plot.lst = list()
  for(i in 1:3) {
    geneids = net.list[[mod.nms[i]]]
    dat = colMeans(exprs(eset.s[featureNames(eset.s)%in%geneids,]))
    fac = factor(as.character(eset.s$Outcome), levels = c("control", "survivor", "nonsurvivor"))
    plot.dat = data.frame("name" =fac, "value" = dat)
    titlestr = paste0(mod.nms[i], "--", names(ss.list)[j])

    # drawing the plot
    plot.lst[[i]] = ggplot(plot.dat, aes(x=name, y=value, fill=name)) +
      geom_boxplot(show.legend = T, width=0.55, outlier.shape = NA)+ theme(legend.title=element_blank())+ xlab("")+ylab("Gene expression (log2)")+ scale_fill_manual(values = c("#619CFF", "#00BA38", "#F8766D"))+theme(aspect.ratio=1)+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+theme(axis.text.y = element_text(face="bold", color="black", size=15))
  }
  comb.plot.lst[[j]] = plot.lst
}

# p = ggpubr::ggarrange(
#   comb.plot.lst[[2]][[1]], comb.plot.lst[[2]][[2]], comb.plot.lst[[2]][[3]],
#   comb.plot.lst[[3]][[1]], comb.plot.lst[[3]][[2]],comb.plot.lst[[3]][[3]], common.legend = T, legend = "bottom",align = "v", labels = c("Immunosuppression", "    Coagulation", "    Inflammation"), hjust = -0.3)

# ggsave("Supplementary_Figure_S3.jpeg", p, height = 7, width = 9)


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

