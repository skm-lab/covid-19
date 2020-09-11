# 07_plot_a_gene.R
# Saroj Kant Mohapatra
# 01-07-2020
####################################

library(Biobase)
library(gplots)
#library(ggplot2)
#library(ggpubr)
library(org.Hs.eg.db)

##################################
# Get the data files
##################################
if(1) {
  rm(list=ls())
  # Get Duke-NUS data
  source("Rcode/02_get_DukeNUS_data.R")

  # Get Wuhan_PBMC data
  source("Rcode/03_get_CRA002390_PBMC_data.R")

  # Get Wuhan_BALF data
  load(file = "Wuhan_BALF/Data/Wuhan_BALF.rda")

  # Get Blanco_melo lung epitheliam infected with
  # covid-19  data
  load(file = "GSE147507/Data/esetNYlung.rda")

  # dat.list.cov$rv_asthama = eset

  #rm("dds", "eset.Wuhan.Balf", "esetDN", "esetNYlung", "esetNYlung", "esetWuhan")
}

# create a function
# Given an expression set onject
# it subtracts median control gene
# expression for each gene
# and subtracts that value from all cases
# returns a matrix of only genes by cases
#------------------------------------------
getCtrlNormalised = function(eset =eset) {
  ctrl.mdn = rowMedians(exprs(eset[,eset$Group=="Zcontrol"]))
  case.mat = exprs(eset[,eset$Group!="Zcontrol"]) - ctrl.mdn
  # return  control normalised case matrix
  return(case.mat)
}

# Color scheme for Duke-NUS data
mycols=c(rep("red",9),rep("blue",13))

# RANTES | Entrez ID: 6352
eg = "6352"
gname = paste(get(eg, org.Hs.egGENENAME), " | ", get(eg, org.Hs.egSYMBOL), " | ", "RANTES")
gxp=getCtrlNormalised(esetDN[eg,])
x = (gxp[1,])
par(mar=c(6,6,4,2))
barplot(x, las=2, col=mycols,
        ylab="Gene experssion level\nNormalized to control",
        main=gname)
legend(x="topright",legend = c("Severe Illness","Moderate Illness"),
       bty="n", fill=c("Red","Blue"))
rm(eg, gxp, x)


# CD18 | Entrez ID: 3689
eg = "3689"
gname = paste(get(eg, org.Hs.egGENENAME), " | ", get(eg, org.Hs.egSYMBOL), " | ", "CD18")
gxp=getCtrlNormalised(esetDN[eg,])
x = (gxp[1,])
par(mar=c(6,6,4,2))
barplot(x, las=2, col=mycols,
        ylab="Gene experssion level\nNormalized to control",
        main=gname)
legend(x="topright",legend = c("Severe Illness","Moderate Illness"),
       bty="n", fill=c("Red","Blue"))
rm(eg, gxp, x)


# ITGAM integrin subunit alpha M [ Homo sapiens (human) ]
# CD11b | Entrez ID: 3684
eg = "3684"
gname = paste(get(eg, org.Hs.egGENENAME), " | ", get(eg, org.Hs.egSYMBOL), " | ", "CD11b")
gxp=getCtrlNormalised(esetDN[eg,])
x = (gxp[1,])
par(mar=c(6,6,4,2))
barplot(x, las=2, col=mycols,
        ylab="Gene experssion level\nNormalized to control",
        main=gname)
legend(x="topright",legend = c("Severe Illness","Moderate Illness"),
       bty="n", fill=c("Red","Blue"))
rm(eg, gxp, x)



# CD66b | Entrez ID: 1088
eg = "1088"
gname = paste(get(eg, org.Hs.egGENENAME), " | ", get(eg, org.Hs.egSYMBOL), " | ", "CD66b")
gxp=getCtrlNormalised(esetDN[eg,])
x = (gxp[1,])
par(mar=c(6,6,4,2))
barplot(x, las=2, col=mycols,
        ylab="Gene experssion level\nNormalized to control",
        main=gname)
legend(x="topright",legend = c("Severe Illness","Moderate Illness"),
       bty="n", fill=c("Red","Blue"))
rm(eg, gxp, x)
##################################################
#
# Get the NETosis module from Cynthiya_nature_2018
#
###################################################
gn.list = read.table("Pubs/pubmed/NET/NETosis_module_cynthiya_et_al_2018.txt", sep="\t", header = F)

netosis.egids = as.character(unlist(mget(gn.list[,1], org.Hs.egSYMBOL2EG, ifnotfound = NA)))
rm(gn.list)

# create a function for barplot
drawLFCbarPlot = function(eg = "1088", mycols = c(rep("red",9),rep("blue",13))){
  gname = paste(get(eg, org.Hs.egGENENAME),
                "\n | ",
                get(eg, org.Hs.egSYMBOL),
                " | ",
                "CD66b")
  gxp = getCtrlNormalised(esetDN[eg,])
  x = (gxp[1,])
  par(mar = c(6,6,4,2))
  barplot(x,
          las = 2,
          col = mycols,
          ylab = "Gene experssion level\nNormalized to control",
          main = gname)
  legend(x ="topright",
         legend = c("Severe Illness","Moderate Illness"),
         bty = "n",
         fill = c("Red","Blue"))
}

# Duke data
#--------------
selgns = intersect(netosis.egids, featureNames(esetDN))

pdf("Pubs/pubmed/NET/NETosis_gene_Barplots_DUKE_Cohort.pdf", onefile= T)
for(egid in selgns){
  drawLFCbarPlot(eg = egid, mycols = c(rep("red",9),rep("blue",13)))
#  scan()#
}
dev.off()
# Wuhan PBMC
#--------------
selgns = intersect(netosis.egids, featureNames(esetWuhan))

for(egid in selgns){
  drawLFCbarPlot(eg = egid, mycols = rep("red1", ) )
  scan()
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
