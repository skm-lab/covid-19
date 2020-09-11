# checking corrplot for
# cytokine and netosis
# genes in DUKe NUS data
# Samanwoy:
library(Biobase)
library(gplots)
#library(ggplot2)
library(org.Hs.eg.db)

# Get Duke-NUS data
source("Rcode/02_get_DukeNUS_data.R")

#### Load the Cytokine and NETosis genes ####
#---------------------------------------------------
egs = read.csv(file="metadata/input_file_for_netosis_c_storm_gns.csv")
egs.n = unique(as.character(egs[egs$class=="netosis","egids"]))

egs.n = intersect(egs.n, featureNames(esetDN))

egs.c = unique(as.character(egs[egs$class=="cStorm","egids"]))
egs.c = intersect(egs.c, featureNames(esetDN))
rm(egs)


getCorrelogram = function(caseid = "Case1"){

  #### create case1 datatframe ####
  #--------------------------------
  datf.case1 = exprs(esetDN)[c(egs.c, egs.n), esetDN$ID== caseid]
  rownames(datf.case1) = as.character(links(org.Hs.egSYMBOL[rownames(datf.case1)])[,"symbol"])

  #### Convert egids to gene symbols#####
  #---------------------------------------------------
  gsyms.c = as.character(links(org.Hs.egSYMBOL[egs.c])[,"symbol"])
  gsyms.n = as.character(links(org.Hs.egSYMBOL[egs.n])[,"symbol"])

  # add a token of class-identity to the gene symbols
  #---------------------------------------------------
  rownames(datf.case1)[rownames(datf.case1)%in%gsyms.c] = paste0("Cytokines_", intersect(rownames(datf.case1),gsyms.c))
  rownames(datf.case1)[rownames(datf.case1)%in%gsyms.n] = paste0(intersect(rownames(datf.case1),gsyms.n), "__NETosis")

  # create the gene by gene correlation matrix
  #---------------------------------------------------
  cor.dat = cor(t(datf.case1)) # for case-1 (Severe)
  rows.to.keep = grep("Cytokines", rownames(cor.dat))
  cols.to.keep = grep("NETosis", colnames(cor.dat))
  cor.plotdat = cor.dat[rows.to.keep, cols.to.keep]
  rm(cor.dat, rows.to.keep, cols.to.keep)


  # Draw the correlogram
  library(RColorBrewer)
  #library(corrplot)
  #par(ps=12)
  #corrplot(corr = cor.plotdat, order = "hclust", method = "square", #col=brewer.pal(n=8, name="RdBu"), tl.cex = 0.65, tl.col = "black")
  library(gplots)

  graphics.off()
  heatmap.2(t(cor.plotdat),
            cellnote= formatC(t(cor.plotdat), digits = 1),
            notecex=1.0,
            notecol="yellow2",
            na.color=par("bg"),
            dendrogram = "none",
            trace="none",
            density.info="none",
            key = TRUE,
            keysize = 1.91,
            key.title = "Correlation",
            key.xlab = "Cor-Coef",
            key.ylab = NULL,
            col = redblue(256),
            margins = c(8, 9),
            main = paste0("Correlogram: ", caseid),
            cexCol = 1,
            sepwidth = c(0.05,0.05),
            colsep = 1:ncol(t(cor.plotdat)),
            sepcolor = "white",
            rowsep = 1:nrow(t(cor.plotdat))
  )

return(cor.plotdat)
}





