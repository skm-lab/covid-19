#####################################################
# geneset heatmap
#####################################################
# Inputs: pathid = "hsa04610" # KEGG PATHWAY ID
#         eset =  esetDN # EXPRESSION SET OBJECT
#         ColSideColors =  A Color vector based on the control and case groups

#-----------------------------
graphics.off()
#source("Rcode/funcDefs/drawGenesetHeatmap.R")
library("org.Hs.eg.db")
library("gplots")
library("RColorBrewer")

load("Metadata/gbyp.rda")

# Define the function based on Duke Transcriptome data
#------------------------------------------------------
#pathid = "hsa04610" # complement and coagulation cascade
#eset = esetDN # duke data
#ColSideColors =   c(rep("green2",10),rep("darkblue",9),rep("red",9),rep("brown",4)) # duke data



# Function Definition
#---------------------
drawGenesHeatmap = function(pathid = pathid,
                       eset = eset){
  egs = intersect(genes.by.pathway[[pathid]], featureNames(eset))
  gnames = AnnotationDbi::select(org.Hs.eg.db,
                  keys=egs,
                  column=c("SYMBOL","GENENAME"))
  labRowStr = with(gnames,
                   paste0(GENENAME, " | ", ENTREZID))
  # Define column side colors
  if(identical(eset$ID, esetDN$ID)){
    ColSideColors =   c(rep("green2",10),
                        rep("darkblue",9),
                        rep("red",9),rep("brown",4)) # duke data
  }else{
    ColSideColors = rep("green2", ncol(eset))
    ColSideColors[which(eset$Group!="Zcontrol")] =
      brewer.pal(n=length(which(eset$Group!="Zcontrol")), "Spectral")
  }
  if(identical(eset$ID, esetDN$ID)){
    colsep = c(10, 19, 28, 32)
  }else{colsep = 1:ncol(eset)}
  ttl.str = sapply(strsplit(pathways.list[[paste0("path:",
                                                  pathid)]],
                            " - Homo sapiens (human)", fixed = T),
                   "[[", 1)

  x11()
  par(cex.main=0.65)
  heatmap.2(exprs(eset[egs,]),
            Rowv=FALSE,
            Colv=FALSE,
            dendrogram="none",
            scale="row",
            trace="none",
            density.info="none",
            keysize = 0.91,
            key.title = " ",
            key.xlab = "Gene expression\n (log2)",
            key.ylab = NULL,
            col = bluered(256),
            margins = c(5,25),
            main = paste0(ttl.str, "\n(", pathid, ")"),
            labRow = labRowStr,
            ColSideColors = ColSideColors,
            sepwidth = c(0.1,0.1),
            colsep = colsep,
            sepcolor = "white",
            rowsep = 1:nrow(eset))
}

