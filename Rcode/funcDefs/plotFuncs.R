# functions for drawing plots
# drawGeneLinePlot: input - egid, temporal plots for Duke-NUS data
# drawGeneBarplot:  inout - egid, barplot for five studies
# drawGenesetBoxplot: input - vector of egids, 
#    boxplot for five studies (control/case, LFC, 
#    or sample-level)
# drawGenesetHeatmap: input - vector of egids,
#    heat map with genes as rows and all samples

# Author: Saroj Kant Mohapatra (and others)
# Address: National Institute of Biomedical Genomics, India

##########################################################
# drawGenesetHeatmap
##########################################################
# For drawing heat map of the gene set
# Mean gene expression in control and nCOVID
# are plotted separately
drawGenesetHeatmap = function(egs, ...) {
  commong = which(table(unlist(sapply(ncov, function(eset) 
    intersect(egs, featureNames(eset)))))==5) %>% names
  gxp.lst = sapply(ncov, function(eset) {
    gset = eset[commong,]
    xcon = rowMeans(exprs(gset[, gset$Group=="Zcontrol"]))
    gxp = exprs(gset[,gset$Group!="Zcontrol"])-xcon
  })
  gxp.df = gxp.lst[[1]]
  for(i in 2:length(gxp.lst)) {
    gxp.df = cbind(gxp.df, gxp.lst[[i]])
  }
  gsyms = links(org.Hs.egSYMBOL[egs])[,"symbol"]
  gnames = links(org.Hs.egGENENAME[egs])[,"gene_name"]
  genestr = paste(gnames, egs, gsyms, sep=" | ")
  
  heatmap.2(gxp.df, trace="none", labRow=genestr,
            margins = c(10,35), key=0.2,
            col=colorRampPalette(c("darkgreen","white","red"))(256), ...)
} 

##########################################################
# drawGenesetBoxplot
##########################################################
# For drawing boxplot of the gene sets
# Mean gene expression in control and nCOVID 
# or sample-level box plots

drawGenesetBoxplot = function(egs, plotLFC=F, plotSample=F, ...) {
  if(plotSample==TRUE) {
    # for each study, data are shifted to median control value
    plotdat = list()
    counter = 0
    mycols = c()
    for(i in 1:5) {
      eset = ncov[[i]]
      sel = intersect(egs, featureNames(eset))
      eset = eset[sel,]
      mycols = c(mycols, c("red","green")[as.numeric(eset$Group=="Zcontrol")+1])
      mycols = c(mycols,"white")
      xcon = exprs(eset[, eset$Group=="Zcontrol"]) %>%
        rowMedians %>% median
      gxp = exprs(eset)
      gxp = gxp-xcon
      for(j in 1:as.numeric(ncol(eset))) {
        counter = counter+1
        plotdat[[counter]] = gxp[,j]
        names(plotdat)[counter] = colnames(eset)[j]
      }
      counter = counter+1
      plotdat[[counter]]=NA
    }
    mynames = names(plotdat)
    mynames = gsub("Whole.","",mynames)
    par(mar=c(13,4,4,2))
    boxplot(plotdat, names=mynames, 
            ylab="Log (Gene expression)",
            las=2, col=mycols, ...)
    
  } else {
    res = lapply(ncov, function(eset) {
      sel = intersect(egs, featureNames(eset))
      gset = eset[sel,]
      x = rowMeans(exprs(gset[, gset$Group=="Zcontrol"]))
      y = rowMeans(exprs(gset[, gset$Group=="nCOV"]))
      return(data.frame("Con"=x, "nCov"=y))
    })
    if(plotLFC==TRUE) {
      plotdat=sapply(res, function(x) {lfc=x$nCov-x$Con; return(lfc)})
      boxplot(plotdat, ylab="Log (Gene expression level)", ...)
      abline(h=0)
    } else {
      plotdat = vector("list",15)
      plotNames = rep(c("Control","nCovid",""), 5)
      for(i in 1:5) {
        plotdat[[3*i-2]]=res[[i]]$Con
        plotdat[[3*i-1]]=res[[i]]$nCov
      }
      boxplot(plotdat, names=plotNames, col=c("green","red","white"), 
              las=2, ylab="Log (Gene expression level)", ...)
      mtext(c("DN.Blood","WU.PBMC","WU.BALF","NY.LungEpith","NY.LungBiops"), at=seq(1.5, 15, by=3), side=1, line=4)
    }    
  }
}

##########################################################
# drawGeneBarplot
##########################################################
# bar plot for mean gene expression for a single gene
drawGeneBarplot = function(eg) {
  if(eg%in%featureNames(ncov[[1]]) | eg%in%featureNames(ncov[[2]])) {
    lfc=sapply(ncov, function(eset){
      res = c("Con"=NA,"nCOV"=NA)
      if(eg%in%featureNames(eset)) {
        is.control = eset$Group=="Zcontrol"
        is.ncov = eset$Group=="nCOV"
        xcon=mean(exprs(eset[eg,is.control]), na.rm=T)
        xcase=mean(exprs(eset[eg,is.ncov]), na.rm=T)
        res = c("Con"=xcon, "nCOV"=xcase)
      }
      return(res)
    })
    gsym = get(eg, org.Hs.egSYMBOL)
    gname = get(eg, org.Hs.egGENENAME)
    titlestr = paste(eg, gsym, gname, sep=" | ")
    barplot(lfc, main=titlestr, beside=TRUE, 
            ylab="Relative gene expression (log scale)", 
            col=c("green","red"))
  } else {
    cat("Gene expression not measured in studies on blood or PBMC.\n")
  }
}

##########################################################
# drawGeneLinePlot
##########################################################
# Draw line plot for DUME-NUS COVID-19 data
# It is assumed that an expression set called esetDN
# is already available in the workspace
# Author: Samanwoy and SKM
# 02 Jul 2020

drawGeneLinePlot = function( egid = "3458"){
  gsym = as.character(unlist(mget(egid, org.Hs.egSYMBOL)))
  
  # Get control mean data
  ctrl.mdn = median(exprs(esetDN[egid,
                                 esetDN$Group=="Zcontrol"]))
  
  case.vec.1 = exprs(esetDN[egid,
                            esetDN$ID=="Case1"])
  day.vec.1 = sapply(strsplit(colnames(case.vec.1), "_"), "[[", 2)
  day.vec.1 = as.numeric(substr(day.vec.1, 2, nchar(day.vec.1)))
  
  case.vec.2 = exprs(esetDN[egid,
                            esetDN$ID=="Case2"])
  day.vec.2 = sapply(strsplit(colnames(case.vec.2), "_"), "[[", 2)
  day.vec.2 = as.numeric(substr(day.vec.2, 2, nchar(day.vec.2)))
  
  case.vec.3 = exprs(esetDN[egid,
                            esetDN$ID=="Case3"])
  day.vec.3 = sapply(strsplit(colnames(case.vec.3), "_"), "[[", 2)
  day.vec.3 = as.numeric(substr(day.vec.3, 2, nchar(day.vec.3)))
  
  #----------------
  # Case 1
  #----------------
  x <- day.vec.1
  y <- case.vec.1
  
  #------------------
  # Plot
  #------------------
  titlestr = paste0(get(egid, org.Hs.egGENENAME), " | ", 
                    get(egid, org.Hs.egSYMBOL), " | ", egid)
  par(mar=c(6.5,5,4,2))
  plot(x, y,
       type="n",
       pch=19,
       col = "red1",
       xaxt = "n",
       main= titlestr,
       xlab="Day (post-onset)",
       ylab="Log (gene expression)",
       ylim= range(c(case.vec.1, case.vec.2, case.vec.3)) +c(-0.2, 0.3),
       xlim=c(min(x), 20))
  title(sub="Data from Duke-NUS center (Ong 2020)", line=5)
  rect(xleft=0,
       ybottom=ctrl.mdn-0.1,
       xright= 28+max(x),
       ytop=ctrl.mdn+0.1,
       col="gray70",
       border=NA)
  abline(h= ctrl.mdn)
  lines(x=x,
        y=y,
        type="o",
        col="red1",
        pch=19,
        cex=1.2,
        lwd=3)
  xtick<-seq(0,
             28,
             by=1)
  axis(side=1,
       at=xtick,
       labels = 0:28)
  box()
  #------------------
  # Case 2
  #------------------
  x = day.vec.2
  y= case.vec.2
  lines(x=x,
        y=y,
        type="o",
        col="blue",
        pch=19,
        cex=1.2,
        lwd=3)
  #------------------
  # Case 3
  #------------------
  x = day.vec.3
  y= case.vec.3
  lines(x=x,
        y=y,
        type="o",
        col="darkgreen",
        pch=19,
        cex=1.2,
        lwd=3)
  
  #------------------
  # Add legend
  #------------------
  legend("topright",
         legend = c("1 (Severe)", "2 (Moderate)", "3 (Moderate)"),
         col = c("red1", "blue", "darkgreen"),
         bty = "n",
         lty=1,
         lwd = 4,
         cex = 1,
         horiz = T)
}

