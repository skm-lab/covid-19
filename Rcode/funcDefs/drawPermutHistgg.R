# This is a function that drwas the histogram
# of simulated enrichment scores by permutaion GSEA
# samanwoy, "Thu Jun 11 19:35:02 2020"
#-------------------------------------
library(genefilter)
library(ggplot2)
drawPermutHistgg <- function(keggid = "hsa04064", eset = eset, direction=NULL,
                                  fac = factor(as.character(eset$Group)),
                                  geneids=NULL,
                                  titlestr = NULL,
                                  nperm=10000,
                                  drawHisto=TRUE,
                                  lwd=3,
                                  cex=1.5) {
  if(is.null(geneids)) {
    pathegs = genes.by.pathway[[keggid]]
    #titlestr = KEGGPATHID2NAME[[strsplit(keggid, "hsa")[[1]][2]]]
    titlestr = titlestr
    } else {
      pathegs = geneids
      titlestr = titlestr
      }
  pathegs = intersect(pathegs, featureNames(eset))
  eset.sub = eset[pathegs,]
  rttStat = rowttests(eset.sub, fac, tstatOnly = TRUE)[["statistic"]]
  zobs = sum(rttStat)/sqrt(length(pathegs))

  permVec <- rep(0, nperm)
  set.seed(123)
  j <- 1L
  while (j < (nperm + 1)) {
    p1 <- sample(fac)
    rstatSim <- rowttests(eset.sub, p1, tstatOnly = TRUE)[["statistic"]]
    permVec[j] = sum(rstatSim)/sqrt(length(pathegs))
    j <- j + 1L
  }
  #permut_pval = ifelse(zobs<0, sum(permVec<zobs), sum(permVec>zobs))
  if(is.null(direction)) {
    permut_pval = ifelse(zobs<0, sum(permVec<zobs), sum(permVec>zobs))
  } else if(direction=="Up") {
    permut_pval = sum(permVec>zobs)
  } else if(direction=="Down") {
    permut_pval = sum(permVec<zobs)
  } else {
    cat("Direction must be Up, Down or NULL.\n")
  }
  permut_pval = permut_pval/nperm
  res <- c("p.val"=permut_pval,"zobs"=zobs)
  df <- data.frame("sim_score"= rep("simul_score", length(permVec)), "perm" = permVec)
  theme_update(plot.title = element_text(hjust = 0.5))
  p = ggplot(df, aes(x=perm)) +
    geom_histogram(color="black", fill="gray", binwidth = 0.75) +
    ylab("Frequency") +
    xlab("Simulated Pathway Score") +
    ggtitle(titlestr) +
    theme(legend.position="right",
          plot.margin = margin(2, 2, 2, 2, "cm")) +
    theme(plot.title = element_text(color="black",
                                    size=18, face="bold"),
          axis.title.x = element_text(color="steelblue",
                                      size=12, face="bold"),
          axis.title.y = element_text(color="steelblue",
                                      size=12, face="bold")) +
    theme(axis.text.y = element_text(size=10)) +
    geom_vline(xintercept=zobs, linetype="solid",
               color = "red",
               size=1.2) +
    geom_text(x=zobs+sign(zobs/2),
              angle= 90,
              y=0.2*max(hist(permVec,
                             plot=FALSE)$counts),
              label=c(paste0("p = ",
                             formatC(permut_pval, digits=1))),
              #labels=c("Observed Score", "Null distribution"),
              col=c("red"), size= 5,
              fontface="italic")
  #if(drawHisto==TRUE){
  #    print(permut_pval)
  #  } else {
  #    return(p)
  #  }
  return(p)
  }
