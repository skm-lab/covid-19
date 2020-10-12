# Code for re-analysis of published host
# transcriptome response to SARS-CoV-2
#############################################################

# Instruction for running the code
#---------------------------------
# After you clone git repository
# Location: /covid-19

# Load the libraries
#---------------------
library(Biobase)
library(gplots)
library(genefilter)
library(dplyr)
library(org.Hs.eg.db)

# function definitions
source("Rcode/funcDefs/plotFuncs.R")
source("Rcode/funcDefs/drawPermutHistgg.R")

# Load data and metadata
source("Rcode/getData.R")
eset = ncov[[1]]
snames = sampleNames(eset)
snames = sapply(strsplit(snames, split="_Whole"),
                function(x) x[[1]][1])
sampleNames(eset) = snames
rm(snames)

# Get the genes of cytokine and NETosis
#-------------------------------------------
egs = read.csv(file="metadata/input_file_for_netosis_c_storm_gns.csv")
egs.n = unique(as.character(egs[egs$class=="netosis","egids"]))
egs.n = intersect(egs.n, featureNames(eset))
egs.c = unique(as.character(egs[egs$class=="cStorm","egids"]))
egs.c = intersect(egs.c, featureNames(eset))
rm(egs)


#-----------------------------------------------------------
# Result 01: Figure 1A, 1B: NETosis is up-regulated in COVID-19
#-----------------------------------------------------------
# Permutation histogram showing up-regulation
# of NETosis genes in blood of COVID-19 patients
# but not cytokine genes
p1 = drawPermutHistgg(
  geneids = egs.n,
  eset = ncov[[1]],
  titlestr = "NETosis Gene Set Expresson in Blood",
  direction = "Up"
)

p2 = drawPermutHistgg(
  geneids = egs.c,
  eset = ncov[[1]],
  titlestr = "Cytokine Gene Set Expresson in Blood",
  direction = "Up"
)
p3 = ggpubr::ggarrange(p1, p2, nrow = 2)

#-----------------------------------------------------------
# Result 02: Figure 2A, NETosis Heatmap
#-----------------------------------------------------------
rtt = rowttests(eset[egs.n,], "Group")

which.control = which(!eset$ID%in%c("Case1","Case2","Case3"))
which.case1 = which(eset$ID=="Case1")
which.case2 = which(eset$ID=="Case2")
which.case3 = which(eset$ID=="Case3")
sel = c(which.control, which.case2, which.case3, which.case1)
ColSideColors=c(rep("green",length(which.control)),
                rep("cyan", length(c(which.case2, which.case3))),
                rep("red", length(which.case1)))
o = order(rtt[,"dm"])
plotdat = exprs(eset[egs.n[o],sel])
labRowStr = AnnotationDbi::select(org.Hs.eg.db,
      keys=rownames(plotdat), columns=c("SYMBOL"))[,"SYMBOL"]
heatmap.2(plotdat, Rowv=F, Colv = F,
          dendrogram = "none", scale="row",
          trace="none",
          density.info="none",
          keysize = 1.5,
          key.title = " ",
          key.ylab = NULL,
          col = bluered(256),
          margins = c(5,8),
          main = "Expression of NETosis genes",
          labRow = labRowStr,
          ColSideColors = ColSideColors,
          sepwidth = c(0.1,0.1),
          sepcolor = "white",
          rowsep = 1:nrow(plotdat))
rm(o,rtt,plotdat, labRowStr, ColSideColors, which.control, which.case1,
   which.case2, which.case3, sel)
graphics.off()


#-----------------------------------------------------------
# Result 03: Figure 3; NETosis in Lung
#-----------------------------------------------------------
# Validation of up-regulation in lung
library(RColorBrewer)
sel.egs = c("1378", "5199", "1511", "7097",
            "7099", "6352", "3689", "1088")
esetLung = ncov[[5]][sel.egs,]
gsyms = as.character(links(org.Hs.egSYMBOL[sel.egs])[, "symbol"])

featureNames(esetLung) = gsyms
xcon = rowMeans(exprs(esetLung[, esetLung$Group == "Zcontrol"]))
xcas = rowMeans(exprs(esetLung[, esetLung$Group == "nCOV"]))

# Calculate p-value of t-test between mild and severe illness
# for every NETosis gene
pvals = apply(exprs(esetLung), 1, function(x) {
  xcon = x[esetLung$Group != "nCOV"]
  xcas = x[esetLung$Group == "nCOV"]
  t.test(xcas, xcon,
         alternative = "greater")$p.value
})

# bar graph of gene expression
which.con = which(esetLung$Group == "Zcontrol")
which.cas = which(esetLung$Group == "nCOV")
x.con = rowMeans(exprs(esetLung[, which.con]))
x.cas = rowMeans(exprs(esetLung[, which.cas]))
sd.con = rowSds(exprs(esetLung[, which.con]))
sd.cas = rowSds(exprs(esetLung[, which.cas]))
plotdat = data.frame(x.con, x.cas)
errordat = data.frame(sd.con, sd.cas)


# Calculate p-value of t-test between mild and severe illness
# for every NETosis gene
pvals = apply(exprs(esetLung), 1, function(x) {
  xcon = x[esetLung$Group != "nCOV"]
  xcas = x[esetLung$Group == "nCOV"]
  t.test(xcas, xcon,
         alternative = "greater")$p.value
})

#re order the variables and plotdat
#------------------------------------
o = order(xcon)
pvals = pvals[o]
errordat = errordat[o,]
plotdat = plotdat[o,]

# Draw the plot
#--------------
b = barplot(
  t(plotdat),
  beside = T,
  col = c("#FEE6CE", "#E6550D"),
  las = 2,
  ylim = c(0, max(plotdat) + 5.5),
  main = "NETosis gene expression in Lung",
  ylab = " Log (gene expression)"
)
legend(
  x = "topleft",
  legend = c("Healthy Lung", "COVID-19 Lung"),
  fill = c("#FEE6CE", "#E6550D"),
  bty = "n"
)
box()
pval.str = ifelse(pvals < 0.001, "***",
                  ifelse(pvals < 0.01, "**",
                         ifelse(pvals < 0.1, "*", "")))

for (i in 1:nrow(plotdat)) {
  text(
    x = mean(b[, i]),
    y = plotdat[i, "x.cas"] + 1,
    labels = pval.str[i],
    adj = 0.5
  )
  xpos.con = b[1, i]
  ymin.con = plotdat[i, "x.con"]
  ymax.con = plotdat[i, "x.con"] + errordat[i, "sd.con"]
  lines(x = c(xpos.con, xpos.con),
        y = c(ymin.con, ymax.con))
  lines(
    x = c(xpos.con - 0.25, xpos.con + 0.25),
    y = c(ymax.con, ymax.con)
  )

  xpos.cas = b[2, i]
  ymin.cas = plotdat[i, "x.cas"]
  ymax.cas = plotdat[i, "x.cas"] + errordat[i, "sd.cas"]
  lines(x = c(xpos.cas, xpos.cas),
        y = c(ymin.cas, ymax.cas))
  lines(
    x = c(xpos.cas - 0.25, xpos.cas + 0.25),
    y = c(ymax.cas, ymax.cas)
  )
}
rm(
  i,
  pval.str,
  pvals,
  sd.cas,
  which.cas,
  x.cas,
  xpos.con,
  xpos.cas,
  ymin.cas
)
graphics.off()

#-----------------------------------------------------------
# Result 04: Figure 5: NETosis Temporal line plot with O2 saturation
#-----------------------------------------------------------
# Association with respiratory function
# % O2 sat for case 1
perc_o2 = read.table("metadata/o2_sat_case_1.txt",
              header=T, sep="\t")
rownames(perc_o2) = perc_o2$day

is.case1 = eset$ID=="Case1" & eset$Day%in%c(2:11)
eset.sub = eset[,is.case1]
day = eset.sub$Day
o2 = perc_o2[as.character(day), "percent_o2_saturation"]

# selected NETosis genes with negative correlation
# with oxygen saturation
egs.sel = c("1378","3684","629",
            "5199","7097","7099",
            "6352","3689","1088")
x11()
par(mfrow=c(3,3))
for(eg in egs.sel) {
  x=as.numeric(exprs(eset.sub[eg,]))
  rval = cor(o2,x)
  par(mar = c(5, 6, 4, 4) + 0.5)
  plot(x=day,
       y=o2,
       pch = 16,
       col = "red",
       xlab="Day post-onset of illness",
       ylab="",
       cex.lab=1,
       cex.main=2,
       main=paste0(get(eg,org.Hs.egSYMBOL), ": ",
                   paste0(" (r = ", formatC(rval, digits=2)), ")"),
       type="o",
       lwd=2,
       cex=1.5,
       axes=F)
  #legend("topright", bty="n",
  #  legend=paste0("r = ", formatC(rval, digits=2)))
  box(lwd=2)
  axis(1, lwd=2)
  axis(2,
       col="red",
       lwd=2,
       cex=1.1,
       las=2)
  mtext("Oxygen saturation (%)",
        side = 2,
        line = 3,
        col="red",
        cex=0.8)
  par(new = TRUE)
  plot(x=day,
       y=x,
       pch = 15,
       col = "blue",
       cex=1.5,
       lwd=2,
       type="o",
       axes = FALSE,
       xlab = "",
       ylab = "")
  axis(side = 4,
       at = pretty(range(x)),
       las=2,
       col="blue",
       lwd=2,
       cex=1.1)
  mtext("Log (gene expression)",
        side = 4,
        line = 3,
        col="blue",
        cex=0.8)
  #scan()
  #graphics.off()
}
rm(eg,
   egs.sel,
   rval,
   x,
   perc_o2,
   is.case1,
   day,
   eset.sub,
   o2)


#-----------------------------------------------------------
# Result 05: Figure 4A, 4B: Sample level
# NETosis expression with NLR
#-----------------------------------------------------------
x11();par(mfrow=c(2,1), mar=c(7,5,4,2), lend=2)

# bar graph of NLR; do not draw
# only to get the position of the barplot
plotLogical = FALSE
source("Rcode/analyse_NLR_in_Duke.R")
#mtext(side=3, text="A", cex=3, at=-3.5)

# boxplot of get mean control gene expression
xcon = rowMeans(exprs(eset[egs.n, eset$Group=="Zcontrol"]))
gexp = exprs(eset[egs.n, eset$Group=="nCOV"])
gexp = gexp-xcon
disease.severity = eset[,colnames(gexp)]$Severity
mycols = c("blue","red")[as.numeric(disease.severity=="Severe")+1]
boxplot(gexp,
        col=mycols,
        las=2,
        outline = F,
        #main="NETosis gene experssion in COVID-19 patients",
        ylab="Average  Log2 gene \nexpression of NETosis",
        at = as.vector(b))
abline(h=0, lwd=2, lty=3)
legend(x="topright", legend=c("Moderate Illness","Severe Illness"),
       fill=c("blue","red"), bty="n")
mtext(side=3, text="A", cex=3, at=-2, line=1.5)
rm(xcon, gexp, disease.severity, mycols)

plotLogical = TRUE
source("Rcode/analyse_NLR_in_Duke.R")
box()
mtext(side=3, text="B", cex=3, at=-2, line=1.5)
#graphics.off()


#-----------------------------------------------------------
# Result 06: Figure 2B: Sample level
# NETosis expression Association with severity
#-----------------------------------------------------------
which.mild = which(eset$Severity == "Moderate")
which.severe = which(eset$Severity == "Severe")
eset.sub = eset[egs.n, c(which.mild, which.severe)]
rm(which.mild, which.severe)

# ORder the genes by increasing gene expression
o = order(rowMeans(exprs(eset.sub)))
eset.sub = eset.sub[o,]
rm(o)

# rename the mild samples as Zmild (instead of simply mild)
# (helps with rowttests)
eset.sub$Severity[eset.sub$Severity == "Moderate"] = "Zmild"

# assign gene symbols as feature names
gsyms.n = select(org.Hs.eg.db,
                 keys = featureNames(eset.sub),
                 columns = c("SYMBOL"))[, "SYMBOL"]
featureNames(eset.sub) = gsyms.n
rm(gsyms.n)

# Calculate p-value of t-test between mild and severe illness
# for every NETosis gene
pvals = apply(exprs(eset.sub), 1, function(x) {
  xmild = x[eset.sub$Severity != "Severe"]
  xsevere = x[eset.sub$Severity == "Severe"]
  t.test(xsevere, xmild,
         alternative = "greater")$p.value
})

# bar graph of gene expression
which.mild = which(eset.sub$Severity == "Zmild")
which.severe = which(eset.sub$Severity == "Severe")
x.mild = rowMeans(exprs(eset.sub[, which.mild]))
x.severe = rowMeans(exprs(eset.sub[, which.severe]))
sd.mild = rowSds(exprs(eset.sub[, which.mild]))
sd.severe = rowSds(exprs(eset.sub[, which.severe]))
plotdat = data.frame(x.mild, x.severe)
errordat = data.frame(sd.mild, sd.severe)
b = barplot(
  t(plotdat),
  beside = T,
  col = c("#FEE6CE", "#E6550D"),
  las = 2,
  ylim = c(0, max(plotdat) + 1.5),
  main = "Association of NETosis gene expression with disease severity",
  ylab = " Log (gene expression)"
)
legend(
  x = "topleft",
  legend = c("Moderate Illness", "Severe Illness"),
  fill = c("#FEE6CE", "#E6550D"),
  bty = "n"
)
box()
pval.str = ifelse(pvals < 0.001, "***", ifelse(pvals < 0.01, "**", ifelse(pvals <
                                                                            0.1, "*", "")))

for (i in 1:length(egs.n)) {
  text(
    x = mean(b[, i]),
    y = plotdat[i, "x.severe"] + 1,
    labels = pval.str[i],
    adj = 0.5
  )
  xpos.mild = b[1, i]
  ymin.mild = plotdat[i, "x.mild"]
  ymax.mild = plotdat[i, "x.mild"] + errordat[i, "sd.mild"]
  lines(x = c(xpos.mild, xpos.mild),
        y = c(ymin.mild, ymax.mild))
  lines(
    x = c(xpos.mild - 0.25, xpos.mild + 0.25),
    y = c(ymax.mild, ymax.mild)
  )

  xpos.severe = b[2, i]
  ymin.severe = plotdat[i, "x.severe"]
  ymax.severe = plotdat[i, "x.severe"] + errordat[i, "sd.severe"]
  lines(x = c(xpos.severe, xpos.severe),
        y = c(ymin.severe, ymax.severe))
  lines(
    x = c(xpos.severe - 0.25, xpos.severe + 0.25),
    y = c(ymax.severe, ymax.severe)
  )
}

rm(i, pval.str, pvals, sd.mild, sd.severe,
   which.mild, which.severe, x.mild, x.severe,
   xpos.mild, xpos.severe, ymin.mild, ymin.severe,
   errordat, eset.sub, b)
graphics.off()

#-----------------------------------------------------------
# Result 07: Figure 2C: Temporal change in gene expression
#--------------------------------------------------------
plotdat = data.frame("D6.Moderate"=as.numeric(exprs(eset[egs.n,"Case2_d6"])),
                     "D6.Severe"=as.numeric(exprs(eset[egs.n,"Case1_d6"])),
                     NA,
                     "D7.Moderate"=as.numeric(exprs(eset[egs.n,"Case2_d7"])),
                     "D7.Severe"=as.numeric(exprs(eset[egs.n,"Case1_d7"])),
                     NA,
                     "D8.Moderate"=as.numeric(exprs(eset[egs.n,"Case2_d8"])),
                     "D8.Severe"=as.numeric(exprs(eset[egs.n,"Case1_d8"])),
                     NA,
                     "D9.Moderate"=as.numeric(exprs(eset[egs.n,"Case2_d9"])),
                     "D9.Severe"=as.numeric(exprs(eset[egs.n,"Case1_d9"])),
                     NA)
boxplot(plotdat, names=rep("",length(plotdat)),
          ylim=c(2,15), ylab="Log (gene expression level)",
          col=rep(c("blue","red","white"), 4), axes=F)
axis(2); box()
legend(x="topleft", legend=c("Moderate Illness","Severe Illness"),
       fill=c("blue","red"), bty="n")
for(i in 1:4) {
  mtext(text=paste0("Day_",i+5), side=1, line=1,
        at=c(2,5,8,11)[i])
}
title(sub="Days post-onset of COVID-19 illness", line=2.5)


d6 = data.frame("Moderate"=as.numeric(exprs(eset[egs.n,"Case2_d6"])),
                "Severe"=as.numeric(exprs(eset[egs.n,"Case1_d6"])))
d7 = data.frame("Moderate"=as.numeric(exprs(eset[egs.n,"Case2_d7"])),
                "Severe"=as.numeric(exprs(eset[egs.n,"Case1_d7"])))
d8 = data.frame("Moderate"=as.numeric(exprs(eset[egs.n,"Case2_d8"])),
                "Severe"=as.numeric(exprs(eset[egs.n,"Case1_d8"])))
d9 = data.frame("Moderate"=as.numeric(exprs(eset[egs.n,"Case2_d9"])),
                "Severe"=as.numeric(exprs(eset[egs.n,"Case1_d9"])))
gsyms.n = select(org.Hs.eg.db,
    keys=egs.n, columns=c("SYMBOL"))[,"SYMBOL"]
rownames(d6) = rownames(d7) = rownames(d8) = rownames(d9) = gsyms.n
res = sapply(list(d6,d7,d8,d9), function(x)
  apply(x, 1, function(y) y[2]-y[1]))
head(res[order(rowMedians(res), decreasing=TRUE),])


rm(res, d6,d7,d8,d9, plotdat, i, gsyms.n)

# Earlier code for line plots. Once the figures for the manuscript are fixed
# this chunk may be removed
# plotdat = list(d6,d7,d8,d9)
#
# o = order(rowMeans(data.frame(d6,d7,d8,d9)))
# d6=d6[o,]
# d7=d7[o,]
# d8=d8[o,]
# d9=d9[o,]
#
# par(mar=c(6,5,4,2))
# plot(x=1:2, y=c(d6[1,]), xlim=c(1,8), ylim=range(data.frame(d6,d7,d8,d9)),
#      type="o", pch=16, lwd=2, axes=F, xlab="", ylab="Log (Gene expression)",
#      main="Up-regulation of NETosis gene expression with severe illness")
# axis(2)
# lines(x=3:4, y=c(d7[1,]), type="o", pch=16, lwd=2)
# lines(x=5:6, y=c(d8[1,]), type="o", pch=16, lwd=2)
# lines(x=7:8, y=c(d9[1,]), type="o", pch=16, lwd=2)
# for(i in 2:length(egs.n)) {
#   lines(x=1:2, y=c(d7[i,]), type="o", pch=16, lwd=2)
#   lines(x=3:4, y=c(d7[i,]), type="o", pch=16, lwd=2)
#   lines(x=5:6, y=c(d8[i,]), type="o", pch=16, lwd=2)
#   lines(x=7:8, y=c(d9[i,]), type="o", pch=16, lwd=2)
# }
# box()
# axis(1, tick=F, labels = F)
# mtext(text=rep(c("Moderate","Severe"), 4), side=1, at=1:8, las=2, line=0.5)
# mtext(text=paste0("Day ", 6:9), side=1, at=seq(from=1.5, to=7.5, by=2),
#       line=4)
#
# pvals=rep(NA,4)
# names(pvals) = c("D6","D7","D8","D9")
# pvals["D6"] = t.test(d6[,1],d6[,2], paired=TRUE)$p.value
# pvals["D7"] = t.test(d7[,1],d7[,2], paired=TRUE)$p.value
# pvals["D8"] = t.test(d8[,1],d8[,2], paired=TRUE)$p.value
# pvals["D9"] = t.test(d9[,1],d9[,2], paired=TRUE)$p.value
# print(pvals)

#-----------------------------------------------------------
# Result 08: Figure 6: Association of complement pathway with severity
#-------------------------------------------------------------
egs.hsa04610 = genes.by.pathway[["hsa04610"]]

gsym.hsa04610 = select(org.Hs.eg.db, keys=egs.hsa04610,
      columns=c("SYMBOL", "GENENAME"))[,"GENENAME"]
which.complement = grep("complement", gsym.hsa04610)
egs.complement = egs.hsa04610[which.complement]
egs.complement = intersect(egs.complement, featureNames(eset))
eset.sub = eset[,eset$Group!="Zcontrol"]
grp = eset.sub$Severity
grp[grp=="Moderate"] = "Zmild"
eset.sub$Group = as.factor(grp)
drawPermutHistgg(eset=eset.sub, direction="Up", geneids=egs.complement,
                 titlestr="Up-regulation of complement gene set in severe COVID-19")

# Earlier code for complement and coagulation pathway together
# Once the figures for the manuscript are fixed
# egs.coag = net.list$coagulation
# egs.coag = genes.by.pathway[["hsa04610"]]
# egs.coag = intersect(egs.coag, featureNames(eset))
# eset.sub = eset[,eset$Group!="Zcontrol"]
# grp = eset.sub$Severity
# grp[grp=="Moderate"] = "Zmild"
# eset.sub$Group = as.factor(grp)
# drawPermutHistgg(eset=eset.sub, direction="Up", geneids=egs.coag,
#       titlestr="Up-regulation of complement and coagulation pathway in severe COVID-19")
# graphics.off()

rm(egs.hsa04610, grp, gsym.hsa04610, which.complement,
   egs.complement, eset.sub)

#---------------------------------------------------------
# Result 09: Figure 7: Which genes are most correlated with cytokine
# IL6, IL8, IL1b, TNFa
#---------------------------------------------------------
select(org.Hs.eg.db,
       keys=egs.c, column=c("SYMBOL"))[,"SYMBOL"]
egs.cytok = c("3576","3553","3569","7124")
gsym.cytok = select(org.Hs.eg.db,
   keys=egs.cytok, column=c("SYMBOL"))[,"SYMBOL"]


#> mget(egs, (org.Hs.egSYMBOL))
# 3576 "CXCL8"
# 3553 "IL1B"
# 3569 "IL6"
# 7124 "TNF"

# For severe illness: Case 1
early.severe = c("Case1_d4","Case1_d5","Case1_d6","Case1_d7")
late.severe = c("Case1_d11","Case1_d12","Case1_d18")

x.early.severe = exprs(eset[,early.severe])
cytok.early.severe = x.early.severe[egs.cytok,]
netos.early.severe = colMeans(x.early.severe[egs.n,])

x.late.severe = exprs(eset[,late.severe])
cytok.late.severe = x.late.severe[egs.cytok,]
netos.late.severe = colMeans(x.late.severe[egs.n,])

cor.early.severe = apply(cytok.early.severe, 1, function(x)
  cor(x, netos.early.severe))
cor.late.severe = apply(cytok.late.severe, 1, function(x)
  cor(x, netos.late.severe))

# For moderate illness: Case 2
early.Case2 = c("Case2_d6","Case2_d7")
late.Case2 = c("Case2_d11","Case2_d12","Case2_d13","Case2_d19")

x.early.Case2 = exprs(eset[,early.Case2])
cytok.early.Case2 = x.early.Case2[egs.cytok,]
netos.early.Case2 = colMeans(x.early.Case2[egs.n,])

x.late.Case2 = exprs(eset[,late.Case2])
cytok.late.Case2 = x.late.Case2[egs.cytok,]
netos.late.Case2 = colMeans(x.late.Case2[egs.n,])

cor.early.Case2 = apply(cytok.early.Case2, 1, function(x)
  cor(x, netos.early.Case2))
cor.late.Case2 = apply(cytok.late.Case2, 1, function(x)
  cor(x, netos.late.Case2))

# For moderate illness: Case 3
late.Case3 = c("Case3_d11","Case1_d12")

x.late.Case3 = exprs(eset[,late.Case3])
cytok.late.Case3 = x.late.Case3[egs.cytok,]
netos.late.Case3 = colMeans(x.late.Case3[egs.n,])

cor.late.Case3 = apply(cytok.late.Case3, 1, function(x)
  cor(x, netos.late.Case3))

# cytokine correlations
cor.cytok = data.frame(cor.early.Case2, cor.late.Case2, cor.late.Case3,
                       cor.early.severe, cor.late.severe)
rownames(cor.cytok) =  select(org.Hs.eg.db,
    keys=rownames(cor.cytok), columns=c("SYMBOL"))[,"SYMBOL"]


# hc with netosis and cytokine genes together
egs.nc = c(egs.n, egs.cytok)
gsyms.nc = select(org.Hs.eg.db,
                  keys=egs.nc, columns=c("SYMBOL"))[,"SYMBOL"]
labRow = paste0(gsyms.nc,
      c(rep("_NETO", length(egs.n)), rep("_CYTO", length(egs.cytok))))

which.case = which(eset$Group=="nCOV")
x11()
x = exprs(eset[egs.nc,which.case])
#rownames(x) = gsyms.nc
rowsep = 1:nrow(x)
heatmap.2(x, Rowv=T, Colv = F,
          dendrogram = "row", scale="row",
          trace="none",
          density.info="none",
          keysize = 1.5,
          key.title = " ",
          key.ylab = NULL,
          col = bluered(256),
          margins = c(5,7),
          main = "",
          labRow = gsyms.nc,
#          ColSideColors = ColSideColors,
          sepwidth = c(0.1,0.1),
          sepcolor = "white",
          rowsep = rowsep)

