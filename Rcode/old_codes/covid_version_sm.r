# This code compares genes common between two COVID datasets
# Then draws comparitive boxplots for theose common genes
# samanwoy "Wed Jun 10 14:43:03 2020"
# Modified by SKM on 03 Jul 2020
# 14 Jul 2020: (1) Group-level testing for up-regulation
#       establishes up-regulation of NETosis, not cytokines
#       (2) Validatation of up-regulation in lung
#       (3) Sample-level plots shows higher level in severe
#       (4) Trajectories: transient or sustained up-regulation
#############################################################

# Location: SKM_Lab/COVID-19
rm(list = ls())
graphics.off()
gc()
options(digits = 2)

# Load the libraries
#--------------------
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
snames = sapply(strsplit(snames, split = "_Whole"),
                function(x)
                  x[[1]][1])
sampleNames(eset) = snames
rm(snames)

#######################################
# Get the genes of cytokine and NETosis
#######################################
egs = read.csv(file = "metadata/input_file_for_netosis_c_storm_gns.csv")
egs.n = unique(as.character(egs[egs$class == "netosis", "egids"]))
egs.n = intersect(egs.n, featureNames(eset))
egs.c = unique(as.character(egs[egs$class == "cStorm", "egids"]))
egs.c = intersect(egs.c, featureNames(eset))
rm(egs)

######################################
# Result 01: NETosis is up-regulated in COVID-19
######################################
# Permutation histogram showing up-regulation
# of NETosis genes in blood of COVID-19 patients
# but not cytokine genes
p1 = drawPermutHistgg(
  geneids = egs.n,
  eset = ncov[[1]],
  titlestr = "NETosis Gene Set Expresson in Blood",
  direction = "Up"
)
graphics.off()

p2 = drawPermutHistgg(
  geneids = egs.c,
  eset = ncov[[1]],
  titlestr = "Cytokine Gene Set Expresson in Blood",
  direction = "Up"
)
p3 = ggpubr::ggarrange(p1, p2, nrow = 2)


# Heat map
rtt = rowttests(eset[egs.n,], "Group")
o = order(rtt[, "dm"])
plotdat = exprs(eset[egs.n[o],])
labRowStr = AnnotationDbi::select(org.Hs.eg.db,
                                  keys = egs.n, column = c("SYMBOL"))[, "SYMBOL"]
ColSideColors = rep("red", ncol(eset))
ColSideColors[eset$Group == "Zcontrol"] = "green"
heatmap.2(
  plotdat,
  Rowv = F,
  Colv = F,
  dendrogram = "none",
  scale = "row",
  trace = "none",
  density.info = "none",
  keysize = 01.9,
  key.title = " ",
  key.ylab = NULL,
  col = bluered(256),
  margins = c(5, 8),
  main = "Expression of NETosis genes",
  labRow = labRowStr,
  ColSideColors = ColSideColors,
  sepwidth = c(0.1, 0.1),
  sepcolor = "white",
  rowsep = 1:nrow(plotdat)
)
rm(o, rtt, plotdat, labRowStr, ColSideColors)
graphics.off()

######################################
# Result 02: Validation of NETosis in Lung
######################################
# Validation of up-regulation in lung
# library(RColorBrewer)
# sel.egs = c("1378", "5199", "1511", "7097",
#             "7099", "6352", "3689", "1088")
# esetLung = ncov[[5]][sel.egs,]
# gsyms = as.character(links(org.Hs.egSYMBOL[sel.egs])[, "symbol"])
# xcon = rowMeans(exprs(esetLung[, esetLung$Group == "Zcontrol"]))
# xcas = rowMeans(exprs(esetLung[, esetLung$Group == "nCOV"]))
# gxp = data.frame("Control" = xcon, "COVID-19" = xcas)
# o = order(xcas)
# par(mar=c(6,6,4,2))
# barplot(
#   t(gxp[o,]),
#   beside = TRUE,
#   col = c("#FEE6CE", "#E6550D"),
#   names = gsyms[o],
#   las=2,
#   main = "NETosis gene expression in Lung",
#   ylab = "Log (gene expression)",
#   ylim = c(0, max(gxp) + 1)
# )
# legend(
#   x = "topleft",
#   legend = c("Control (2 individuals)", "COVID-19 (2 patients)"),
#   fill = c("#FEE6CE", "#E6550D"),
#   bty = "n"
# )
# box()
# rm(esetLung, xcon, xcas, gxp, gsyms, sel.egs, o)
# graphics.off()
######################################
# Result 02: Validation of NETosis in Lung
######################################
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


#############################################
# Association with respiratory function
# % O2 sat for case 1
#############################################
perc_o2 = read.table("metadata/o2_sat_case_1.txt",
                     header = T,
                     sep = "\t")
rownames(perc_o2) = perc_o2$day

is.case1 = eset$ID == "Case1" & eset$Day %in% c(2:11)
eset.sub = eset[, is.case1]
day = eset.sub$Day
o2 = perc_o2[as.character(day), "percent_o2_saturation"]

# selected NETosis genes with negative correlation
# with oxygen saturation
# FAT1 = "2195"
# COX2 = "6775079"
egs.sel = c("3586",
            "6347",
            "1378",
            "3684",
            "629",
            "5199",
            "7097",
            "7099",
            "6352",
            "3689")
for (eg in egs.sel) {
  x = as.numeric(exprs(eset.sub[eg,]))
  rval = cor(o2, x)
  par(mar = c(5, 4, 4, 4) + 0.5)
  plot(
    x = day,
    y = o2,
    pch = 16,
    col = "red",
    xlab = "Day post-illness onset",
    ylab = "",
    cex.lab = 1.5,
    cex.main = 3,
    main = paste0(get(eg, org.Hs.egSYMBOL)),
    type = "o",
    lwd = 4,
    cex = 2,
    axes = F
  )
  legend("topright",
         bty = "n",
         legend = paste0("r = ", formatC(rval, digits = 2)))
  box(lwd = 2)
  axis(1, lwd = 2)
  axis(
    2,
    col = "red",
    lwd = 2,
    cex = 1.1,
    las = 2
  )
  mtext(
    "Oxygen saturation (%)",
    side = 2,
    line = 3,
    col = "red",
    cex = 1.1
  )
  par(new = TRUE)
  plot(
    x = day,
    y = x,
    pch = 15,
    col = "blue",
    cex = 2,
    lwd = 4,
    type = "o",
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  axis(
    side = 4,
    at = pretty(range(x)),
    las = 2,
    col = "blue",
    lwd = 2,
    cex = 1.1
  )
  mtext(
    "Log (gene expression)",
    side = 4,
    line = 3,
    col = "blue",
    cex = 1.1
  )
  scan()
  #graphics.off()
}
#rm(eg, egs.sel, rval, x, perc_o2, is.case1, day, eset.sub, o2)





# selected NETosis genes with negative correlation
# with oxygen saturation
egs.sel = c("1378", "3684", "629",
            "5199", "7097", "7099",
            "6352", "3689")
# egs.sel = egs.n
# get Sample Medians for
# all the NETosis
# genes for that time point
#----------------------------
library(plot3D)

eg = "3569" #IL6
eg = "3586" #IL10

x = as.numeric(apply(exprs(eset.sub[egs.sel,]), 2, median))
y = as.numeric(exprs(eset.sub[eg,]))
z = o2

graphics.off()
par(mar = c(2, 2, 2, 2))
scatter3D(
  x,
  y,
  z,
  colvar = eset.sub$Day,
  phi = 20,
  theta = 40,
  cex = 4.5,
  col = heat.colors(length(x)),
  add = FALSE,
  pch = 20,
  bty = "g",
  expand = 0.51,
  type = "s",
  #main = "IFIT2",
  xlab = "NETosis",
  ylab = links(org.Hs.egSYMBOL[eg])[, "symbol"],
  zlab = "O2",
  ticktype = "detailed",
  clab = c("Day", "Post-infection", "")
)

library("plot3Drgl")
plotrgl()
movie3d(spin3d(axis = c(0, 0, 1)), duration = 6,
        dir = "Results/Cytokines_leading_to_NETosis/")

############################################
#  Heatmap for O2 saturaion - IL6  - NETosis
############################################
graphics.off()
# Add color pallet
library(RColorBrewer)
#set.seed(123)
colVec = brewer.pal(3, "Oranges")

library(ComplexHeatmap)

eg = "3569" # IL6

eg = "3586" # IL10

eg = "3433" # IFIT2

# Create Matrix
mat = exprs(eset.sub[egs.sel,]) # NETosis Genes
#mat = exprs(eset.sub[c("6348", "6351", "2919", "3627", "2920", "3553", "3557"),]) # Other Cytokines

# Normalize the matrix with Day-4 gene expression
#------------------------------------------------
d4 = mat[, "Case1_d4"]
names(d4) = rownames(mat)
mat = (mat-d4)+1


# Add rownames as Gene symbol
rownames(mat) = as.character(unlist(mget(rownames(mat),
                                         org.Hs.egSYMBOL)))
# Add column names
colnames(mat) = paste0("Day-",
                       eset.sub$Day)

get(eg, org.Hs.egSYMBOL)

# Add column Barplot
column_ha = HeatmapAnnotation(
  #O2_Saturation = o2,
  gp = gpar(col = "black"),
  IL6 = anno_barplot(as.numeric(exprs(eset.sub[eg,]))),
  height = unit(3, "cm")
)




# draw Heatmap
ht_opt(
  heatmap_column_names_gp = gpar(fontface = "italic"),
  heatmap_column_title_gp = gpar(fontsize = 10),
  legend_border = "black",
  heatmap_border = TRUE,
  annotation_border = TRUE
)
Heatmap(
  mat,
  name = "Gexp",
  top_annotation = column_ha,
  cluster_columns = FALSE,
  col = colVec,
  rect_gp = gpar(col = "black"),
  row_title = "Genes",
  column_title = "Gene Expression Matrix",
  column_title_side = "bottom",
  row_title_side = "left",
)


######################################
# Result 03: Sample level testing
#     Is NETosis gene set up in the sample
######################################
# get mean control gene expression
xcon = rowMeans(exprs(eset[egs.n, eset$Group == "Zcontrol"]))
gexp = exprs(eset[egs.n, eset$Group == "nCOV"])
gexp = gexp - xcon
disease.severity = eset[, colnames(gexp)]$Severity
mycols = c("#FEE6CE", "#E6550D")[as.numeric(disease.severity == "Severe") + 1]
x11()
par(mar = c(6, 4, 3, 2), lend = 2, mfrow = c(2, 1))
boxplot(
  gexp,
  col = mycols,
  outline = F,
  las = 2,
  main = "NETosis gene experssion in COVID-19 patients",
  ylab = "\nLog (normalized gene expression)"
)
abline(h = 0, lwd = 1, lty = 3)
legend(
  x = "topright",
  legend = c("Moderate", "Severe"),
  fill = c("#FEE6CE", "#E6550D"),
  bty = "n"
)
rm(xcon, gexp, disease.severity, mycols)
#graphics.off()

library(openxlsx)
# load the data
dat = read.xlsx("metadata/Input_data_for_correlation_analysis_CIBERSORTx_Job9_Results_duke_data.xlsx",
                rowNames = T)

# Find the columns which are lymphocytes
#----------------------------------------
which.lymph = grep("_Lymphocyte", colnames(dat))
which.case = which(dat$Group!="Control")

# Create plot dat
#-----------------
nlr.vec = dat[which.case,]$Neutrophils/
  rowMeans(data.matrix(dat[which.case,which.lymph]))
#x11()
par(mar=c(8,4,2,2))
barplot(nlr.vec,
        col = c("#FEE6CE", "#E6550D")[factor(dat[which.case,]$Group)],
        las=2,
        main = "Neutrophil to Leukocyte Ratio in Duke cohort",
        ylab = "Relative fraction of \nNLR obtained by Deconvolution",
        names.arg = paste0(dat[which.case,]$PTID,
                           "_",
                           dat[which.case,]$Time))
box()
legend("topright",
       c("Moderate", "Severe"),
       fill = c("#FEE6CE", "#E6550D"),
       bty = "n"
)



######################################
# Result 04: Association with severity
######################################
# Create a small expression set with disease samples only
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
                 column = c("SYMBOL"))[, "SYMBOL"]
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
rm(
  i,
  pval.str,
  pvals,
  sd.mild,
  sd.severe,
  which.mild,
  which.severe,
  x.mild,
  x.severe,
  xpos.mild,
  xpos.severe,
  ymin.mild,
  ymin.severe
)
graphics.off()


######################################
# Result 05: Temporal change in gene expression
######################################
d6 = data.frame("Mild" = as.numeric(exprs(eset[egs.n, "Case2_d6"])),
                "Severe" = as.numeric(exprs(eset[egs.n, "Case1_d6"])))
d7 = data.frame("Mild" = as.numeric(exprs(eset[egs.n, "Case2_d7"])),
                "Severe" = as.numeric(exprs(eset[egs.n, "Case1_d7"])))
d8 = data.frame("Mild" = as.numeric(exprs(eset[egs.n, "Case2_d8"])),
                "Severe" = as.numeric(exprs(eset[egs.n, "Case1_d8"])))
d9 = data.frame("Mild" = as.numeric(exprs(eset[egs.n, "Case2_d9"])),
                "Severe" = as.numeric(exprs(eset[egs.n, "Case1_d9"])))
rownames(d6) = rownames(d7) = rownames(d8) = rownames(d9) = gsyms.n
o = order(rowMeans(data.frame(d6, d7, d8, d9)))
d6 = d6[o,]
d7 = d7[o,]
d8 = d8[o,]
d9 = d9[o,]

plot(
  x = 1:2,
  y = c(d6[1,]),
  xlim = c(1, 8),
  ylim = range(data.frame(d6, d7, d8, d9)),
  type = "o",
  pch = 16,
  lwd = 2,
  axes = F,
  xlab = "",
  ylab = "Log (Gene expression)",
  main = "Up-regulation of NETosis gene expression with severe illness"
)
axis(2)
lines(
  x = 3:4,
  y = c(d7[1,]),
  type = "o",
  pch = 16,
  lwd = 2
)
lines(
  x = 5:6,
  y = c(d8[1,]),
  type = "o",
  pch = 16,
  lwd = 2
)
lines(
  x = 7:8,
  y = c(d9[1,]),
  type = "o",
  pch = 16,
  lwd = 2
)
for (i in 2:length(egs.n)) {
  lines(
    x = 1:2,
    y = c(d7[i,]),
    type = "o",
    pch = 16,
    lwd = 2
  )
  lines(
    x = 3:4,
    y = c(d7[i,]),
    type = "o",
    pch = 16,
    lwd = 2
  )
  lines(
    x = 5:6,
    y = c(d8[i,]),
    type = "o",
    pch = 16,
    lwd = 2
  )
  lines(
    x = 7:8,
    y = c(d9[i,]),
    type = "o",
    pch = 16,
    lwd = 2
  )
}
box()
axis(1, tick = F, labels = F)
mtext(
  text = rep(c("Mild", "Severe"), 4),
  side = 1,
  at = 1:8,
  las = 2,
  line = 0.5
)
mtext(
  text = paste0("Day ", 6:9),
  side = 1,
  at = seq(from = 1.5, to = 7.5, by = 2),
  line = 3
)

######################################
# Result 06: Association of coagulopathy with severity
######################################
egs.coag = net.list$coagulation
egs.coag = intersect(egs.coag, featureNames(eset))
eset.sub = eset[, eset$Group != "Zcontrol"]
grp = eset.sub$Severity
grp[grp == "Moderate"] = "Zmild"
eset.sub$Group = as.factor(grp)
drawPermutHistgg(
  eset = eset.sub,
  direction = "Up",
  geneids = egs.coag,
  titlestr = "Up-regulation of coagulopathy module in severe COVID-19"
)
graphics.off()
