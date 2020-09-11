# This code performs survivor and non survivor
# barplots within each sepsis cohort
# for netosis and cytokine response genes
# Samanwoy "Fri Jul 17 20:46:08 2020"
####################################
rm(list=ls())
graphics.off()
#################################################
# load the first set of data for global COHORT-1
#################################################
library(ssgeosurv)
library(Biobase)
library(org.Hs.eg.db)
library(plotrix)

data(ss.list)
sel.studies = c(1:4,6,8)
ss.list = ss.list[sel.studies]

eset = ss.list[[1]]
study.id = names(ss.list)[1]

#### Load the Cytokine and NETosis genes ####
egs = read.csv(file="metadata/input_file_for_netosis_c_storm_gns.csv")
egs.n = unique(as.character(egs[egs$class=="netosis","egids"]))


egs.n = intersect(egs.n, featureNames(eset))

egs.c = unique(as.character(egs[egs$class=="cStorm","egids"]))
egs.c = intersect(egs.c, featureNames(eset))
rm(egs)

# create datatframe
datf = exprs(eset)[c(egs.c, egs.n),]

# create the vector for outcome
outcome.vec = as.character(ss.list[[study.id]]$Outcome)
names(outcome.vec) = sampleNames(ss.list[[study.id]])
outcome.vec.nw = outcome.vec[colnames(datf)]
rm(outcome.vec)


# subset the survivor and nonsurvivor data
dat.mdn.surv = rowMeans(as.data.frame(datf[,which(outcome.vec.nw=="survivor")]))
# names(dat.mdn.surv) = paste0("Surv_", names(dat.mdn.surv))

dat.mdn.nonsurv = rowMeans(as.data.frame(datf[,which(outcome.vec.nw=="nonsurvivor")]))
# names(dat.mdn.nonsurv) = paste0("NonSurv_", names(dat.mdn.nonsurv))

dat.mdn.ctrl = rowMeans(as.data.frame(datf[,which(outcome.vec.nw=="control")]))
# names(dat.mdn.ctrl) = paste0("Control_", names(dat.mdn.ctrl))

plotdat = cbind(dat.mdn.ctrl, dat.mdn.surv, dat.mdn.nonsurv)
gsyms.c = as.character(links(org.Hs.egSYMBOL[egs.c])[,"symbol"])

gsyms.n = as.character(links(org.Hs.egSYMBOL[egs.n])[,"symbol"])

x11()
par(mfrow=c(2, 1))
gap.barplot(t(plotdat[egs.c,]), gap = c(2,6), beside=TRUE, col=c("green", "blue","red"),
        names= gsyms.c, main= paste0("Cytokine gene \nexpression -", study.id),
        ylab="Log (gene expression)",
        #ylim=c(0,2),
        las=2)
legend(x="topleft", legend=c("Control","Survivor", "Nonsurvivor"), horiz = T,
       fill=c("green", "blue", "red"), bty="n")
box()

gap.barplot(t(plotdat[egs.n,]), gap = c(2,6), beside=TRUE, col=c("green", "blue","red"),
        names=gsyms.n, main= paste0("NETosis gene \nexpression -", study.id),
        ylab="Log (gene expression)",
        #ylim=c(0,2),
        las=2)
legend(x="topleft", legend=c("Control","Survivor", "Nonsurvivor"), horiz = T,
       fill=c("green", "blue", "red"), bty="n")
box()



