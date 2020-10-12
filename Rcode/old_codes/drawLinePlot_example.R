graphics.off()
rm(list = ls())
gc()
setwd("D:/Dropbox/SKM_LAB/COVID19")

net.dat <- read.table("D:/Dropbox/SKM_LAB/nibmg2020/Manuscript.adult.sepsis/Analysis/metadata/network_input.csv", sep="\t", header=T)
net.list = with(net.dat, split(x=Egid, f=Process))

pathid = "coagulation"
#pathid = "immunosuppression"
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
# Get the intersecting genes between study and module
cmn.gns = as.character(intersect(featureNames(esetDN), net.list[[pathid]]))


egid = "695"
gsym = as.character(unlist(mget(egid, org.Hs.egSYMBOL)))
ttlstr = paste0("Module: ", pathid, " Gene - ", gsym)


# Get control mean data
ctrl.mdn = median(exprs(esetDN[egid, esetDN$Group=="Zcontrol"]))

case.vec.1 = exprs(esetDN[egid,esetDN$ID=="Case1"])
day.vec.1 = as.numeric(sapply(strsplit(colnames(case.vec.1), "_d"), "[[", 2))

case.vec.2 = exprs(esetDN[egid,esetDN$ID=="Case2"])
day.vec.2 = as.numeric(sapply(strsplit(colnames(case.vec.2), "_d"), "[[", 2))

case.vec.3 = exprs(esetDN[egid,esetDN$ID=="Case3"])
day.vec.3 = as.numeric(sapply(strsplit(colnames(case.vec.3), "_d"), "[[", 2))




# Version 2
###############################################
#rm(list=ls())

# some data: copying approximately by looking
# at the box plot for immunosuppression
#----------------
# Case 1
#----------------
x <- day.vec.1
y <- case.vec.1

# Get the last x-value (i.e., day of return to base line)
# Use a linear model from the point of max logfc
get.slope = function(x,y, direction="down") {
  if(direction=="down") {
    which.peak = which.min(y)
  } else {
    which.peak = which.max(y)
  }
  x1=x[which.peak:length(x)]
  y1=y[which.peak:length(y)]
  slope = as.numeric(lm(y1~x1)$coefficients[2])
  return(slope)
}
#------------------
# Plot
#------------------
slope = get.slope(x,y)
plot(x, y,
     type="n",
     pch=19,
     col = "red1",
     xaxt = "n",
     main= paste0("Gene expression trajectory\n", ttlstr),
     xlab="Day",
     ylab="Level of gene expression (log2)",
     ylim=range(y)+c(-0.2, 1),
     xlim=c(min(x), 28))
rect(xleft=0,
     ybottom=ctrl.mdn-0.1,
     xright= 28+max(x),
     ytop=ctrl.mdn+0.1,
     col="gray90",
     border=NA)
abline(h= ctrl.mdn)
lines(x=x,
      y=y,
      type="o",
      col="red1",
      pch=19)
xtick<-seq(0,
           28,
           by=1)
axis(side=1,
     at=xtick,
     labels = 0:28)
box()
# calculate the y2 value
y2 = slope*(28 - x[length(x)]) + y[length(y)]
#-------------------------------------------
# Draw the 28 day mortality prediction line
#-------------------------------------------
lines(x=c(x[length(x)], 28),
      y=c(y[length(y)], y2),
      lty=2,
      col="red1",
      lwd=2)
cat(paste0("Days to base line (case 1): ", formatC(slope, width=2),"\n"))
#------------------
# Case 2
#------------------
x = day.vec.2
y= case.vec.2
slope = get.slope(x,y)
lines(x=x,
      y=y,
      type="o",
      col="blue",
      pch=19)
y2 = slope*(28 - x[length(x)]) + y[length(y)]
#-------------------------------------------
# Draw the 28 day mortality prediction line
#-------------------------------------------
lines(x=c(x[length(x)], 28),
      y=c(y[length(y)], y2),
      lty=2,
      col="blue",
      lwd=2)
cat(paste0("Days to base line (case 2): ", formatC(slope, width=2),"\n"))
#------------------
# Case 3
#------------------
x = day.vec.3
y= case.vec.3
slope = get.slope(x,y)
lines(x=x,
      y=y,
      type="o",
      col="green3",
      pch=19)

y2 = slope*(28 - x[length(x)]) + y[length(y)]
#-------------------------------------------
# Draw the 28 day mortality prediction line
#-------------------------------------------
lines(x=c(x[length(x)], 28),
      y=c(y[length(y)], y2),
      lty=2,
      col="green4",
      lwd=2)
cat(paste0("Days to base line (case 3): ",
           formatC(slope, width=2),
           "\n"))
#------------------
# Add legend
#------------------
legend("topleft",
       legend = c("Case-1", "Case-2", "Case-3"),
       col = c("red1", "blue", "green4"),
       bty = "n",
       lty=1,
       lwd = 3,
       cex = 0.57,
       horiz = F)
