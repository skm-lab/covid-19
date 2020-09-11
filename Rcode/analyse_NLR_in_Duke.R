# code to analyse the NLR ratio from
# Duke COVID-19 cohort
# of Cybersort result
# Samanwoy "Mon Aug 24 11:04:55 2020"
#----------------------
library(openxlsx)

# load the data
dat = read.xlsx("metadata/Input_data_for_correlation_analysis_CIBERSORTx_Job9_Results_duke_data_SKM.xlsx",
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
#par(mar=c(8,6,2,2))
b = barplot(nlr.vec,
        plot = plotLogical,
        ylim=c(0,40),
        col = c("blue", "red")[factor(dat[which.case,]$Group)],
        las=2,
        #main = "Neutrophil to Leukocyte Ratio in Duke cohort",
        #ylab = "Relative fraction of \nNLR obtained by Deconvolution",
        ylab="Neutrophil to \nLeukocyte Ratio (NLR)",
        names.arg = paste0(dat[which.case,]$PTID,
                           "_",
                           dat[which.case,]$Time))
#legend(x="topright", legend=c("Moderate Illness","Severe Illness"),       fill=c("blue","red"), bty="n")
#box()
rm(list=c("dat", "nlr.vec", "which.case","which.lymph"))

