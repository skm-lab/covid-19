# code for barplot for neutrophil
# fractions in  CIBERSORTx results
# in Duke NUS data
# samanwoy
# "Wed Jul 15 11:32:10 2020"
#------------------------------
#### Get Duke-NUS Cibersort output ####
dat = read.csv("cibersort_COVID-19/Duke_Cohort/Duke_CIBERSORTx_Job9_Results.csv", row.names = 1)
dat = dat[-c(grep("HC", rownames(dat))),]

#### Draw a barplot ######
# Load ggplot2
library(ggplot2)

# create dummy data
data <- data.frame(
  name = factor(rownames(dat), levels = rownames(dat)),
  value = dat$Neutrophils,
  std.dev = sd(dat$Neutrophils)
)

colvec = c(rep("red", 9), rep("blue", 13))
# Most basic error bar
 p = ggplot(data) +
  geom_bar( aes(x=name, y=value), stat="identity", fill= colvec,
            alpha=0.7) + xlab("")+ ylab("Percentage of Neutrophils")+  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
