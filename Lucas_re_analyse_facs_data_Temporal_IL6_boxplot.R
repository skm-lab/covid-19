# this is the code to analyse FACs data for Cytokine levels
# Longitudinal analyses reveal immunological misfiring in severe COVID-19
# Link to Article: https://pubmed.ncbi.nlm.nih.gov/32717743/
# Article Published: 27 July 2020
# samanwoy
#  "Tue Jul 28 10:57:21 2020"
#------------------------------------
# Figure 8: Temporal boxplot of IL6
# from Lucas et. al.2020
#------------------------------------
library("readxl")
library("ggplot2")

#### Load the data ####
#----------------------
#dat = read_xlsx("C:/Users/Samanwoy/Dropbox/SKM_LAB/COVID-19/
#Cytokine_levels_other_than_Gexp_Studies_2020/
#Data_for_Reanalysis/41586_2020_2588_MOESM3_ESM.xlsx", sheet = 1, skip = 26)

dat = read_xlsx("C:/Users/Samanwoy/Dropbox/SKM_LAB/COVID-19/Data/41586_2020_2588_MOESM3_ESM.xlsx", sheet = 1, skip = 26)

# > dat$ID[is.na(dat$`Clinical score`)]
# [1] "Pt055.3" "Pt063.5"
to.remove = which(dat$ID%in%c("Pt055.3", "Pt063.5"))
dat = dat[-c(to.remove),]
rm(to.remove)

Group = rep("Control", nrow(dat))
Group[which(dat$`Clinical score`=="1")] = "Moderate"
Group[which(dat$`Clinical score`=="2")] = "Moderate"
Group[which(dat$`Clinical score`=="3")] = "Moderate"

Group[which(dat$`Clinical score`=="4")] = "Severe"
Group[which(dat$`Clinical score`=="5")] = "Severe"


# Construct a day vector
#------------------------
day = rep("0", nrow(dat))
token.dat = do.call(rbind, strsplit(dat$ID, ".", fixed = T))
day[which(sapply(token.dat[,2], function(x){nchar(x)}) == "1")] = token.dat[,2][which(sapply(token.dat[,2], function(x){nchar(x)}) == "1")]
rm(token.dat)

#### temporal line plot for IL6 gene ####
#---------------------------------
p4 = ggplot(dat, aes(x = day,
                     y = IL6,
                     fill=Group)) +
  geom_boxplot(lwd = 0.05, width=0.99) +
  theme_bw(base_size = 14)+
  coord_cartesian(ylim = c(-2, 4))+
  ggtitle("Temporal Boxplot of IL6 (Lucas2020)")+
  scale_fill_brewer(palette="Oranges")+
  ylab("log10 Plasma IL6 \n(pg/mL)")+
  theme(legend.title = element_blank())
# Print the plot
print(p4)
