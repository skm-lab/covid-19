# A single-cell atlas of the peripheral immune response in patients with severe COVID-19
# Wilke et.al.  2020, Nature Medicine

# Link: https://www.nature.com/articles/s41591-020-0944-y
# Samanwoy
# "Tue Aug 11 1:16:21 2020"
#################################
## covid_combined.nc <- readRDS(url("https://hosted-matrices-prod.s3-us-west-2.amazonaws.com/Single_cell_atlas_of_peripheral_immune_response_to_SARS_CoV_2_infection-25/blish_covid.seu.rds"))

## covid_combined.nc = readRDS("D:/not_junk#/Wilke_Covid_nature_medicine/blish_covid.seu.rds")

#################################################
# start from heare for analysis with bioconductor
# SingleCellClass
#################################################
library(ggpubr)
library(RColorBrewer)

if(file.exists("Data/Wilk_agg_sce.rda")){
  load("Data/Wilk_agg_sce.rda")
  cat("Loading data from file..")
}else{
  load("C:/Users/Samanwoy/Documents/Wilk_data/Wilke_sce_calss_object_converted_from_seurat.rda")

  # Aggregation across groups or clusters
  #-----------------------------------------
  #sce$Status = factor(sce$Status, levels=c("Healthy", "COVID"))
  scepheno = colData(sce)
  agg_sce <- aggregateAcrossCells(sce, ids=sce$orig.ident)
  agg_sce <- logNormCounts(agg_sce)
  save(agg_sce, file = "C:/Users/Samanwoy/Documents/Wilk_agg_sce.rda")
  }


# Get the genes of cytokine and NETosis
#-----------------------------------------
egs = read.csv(
  "metadata/input_file_for_netosis_c_storm_gns.csv"
)
egs.n = unique(as.character(egs[egs$class == "netosis", "egids"]))
#egs.n = intersect(egs.n, featureNames(eset))
egs.c = unique(as.character(egs[egs$class == "cStorm", "egids"]))
#egs.c = intersect(egs.c, featureNames(eset))
rm(egs)

library(org.Hs.eg.db)
gsyms.netosis = links(org.Hs.egSYMBOL[egs.n])[, "symbol"]
gsyms.netosis = intersect(rownames(agg_sce), gsyms.netosis)


gsyms.cstorm = links(org.Hs.egSYMBOL[egs.c])[, "symbol"]
gsyms.cstorm = intersect(rownames(agg_sce), gsyms.cstorm)

Group = rep("Healthy", ncol(agg_sce))
Group[which(agg_sce$Ventilated=="2")] = "NonVent"
Group[which(agg_sce$Ventilated=="3")] = "ARDS"




#---------------------------------------------
# Figure 1C: NETosis gene expression boxplots
#---------------------------------------------
par(mar=c(4,6,4,1))
ids = which(rownames(agg_sce)%in%gsyms.netosis)
plotdat = logcounts(agg_sce)[ids,]

plotdat.n = data.frame("Healthy"=rowMeans(plotdat[, which(Group=="Healthy")]),
                       "COVID19" = rowMeans(plotdat[, which(Group!="Healthy")]))

boxplot(plotdat.n, col = brewer.pal(3, "Oranges"),
        ylab = "Relative Gene Expression \n(Aggregate Log Counts)",
        boxwex = 0.32,
        main = "NETosis Gene Expression \nWilk_2020",
        names = c("Healthy", "COVID19"))
lines(apply(plotdat.n, 2, median), lty=2, lwd=2)

#-----------------------------------------
# Figure 9: Gene Level boxlot CEACAM8
#-----------------------------------------
gn = "CEACAM8"
plotdat = logcounts(agg_sce)[gn,]

# create dataframe
df = data.frame(group = Group, value = as.numeric(plotdat))
df$group = factor(df$group, levels= c("Healthy", "NonVent", "ARDS"))

# plotting
ggplot(df, aes(x=group, y=value, fill=group))+
  geom_boxplot(width = 0.47)+
  ylab("Relatvie gene exppression\n (aggregate logcounts)")+
  scale_fill_brewer(palette="Reds")+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.title.x = element_blank())+
  #theme(axis.title.y = element_blank())+
  ggtitle(paste0(gn ," in PBMC Wilk_2020\n"))+
  stat_compare_means(method = "anova", vjust = 01.6, label.x = 1.5)

#-----------------------------------------
# Figure 10: Gene Level boxlot DNASE1
#-----------------------------------------
gn = "DNASE1"
plotdat = logcounts(agg_sce)[gn,]

# create dataframe
df = data.frame(group = Group, value = as.numeric(plotdat))
df$group = factor(df$group, levels= c("Healthy", "NonVent", "ARDS"))

# plotting
ggplot(df, aes(x=group, y=value, fill=group))+
  geom_boxplot(width = 0.47)+
  ylab("Relatvie gene exppression\n (aggregate logcounts)")+
  scale_fill_brewer(palette="Blues")+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.title.x = element_blank())+
  #theme(axis.title.y = element_blank())+
  ggtitle(paste0(gn ," in PBMC Wilk_2020\n"))+
  stat_compare_means(method = "anova", vjust = 01.6, label.x = 1.5)

# Cytokine storm boxplot
#-----------------------------------------
ids = which(rownames(agg_sce)%in%gsyms.cstorm)
plotdat = logcounts(agg_sce)[ids,]

plotdat.n = data.frame("Healthy"=rowMeans(plotdat[, which(Group=="Healthy")]),
                       "NonVent" = rowMeans(plotdat[, which(Group=="NonVent")]),
                       "ARDS" = rowMeans(plotdat[, which(Group=="ARDS")]))
rm(plotdat)
library(RColorBrewer)
boxplot(plotdat.n, col = brewer.pal(3, "Oranges"),
        ylab = "Relative Gene Expression \n(Aggregate Log Counts)",
        boxwex = 0.32,
        main = "Cytokine Storm Gene Expression \nWilke2020")
lines(apply(plotdat.n, 2, median), lty=2, lwd=2)

rm(ids, plotdat.n)

# Netosis boxplot
#-----------------------------------------
ids = which(rownames(agg_sce)%in%gsyms.netosis)
plotdat = logcounts(agg_sce)[ids,]

plotdat.n = data.frame("Healthy"=rowMeans(plotdat[, which(Group=="Healthy")]),
                       "NonVent" = rowMeans(plotdat[, which(Group=="NonVent")]),
                       "ARDS" = rowMeans(plotdat[, which(Group=="ARDS")]))
rm(plotdat)

boxplot(plotdat.n, col = brewer.pal(3, "Oranges"),
        ylab = "Relative Gene Expression \n(Aggregate Log Counts)",
        boxwex = 0.32,
        main = "NETosis Gene Expression \nWilk_2020",
        names = c("Healthy", "NonVent", "ARDS"))
lines(apply(plotdat.n, 2, median), lty=2, lwd=2)



