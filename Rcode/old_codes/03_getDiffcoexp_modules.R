# Drawing differential coexpression network for
# 3 biological modules
# samanwoy Fri Jun 19 08:58:58 2020
#--------------------------------------
rm(list = ls())
graphics.off()
gc()
# Load the three biological module genes
#----------------------------------------
net.dat <- read.table("D:/Dropbox/SKM_LAB/nibmg2020/Manuscript.adult.sepsis/Analysis/metadata/network_input.csv", sep="\t", header=T)
net.list = with(net.dat, split(x=Egid, f=Process))

 pathid = "coagulation"
# pathid = "immunosuppression"
# pathid = "inflammation"
# Load the libraries
#---------------------
library(Biobase)
library(org.Hs.eg.db)
library(gplots)
#library(ggplot2)
#library(ggpubr)
library(org.Hs.eg.db)

# Get Duke-NUS data
source("Rcode/02_get_DukeNUS_data.R")
# subset with genes of Coagulation modules
DNmat.coag = exprs(esetDN[featureNames(esetDN)%in%net.list[[pathid]], ])

# Convert the rownames
# from ENTREZID to Gene Symbols
gsyms = as.character(unlist(mget(rownames(DNmat.coag), org.Hs.egSYMBOL)))
rownames(DNmat.coag) = gsyms


# set a vector of groups
cond.group = as.character(esetDN$Group)
names(cond.group) = colnames(esetDN)

#apply the z-score method with Spearman correlations
library(dcanr)
z_scores <- dcScore(emat = DNmat.coag, condition = cond.group, dc.method = 'diffcoex', cor.method = 'spearman')

#perform a statistical test: the z-test is selected automatically
raw_p <- dcTest(z_scores, DNmat.coag, cond.group)
print(raw_p[1:5, 1:5])

#adjust p-values (raw p-values from dcTest should NOT be modified)
adj_p <- dcAdjust(raw_p, f = p.adjust, method = 'fdr')
print(adj_p[1:5, 1:5])

# Create and Plot the Graph object
library(igraph)

#get the differential network
dcnet <- dcNetwork(z_scores, adj_p)
main.comp <- induced_subgraph(
  dcnet, V(dcnet)[components(dcnet)$membership == which.max(components(dcnet)$csize)]
)

# V(main.comp)$color <- "#FDAE61"
# l <-layout_nicely(main.comp)
l <- layout.fruchterman.reingold(main.comp)*0.5
deg = degree(main.comp, mode = "all")

graphics.off()
par(mar = c(0,0,2,2))
plot(main.comp,

     main = "Diff-Coex Coagulation
     gene network: Duke Cohort",

     layout=l,

     vertex.size = deg*4,
     vertex.color= "darkorange",
     vertex.frame.color= "black",
     edge.curved = 0.3,
     vertex.label = V(main.comp)$names,
     vertex.shape="circle",
     vertex.label.dist = 2,
     vertex.label.cex=1,
     vertex.label.font = 2,
     vertex.label.color = "black",

     edge.color = "gray60",
     edge.width=2,                                 # Edge width, defaults to 1
     #edge.arrow.size=2,                            # Arrow size, defaults to 1
     #edge.arrow.width=1,                           # Arrow width, defaults to 1
     edge.lty="solid")


#convert to an adjacency matrix
adjmat <- as_adj(dcnet, sparse = FALSE)
print(adjmat[1:5, 1:5])

#convert to a data.frame
edgedf <- as_data_frame(dcnet, what = 'edges')
print(head(edgedf))

