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
# pathid = "coagulation"
# pathid = "immunosuppression"
 pathid = "inflammation"
# Load the libraries

#---------------------
library(Biobase)
library(gplots)
#library(ggplot2)
#library(ggpubr)
library(org.Hs.eg.db)

# Get Duke-NUS data
source("Rcode/02_get_DukeNUS_data.R")

eset.sub1 = esetDN[featureNames(esetDN)%in%net.list[[pathid]], ]


require("igraph")
# Get the gene symbols
egs = featureNames(eset.sub1)
gsyms = mget(egs, org.Hs.egSYMBOL, ifnotfound=NA)
egs = egs[!is.na(gsyms)]
gsyms = as.character(unlist(mget(egs, org.Hs.egSYMBOL, ifnotfound=NA)))

eset.sub2 = eset.sub1[egs,]
featureNames(eset.sub2) = gsyms

xcon = exprs(eset.sub2[,eset.sub2$Group=="Zcontrol"])
xcase  = exprs(eset.sub2[,eset.sub2$Group=="nCOV"])

# create combined names for nodes
comb.nm = paste0(egs, "-", gsyms)

rownames(xcon) = comb.nm
rownames(xcase) = comb.nm

# get correlation coefficient matrix
con.Cor.Mat = cor(t(xcon))
case.Cor.Mat = cor(t(xcase))

corList = list()
corList[["con"]] = con.Cor.Mat
corList[["case"]] = case.Cor.Mat

# save file
for (i in 1: 2) {
  xlsx::write.xlsx(corList[i], file= paste0(pathid, "_Cormat.xlsx"), sheetName = names(corList)[i],
             col.names=F, row.names= T, append = T)
}


# Corr mat to graph
cor_g <- graph_from_adjacency_matrix(con.Cor.Mat, mode='undirected', weighted = "correlation", diag = F)

g = simplify(cor_g)

# Get the edge lists for which correletaion is > 0.7
cor_edge_list <- as_data_frame(cor_g, 'edges')
only_sig <- cor_edge_list[abs(cor_edge_list$correlation) > 0.7, ]
new_g <- graph_from_data_frame(only_sig, F)


only_sig
# from     to correlation
# 76 CD59  MASP1   0.9089286
# 78  CR1    CR2   0.7405917
# 79  CR1 FCGR2A   0.7849923


#some example
g2 = simplify(cor_g)
g = (cor_g)
E(g2)$weight = sapply(E(g2), function(e) {
  length(all_shortest_paths(g, from=ends(g2, e)[1], to=ends(g2, e)[2])$res) } )
plot(g2,
     vertex.size=15,
     vertex.label.dist=0.5,
     vertex.label.cex=0.8,
     vertex.label.degree=-pi/2,
     edge.width=E(g2)$weight,
     layout=layout_with_kk,
     margin=-0.2)



# subtract case matrix from control mat
diff.mat = abs(case.Cor.Mat)- abs(con.Cor.Mat)
# convert those abs diff >0.9 to 1, rest 0
diff.mat[diff.mat > 0.7] = "diff_coex_nodes"




# # select genes with absolute difference of correlation coefficients > average of all diffs
rcon = c()
rss = c()
gpairs = c()
icounter = 1
for(i in 1:(length(gsyms)-1)) {
  for(j in (i+1):length(gsyms)) {
    print(i)
    g1 = gsyms[i]
    g2 = gsyms[j]
    gp = paste0(g1,"_",g2)
    rcon[icounter] = cor(xcon[i,],xcon[j,])
    rss[icounter] = cor(xss[i,],xss[j,])
    gpairs[icounter] = gp
    icounter = icounter+1
  }
}
names(rcon) = gpairs
names(rss) = gpairs
rm(g1, g2, gp, icounter, i, j)

# select genes with absolute correlation coefficients high (>0.8) in one group
# and low (<0.3) in the other group
is.diffcoexp = apply(cbind(rcon, rss), 1, function(x) sum(abs(x)>0.3)==2 & sum(sign(x))==0)
selgp = names(which(is.diffcoexp))


###################################################
# Draw the gene pair network
###################################################
edges = data.frame("from"=sapply(strsplit(selgp, split="_"), "[[", 1),
                   "to"=sapply(strsplit(selgp, split="_"), "[[", 2))

nodes = unique(as.character(unlist(edges[,1:2])))
g <- graph_from_data_frame(d=edges, vertices=nodes, directed=F)
g <- simplify(g)
wtc <- cluster_louvain(g)
plot(g)
