# Extra 1
# Module-related analysis code
# Was used for checking the data
# may be removed later

############################################################
# Do a tabulation of selected KEGG pathways and 
# overlapping genes in the modules
############################################################
module = "immunosuppression"
imms.ids=c("hsa04612","hsa04660","hsa04650")
# 04612	Antigen processing and presentation
# 04660	T cell receptor signaling pathway
# 04650	Natural killer cell mediated cytotoxicity
imms.res = getTabulated(pathids=imms.ids, module)
#write.csv(imms.res, file="immunosuppression.csv")

# 04620	Toll-like receptor signaling pathway
# 04062	Chemokine signaling pathway
# 04670	Leukocyte transendothelial migration
# 04810	Regulation of actin cytoskeleton
# 04666	Fc gamma R-mediated phagocytosis
module = "inflammation"
infl.ids = c("hsa04620","hsa04062","hsa04670","hsa04810","hsa04666")
infl.res = getTabulated(pathids=infl.ids, module)
#write.csv(infl.res, file="result/inflammation.csv")

# 04610	Complement and coagulation cascades
# 04611	Platelet activation
# 04662	B-cell receptor signaling
module="coagulation"
coag.ids = c("hsa04610", "hsa04611", "hsa04662")
coag.res = getTabulated(pathids=coag.ids, module)
#write.csv(coag.res, file="result/coagulation.csv")

################################################
# Check which genes of immunosuppression modules
# are altered in covid
################################################
# Write gene annotaion to a file for marking
# the direction of gene expression manually on a paper
# keggids = c(imms.ids, infl.ids, coag.ids)
# geneData = data.frame("keggid"=NA, "gene_id"=NA,"symbol"=NA)
# for(keggid in keggids) {
#  egs = genes.by.pathway[[keggid]]
#  saveres = links(org.Hs.egSYMBOL[egs[egs%in%Lkeys(org.Hs.egSYMBOL)]])
#  saveres = data.frame(keggid, saveres)
#  geneData = rbind(geneData,saveres)  
#}
# write.table(geneData, file=paste0("result/geneMapped.csv"),
#        quote=F, sep="\t", row.names=F)
#save(geneData, file="result/geneData.rda")
# load("result/geneData.rda")
library(dplyr)
kegg = infl.ids[1]

for(kegg in infl.ids[2:5]) {
  print(pathways.list[kegg])
  egs = genes.by.pathway[[kegg]]
  egs.inflam = sapply(egs, function(eg) any(sapply(ncov, function(eset) 
    eg%in%featureNames(eset)))) %>% which %>% names
  
  require(genefilter)
  lfc.inflam = sapply(ncov, function(eset) {
    sapply(egs.inflam, function(eg) {
      if(eg%in%featureNames(eset)) {
        res=as.numeric(rowttests(eset[eg,], "Group")[,"dm"])
      } else {
        res=NA
      }
      return(res)
    })
  })
  upg.inflam = apply(lfc.inflam, 1, 
                     function(x) sum(x>0, na.rm=TRUE)==4) %>% which() %>% names()
  
  for(eg in upg.inflam) {
    drawGeneBarplot(eg);scan()}
}


# Load the three biological module genes
#----------------------------------------
# net.dat <- read.table("C:/Users/Samanwoy/Dropbox/SKM_LAB/nibmg2020/Manuscript.adult.sepsis/Analysis/metadata/network_input.csv", sep="\t", header=T)
net.dat <- read.table("C:/Users/Saroj/Dropbox/Work/SKM_LAB/nibmg2020/Manuscript.adult.sepsis/Analysis/metadata/network_input.csv", sep="\t", header=T)
net.list = with(net.dat, split(x=Egid, f=Process))
rm(net.dat)

# do some KEGG mapping with module genes
getStrings = function(egs, colour="green") {
  write.table(as.matrix(paste(egs, colour, sep=" ")), file="temp.txt", quote=F, 
              row.names=F, col.names=F)
}

getTabulated = function(pathids, module="immunosuppression") {
  egs.lst = sapply(pathids, function(id) genes.by.pathway[[id]])
  #agpp = genes.by.pathway[["hsa04612"]]
  #tcr = genes.by.pathway[["hsa04660"]]
  #nkc = genes.by.pathway[["hsa04650"]]
  egs.module = net.list[[module]]
  
  myres = sapply(egs.lst, function(egs) {
    sapply(egs.lst, function(egs1) {
      len2 = length(intersect(egs,egs1))
      len1 = length(intersect(egs.module, intersect(egs, egs1)))
      paste0(len1,"/",len2)
    })
  })
  ids = rownames(myres)
  pathnames = pathways.list[ids]
  rownames(myres) = paste0(ids, " | ", pathnames)
  return(myres)
}
