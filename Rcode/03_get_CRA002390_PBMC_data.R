# gets RNA Seq data from PBMC of COVID19 patients
# format of the data: Exprression set object
# paper: CRA002390_Xiong2020/Github
# samanwoy: "Sat Jun 06 08:56:05 2020"
######################################
eset.fn = "Data/CRA002390_PBMC.rda"
if(file.exists(file=eset.fn)) {
  cat("Reading expression set from file ...")
  load(file=eset.fn)
  cat(" done!\n")
  rm(eset.fn)
} else {
  library("org.Hs.eg.db")
  # read the raw count data
  ##########################
  dat = read.table("Data/pbmc.count.tsv", sep="\t", row.names = 1, header = T)

  #convert the ensemble ID to entrez ID
  ########################################
  ensmb.ids = sapply(strsplit(rownames(dat), ".", fixed = T), "[[", 1)
  # these ids have multiple ENTREZ IDs for single ENSEMBL IDs
  ############################################################
  ids.tobe.rmvd = names(which(sapply(mget(ensmb.ids, org.Hs.egENSEMBL2EG, ifnotfound = NA), length)>1))

  # remove those rows, repeat
  ############################
  dat.nw = dat[-c(which(ensmb.ids%in%ids.tobe.rmvd)),]
  rm(dat, ensmb.ids, ids.tobe.rmvd)
  ensmb.ids.nw = sapply(strsplit(rownames(dat.nw), ".", fixed = T), "[[", 1)
  eg.ids = as.character(unlist(mget(ensmb.ids.nw, org.Hs.egENSEMBL2EG, ifnotfound = NA)))
  rownames(dat.nw) = ensmb.ids.nw

  # keep only the rows with sucessful mapping to ENTREZID
  # also remove duplicated entrez id by variance filter
  ########################################################
  source("Rcode/varianceFilter.R")
  dat.filtered = varianceFilter(x=as.matrix(dat.nw), probeids=rownames(dat.nw), egids=eg.ids, varFunc = "var")
  rm("dat.nw", "eg.ids", "ensmb.ids.nw", "varianceFilter")

  # Normalisation with DESeq2
  ############################
  library(DESeq2)
  # get the targets file
  targets = data.frame(colnames(dat.filtered))
  rownames(targets) = colnames(dat.filtered)
  Group = c("Zcontrol", "Zcontrol", "Zcontrol", "nCOV", "nCOV", "nCOV")
  names(Group) = rownames(targets)
  targets = cbind(targets, Group)

  dds <- DESeqDataSetFromMatrix(countData =dat.filtered,
                                colData = targets,
                                design = ~ Group)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]

  ## Run Differential expression analysis
  dds <- DESeq(dds)
  de.wuhan <- results(dds)

  # Get the normalised count matrix from DESeq2
  normcount <- counts(dds, normalized = TRUE)
  normcount.logged = log2(normcount+1)

  # create the expression set
  esetWuhan = ExpressionSet(assayData = as.matrix(normcount.logged),
                       phenoData = AnnotatedDataFrame(targets),
                       annotation = "org.Hs.eg.db")
  rm("dat.filtered", "Group", "normcount", "normcount.logged", "targets")
  save(esetWuhan, dds, file=eset.fn)
  rm(eset.fn)
}


