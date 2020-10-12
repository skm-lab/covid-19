# get_DukeNUS_data.R
# To get the expression data from E-MTAB-8871 (from Duke-NUS center at Singapore)
# blood immune gene expression in response to SARS-COV-2
#-----------------------------------------------------------
eset.fn = "Data/EMTAB8871.rda"
if(file.exists(file=eset.fn)) {
  cat("Reading expression set from file ...")
  load(file=eset.fn)
  cat(" done!\n")
  rm(eset.fn)
} else {
  # Read meta data
  fn = "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8871/E-MTAB-8871.sdrf.txt"
  metad = read.table(file=fn, sep="\t", header=TRUE)
  Subject = as.character(metad[,"Source.Name"])
  ID = metad[,"Factor.Value.individual."]
  Day = metad[,"Factor.Value.time."]
  Group = as.character(metad[,"Factor.Value.disease."])
  Group[Group=="normal"] = "Zcontrol"
  Group[Group=="COVID-19"] = "nCOV"
  metadf = data.frame(ID, Day, Group)
  rownames(metadf) = Subject
  # remove rep2 from metadata, there is no data
  metadf = metadf[rownames(metadf)!="Case1_d4_rep2",]
  sel = rownames(metadf)=="Case1_d4_rep1"
  rownames(metadf)[sel] = "Case1_d4"
  rm(metad, Day, fn, Group, ID, Subject, sel)


  # Read data; add subject id and time to get sample id
  dat <- read.table(unzip("E-MTAB-8871.processed.1.zip"),
                    row.names=1, header=T, sep="\t")
  day = dat["Time",]
  day = sapply(day, function(x)
    ifelse(x==0, "", paste0("d",as.character(x))))

  id = sapply(colnames(dat), function(x) {
    y = strsplit(x, "\\.")[[1]]
    ifelse(is.na(y[2]), y[1], paste0(y[1],y[2]))
  })
  sample_id=paste0(id,ifelse(day=="","","_"),day)
  colnames(dat) = sample_id
  # align the rownames of metadata with colnames of data
  if(identical(sort(colnames(dat)), sort(rownames(metadf)))) {
    dat = dat[,rownames(metadf)]
  } else {
    cat("Mismatch between data and metadata sample names.\n")
  }
  metad1 = data.frame(colnames(dat), t(dat[1:2,]))
  rownames(metad1) = colnames(dat)
  metadf = data.frame(metadf, metad1)
  dat = dat[-c(1:2),]
  rm(day, id, sample_id, metad1)

  # get the entrez ids as the feature names and rownames of expression data
  fn = "https://www.ebi.ac.uk/arrayexpress/files/A-MTAB-668/A-MTAB-668.adf.txt"
  featdat = read.table(file=fn, skip=16, sep="\t", header=TRUE)
  # keep only those rows with a unique entrez id
  refseqs=as.character(featdat$Reporter.Database.Entry..refseq.)
  refseqs = sapply(refseqs, function(x) strsplit(x, split="\\.")[[1]][1])
  library("org.Hs.eg.db")
  egids = select(org.Hs.eg.db, keys=refseqs,
                 columns = "ENTREZID", keytype="REFSEQ")[,"ENTREZID"]
  featdat = data.frame(featdat, egids)
  sel = which(!duplicated(egids) & !is.na(egids))
  featdat = featdat[sel,]
  rownames(featdat) = featdat[,1]
  rm(egids, fn, refseqs, sel)

  sel = intersect(rownames(featdat), rownames(dat))
  dat = dat[sel,]
  featdat = featdat[sel,]
  if(identical(rownames(dat), rownames(featdat))) {
    rownames(dat) = featdat$egids
    rownames(featdat) = featdat$egids
  } else {
    cat("Mismatch between data rownames and feature names.\n")
  }
  rm(sel)

  # format the data; required because the file contained the character values
  # originally (2 rows); although these were removed later
  dat = apply(dat, 2, as.numeric)

  # create the expression set
  eset = ExpressionSet(assayData = as.matrix(dat),
                       phenoData = AnnotatedDataFrame(metadf),
                       featureData = AnnotatedDataFrame(featdat),
                       annotation = "org.Hs.eg.db")
  rm(dat, metadf, featdat)
  save(eset, file=eset.fn)
  rm(eset.fn)
}

# order by case and day of sampling
o = order(eset$Disease, eset$ID, eset$Day)
esetDN = eset[,o]
rm(o, eset)

# add severity information
# Case 1 - severe, Cases 2,3 - Moderate
severity = rep(NA, ncol(esetDN))
names(severity) = sampleNames(esetDN)
is.severe = esetDN$ID=="Case1"
is.mod = esetDN$ID=="Case2"|esetDN$ID=="Case3"
severity[is.severe]="Severe"
severity[is.mod]="Moderate"
esetDN$Severity=severity
rm(severity, is.mod, is.severe)


