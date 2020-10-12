# code for PCA for severity associated PCA
# in Duke NUS data
# samanwoy
# "Sat Jul 11 09:34:55 2020"
#------------------------------
#### Get Duke-NUS data ####
source("Rcode/02_get_DukeNUS_data.R")

#### Convert the data to DEseq object ######
library(ggplot2)
library(DESeq2)
# get the targets file
targets <- pData(esetDN)
gexp <- 2^exprs(esetDN)

Group = targets$Group
names(Group) = rownames(targets)
dds <- DESeqDataSetFromMatrix(countData = round(gexp),
                              colData = targets,
                              design = ~ Group)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]

## Run Differential expression analysis
dds <- DESeq(dds)
de.MIT <- results(dds)

# Get the normalised count matrix from DESeq2
# normcount <- counts(dds, normalized = TRUE)
# normcount.logged = log2(normcount+1)

# create the expression set
# esetMIT = ExpressionSet(assayData = as.matrix(normcount.logged),
#                         phenoData = AnnotatedDataFrame(targets),
#                         annotation = "org.Hs.eg.db")
# rm("dat.filtered", "Group", "normcount", "normcount.logged", "targets")


#### variance stabilisation #####
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)

#### Sample distance heatmap ####
sampleDists <- dist(t(assay(rld)))

library("pheatmap")
library("RColorBrewer")

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(rld$colnames.dat., sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         main = c("Heatmap of sample-to-sample distances using the variance \nstabilizing transformed values (Method: rlog)"),
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# keep only cases remove healthy controls
rld.case = rld[,c(which(rld$Group!="Zcontrol"))]

pcaData <- plotPCA(rld.case, intgroup = c( "ID", "Day"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = ID, shape = factor(pcaData$Day))) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")






