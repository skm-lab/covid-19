# Get sepsis data from two existing data packages
# These are normalized data with log(Intensity) extracted from NCBI GEO
# GEO Series ID is retained for each data set
# samanwoy
# "Fri Jun 05 19:14:45 2020"
# ssgeosurv package: "GSE26378", "GSE26440", "GSE4607",  "GSE9692",  "GSE48080", "GSE54514", "GSE63042", "GSE95233"
###########################
library("ssgeosurv")
data("ss.list")

# ssnibmg package: "GSE13904", "GSE26378", "GSE26440", "GSE4607",  "GSE8121",  "GSE9692"
###################
library("ssnibmg")
data("ss.eset")

# find the data sets overlpping between two packages
cmn.studies = intersect(names(ss.list), names(ss.eset))
eset.lst.2 = ss.eset[setdiff(names(ss.eset),cmn.studies)]

# Create a new list with 10 studies
####################################
ss.list.10 = c(ss.list, eset.lst.2)

# remove unwanted variables
rm("cmn.studies", "eset.lst.2", "ss.eset","ss.list")


