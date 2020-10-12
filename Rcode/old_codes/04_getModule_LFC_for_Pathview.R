# Convert the sample LFC to a .csv
# file later into kegg map
#load("RDA_files_sample_level_LFC/Coagulation_module_sample_level_LFC.rda")
duke = plot.dat.list$ncov_DN
wu.PBMC = plot.dat.list$ncov_Wuhan_pbmc
wu.balf = plot.dat.list$ncov_Wuhan_balf
ny.epithelium = plot.dat.list$ncov_Lung_epithelium
ny.biopsy = plot.dat.list$ncov_Lung_biopsy

################################
write.csv(duke, "Duke_pbmc_immsupp_LFC.csv")
write.csv(wu.PBMC, "Wu_PBMC_immsupp_LFC.csv")
write.csv(wu.balf, "Wu_BALF_immsupp_LFC.csv")
write.csv(ny.epithelium, "NY_epithelium_immsupp_LFC.csv")
write.csv(ny.biopsy, "NY_Biopsy_immsupp_LFC.csv")

###################################
# Create individual
# sample level plots
# for duke data (Coagulation Module)
#----------------------------------
#load("D:/Dropbox/SKM_LAB/COVID19/RDA_files_sample_level_LFC/Coagulation_module_sample_level_LFC.rda")

load("D:/Dropbox/SKM_LAB/COVID19/RDA_files_sample_level_LFC/Immunosuppression_module_sample_level_LFC.rda")
duke.dat = as.data.frame(plot.dat.list$ncov_DN)

# Draw pathview plots for
# Complement coagulation pathway
# Immuno suppression pathway
# Platelet activation
# T Cell Receptor pathway
#---------------------------------
library(pathview)
for(i in 1:ncol(duke.dat)){
  id = colnames(duke.dat)[i]
  lfc.vec = duke.dat[, id]
  names(lfc.vec) = rownames(duke.dat)
  ttlstr = paste0("Duke_nCOV_", id)
  pv.out <- pathview(gene.data = lfc.vec,
                     pathway.id = "04650",
                     species = "hsa",
                     out.suffix = ttlstr)
  # remove the redundant
  # .XML and .PNG files
  #unlink(c("hsa04610.png", "hsa004612.xml"))
}
#--------------------------------
# Create and move the files
# according to their case number
#-------------------------------
dir.create("Case_1_NK_Cell_Duke")
dir.create("Case_2_NK_Cell_Duke")
dir.create("Case_3_NK_Cell_Duke")
#----------------
# move the files
#-----------------
file.copy(from = list.files(pattern = "_Case1"), to = "Case_1_NK_Cell_Duke/")
#----------------------------
# Create the gif fro case-1
#----------------------------
shell("convert -set delay 80 -loop 0 *.png Duke_Case1_hsa04650_covid.gif")


file.copy(from = list.files(pattern = "_Case2"), to = "Case_2_NK_Cell_Duke/")
shell("convert -set delay 80 -loop 0 *.png Duke_Case2_hsa04650_covid.gif")


file.copy(from = list.files(pattern = "_Case3"), to = "Case_3_NK_Cell_Duke/")
shell("convert -set delay 80 -loop 0 *.png Duke_Case3_hsa04650_covid.gif")

