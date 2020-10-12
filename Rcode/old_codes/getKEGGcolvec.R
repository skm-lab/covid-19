load("RDA_files_sample_level_LFC/Coagulation_module_sample_level_LFC.rda")
wu.pbmc = plot.dat.list$ncov_Wuhan_pbmc

getSampleColr = function(smplVec = wu.pbmc[,1]){
  colvec = rep("blue", length(smplVec))
  colvec[smplVec>0] = "red"
  return(colvec)
}
