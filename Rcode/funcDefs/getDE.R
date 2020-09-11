# This is a code for getting differentiall
# expression results from an expression set
# object by row wise ttests
# (rowttests function from genefilter package)
# author: samanwoy
# "Mon Jun 08 07:26:28 2020"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create a function for DE analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function Inputs:
#~~~~~~~~~~~~~~~~~
# eset  =  an expression set object with the contrast as a FACTOR named "Group"
# direction of fold change depends on the alphabetical order
# always make control as "Zcontrol or likewise for other ariables
# FDR_cutoff = the threshold for "Multiple testing p value"
# sig.gns.only = if "TRUE" returns the result matrix only for FDR corrected significant genes
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getDE = function(eset = eset, FDR_cutoff = 0.05, sig.gns.only = TRUE){
  # load libraries
  ################
  library("genefilter")
  library("org.Hs.eg.db")

  if(sig.gns.only == "TRUE"){
    # Perform Differential expression
    rtt = rowttests(eset, fac = factor(as.character(eset$Group)))
    rownames(rtt) = featureNames(eset)
    # perform correction for multiple testing
    #########################################
    sig.gns.fdr = rownames(rtt)[p.adjust(rtt$p.value, method = "BH")< FDR_cutoff]
    sig.rtt = rtt[rownames(rtt)%in%sig.gns.fdr,]
    # return rtt mat for FDR significant genes
    ##########################################
    gsyms = as.character(unlist(mget(sig.gns.fdr, org.Hs.egSYMBOL, ifnotfound = NA)))
    res.dat = data.frame("ENTREZ_ID"= sig.gns.fdr, "Symbol" = gsyms, "LFC"= sig.rtt$dm, "pvalue_raw" = sig.rtt$p.value, "pvalue_FDR" = p.adjust(sig.rtt$p.value, method = "BH"))
    return(res.dat)
  }else{
    ##########################
    rtt = rowttests(eset, fac = factor(as.character(eset$Group)))
    rownames(rtt) = featureNames(eset)
    # return rtt mat for ALL significant genes
    ##########################################
    gsyms = as.character(unlist(mget(rownames(rtt), org.Hs.egSYMBOL, ifnotfound = NA)))
    res.dat = data.frame("ENTREZ_ID"= rownames(rtt), "Symbol" = gsyms, "LFC"= rtt$dm, "pvalue_raw" = rtt$p.value, "pvalue_FDR" = p.adjust(rtt$p.value, method = "BH"))
    return(res.dat)
  }
}