#FromHuangetal.Feb2020Lancet
glist=c("IL1B","IL1RA","IL2","IL4","IL5","IL6","IL7","IL8","CXCL8",
        "IL9","IL10") # not found: "IL12p70"
glist=  c("IL13","IL15","IL17A") # not found "Eotaxin"
glist = c("CCL11","GCSF","CSF3","GMCSF","CSF2",
        "IFNG")
glist=c("IP-10","CXCL10","MCP1","CCL2","MIP1A","CCL3",
        "MIP1B","CCL4","PDGFB","RANTES","CCL5","TNF","VEGFA")
for(id in glist) {
  drawGeneBarplot(get(id, org.Hs.egALIAS2EG));scan()}

