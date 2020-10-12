# create a function for pathway boxplots
# Given a KEGG pathway Id,
# Expression set object
# This function draws a pathway
# boxplot based on ggplot2
# Samanwoy: Fri Jun 12 07:35:22 2020
#------------------------------------------

draw.boxplot = function(eset=esetDN, paths = "hsa04610", ttlstr = "Plot title"){
  egs = genes.by.pathway[[paths]]
  sel = featureNames(eset)%in% egs
  eset.sel = eset[sel,]

  # get the plot data
  #-------------------
  dat = colMeans(exprs(eset.sel))
  orig.fac = factor(as.character(eset$Group))
  fac = factor(orig.fac,levels = rev(levels(orig.fac)),ordered = TRUE)

  plot.dat = data.frame("name" = fac, "value" = dat)
  # Draw the boxplot
  #-------------------
  p = ggplot(plot.dat, aes(x=name, y=value, fill=name)) +

    geom_boxplot(show.legend = T, width=0.65, outlier.shape = NA)+

    xlab("")+
    theme_bw(base_size=12)+

    ylab("Gene expression\n (log2)")+

    theme(axis.title.y = element_text(face = "bold", colour = c("#619CFF")))+

    scale_fill_manual(values = c("#00BA38", "#F8766D"))+

    theme(aspect.ratio= .99)+

    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+

    ggtitle(ttlstr)+

    theme(plot.title = element_text(face="bold", color="black"))+

    theme(legend.title=element_blank()) +

    theme(legend.key.size = unit(2.95,"line"))+

    theme(legend.text=element_text(size=14))

  return(p)}