# drawGeneLinePlot
# Draw line plot for DUME-NUS COVID-19 data
# It is assumed that an expression set called esetDN
# is already available in the workspace
# Author: Samanwoy and SKM
# 02 Jul 2020

drawGeneLinePlot = function( egid = "3458"){
  gsym = as.character(unlist(mget(egid, org.Hs.egSYMBOL)))

  # Get control mean data
  ctrl.mdn = median(exprs(esetDN[egid,
                                 esetDN$Group=="Zcontrol"]))
  
  case.vec.1 = exprs(esetDN[egid,
                            esetDN$ID=="Case1"])
  day.vec.1 = as.numeric(sapply(strsplit(colnames(case.vec.1), "_d"), "[[", 2))
  
  case.vec.2 = exprs(esetDN[egid,
                            esetDN$ID=="Case2"])
  day.vec.2 = as.numeric(sapply(strsplit(colnames(case.vec.2), "_d"), "[[", 2))
  
  case.vec.3 = exprs(esetDN[egid,
                            esetDN$ID=="Case3"])
  day.vec.3 = as.numeric(sapply(strsplit(colnames(case.vec.3), "_d"), "[[", 2))
  
  #----------------
  # Case 1
  #----------------
  x <- day.vec.1
  y <- case.vec.1
  
  #------------------
  # Plot
  #------------------
  titlestr = paste0(get(egid, org.Hs.egGENENAME), " | ", 
                    get(egid, org.Hs.egSYMBOL), " | ", egid)
  par(mar=c(6.5,5,4,2))
  plot(x, y,
       type="n",
       pch=19,
       col = "red1",
       xaxt = "n",
       main= titlestr,
       xlab="Day",
       ylab="Log (gene expression)",
       ylim= range(c(case.vec.1, case.vec.2, case.vec.3)) +c(-0.2, 0.3),
       xlim=c(min(x), 20))
  title(sub="Data from Duke-NUS center (Ong 2020)", line=5)
  rect(xleft=0,
       ybottom=ctrl.mdn-0.1,
       xright= 28+max(x),
       ytop=ctrl.mdn+0.1,
       col="gray70",
       border=NA)
  abline(h= ctrl.mdn)
  lines(x=x,
        y=y,
        type="o",
        col="red1",
        pch=19,
        cex=1.2,
        lwd=3)
  xtick<-seq(0,
             28,
             by=1)
  axis(side=1,
       at=xtick,
       labels = 0:28)
  box()
  #------------------
  # Case 2
  #------------------
  x = day.vec.2
  y= case.vec.2
  lines(x=x,
        y=y,
        type="o",
        col="blue",
        pch=19,
        cex=1.2,
        lwd=3)
  #------------------
  # Case 3
  #------------------
  x = day.vec.3
  y= case.vec.3
  lines(x=x,
        y=y,
        type="o",
        col="darkgreen",
        pch=19,
        cex=1.2,
        lwd=3)
  
  #------------------
  # Add legend
  #------------------
  legend("topright",
         legend = c("Case-1 (Severe)", "Case-2 (Moderate)", "Case-3 (Moderate"),
         col = c("red1", "blue", "darkgreen"),
         bty = "n",
         lty=1,
         lwd = 4,
         cex = 1,
         horiz = F)
  
}
