setwd("C:/Users/Jieming/Documents/thesis/lynne_work/motifVar_pNets/tpr")

colors <- c("red","blue","green","orange","cyan","pink","purple",
            "brown","black","slategray1","violetred","tan","deeppink","darkgreen", 
            "orchid","darksalmon","antiquewhite3","magenta","darkblue","peru","slateblue",
            "thistle","tomato","rosybrown1","royalblue","olivedrab","olivedrab1","olivedrab2","olivedrab3","olivedrab4",
            "magenta1","magenta2","magenta3","magenta4") ##Set of 34 colours (no red) to use for all plots

filename1 = "combined-1KG.ESP6500.snps.nonmono.smartDomain2gPos.TPR.34aa.noS.txt.jm"
data1 = read.table(filename1,header=T)


## separate into each position
positions = 1:34
position_maf = lapply(positions, function(x) data1$maf[data1$resNum == x])
# position_maf = lapply(positions, function(x) data1$maf[(data1$resNum == x) & (data1$NS == "nonsynonymous")])

## histogram
binSize=100
bins=pretty(0:1,binSize)
# zoomx=c(0,0.008)
# zoomy=c(0.9,1)
zoomx=c(0,1)
zoomy=c(0,1)


x11()
par(cex.axis=1, cex.lab=2, cex.main=2, mar=c(5.1,4.7,4.1,2.1))

for (i in 1:length(positions))
{
  h=hist(position_maf[[i]],breaks=bins,plot=FALSE)
  
  if(i>1){ par(new=TRUE) }
  plot(h$mid,(h$counts/sum(h$counts)), type='b',col=colors[i], ylim=zoomy, xlim=zoomx)
  
#   if((i %% 2)==0)
#   {
#     text((h$mid - 0.0003),(max(h$counts/sum(h$counts))-0.01),labels=i, col=colors[i])
#   }
#   else
#   {
#     text((h$mid - 0.0001),(max(h$counts/sum(h$counts))-0.01),labels=i, col=colors[i])
#   }
}

# legend(0.002,1,
#        positions,
#        colors,
#        text.col = "black", pch = 15, bg = 'white')