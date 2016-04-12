setwd("C:/Users/JM/thesis/mark_work/allele_specificity/enrichment/2-others-gencode")

## print ASE gene struct enrichment
# data = read.table('fishersresults_gencode-asb.txt',stringsAsFactors = F)
# data = read.table('fishersresults_gencode-ase.txt',stringsAsFactors = F)
# data = read.table('gene-asb.txt',stringsAsFactors = F)
# data = read.table('gene-asb.txt',stringsAsFactors = F)
data = read.table('enhancer-fig.txt',stringsAsFactors = F)
x11()

# p val < 0.05 (for vertical barplot)
par(mar=c(5,8,4,4),xpd=TRUE)
bp = barplot(data$V2, names.arg=data$V1, cex.axis=2, 
             cex.names=1.5, cex.lab=2,ylab="odds ratio", ylim=c(0,3))
abline(h=1, col="red", lty="dotted")
for (i in 1:length(bp))
{
  if(data$V3[i] < 0.001)
  {    text(bp[i], data$V2[i], labels = "***", pos=3, cex=2)   }
  else if(data$V3[i] < 0.01)
  {    text(bp[i], data$V2[i], labels = "**", pos=3, cex=2)    }
  else if(data$V3[i] < 0.05)
  {    text(bp[i], data$V2[i], labels = "*", pos=3, cex=2)     }
  text(bp[i], data$V2[i]+0.3, labels = round(data$V3[i],2), pos=3, cex=2)
}

## print MAE enrichment
setwd("C:/Users/JM/thesis/mark_work/allele_specificity/enrichment/mae")
x11()
# data = read.table('fishersresults_mae.txt',stringsAsFactors = F)
data = read.table('fishersresults_mae-asb.txt',stringsAsFactors = F)
## p val < 0.05 (for horiz barplot)
par(mar=c(5,12,4,6),xpd=TRUE)
bp = barplot(data$V2, names.arg=data$V1, cex.axis=2, 
             cex.names=1.5, cex.lab=2, xlab="odds ratio",
             horiz=T, las=1, xlim=c(0,10))
abline(v=1, col="red", lty="dotted")
for (i in 1:length(bp))
{
  if(data$V3[i] < 0.05)
  {
    text(data$V2[i], bp[i], labels = "*", pos=4, cex=2) 
  }
  text(data$V2[i]+0.2, bp[i], labels = round(data$V3[i],2), pos=4, cex=2)
}

## print pgene
# setwd("C:/Users/JM/thesis/mark_work/allele_specificity/enrichment/ase")
setwd("C:/Users/JM/thesis/mark_work/allele_specificity/enrichment/asb")
x11()
data = read.table('fishersresults_pgene-asb.txt',stringsAsFactors = F)
# data = read.table('fishersresults_pgene-ase.txt',stringsAsFactors = F)
## p val < 0.05 (for horiz barplot)
par(mar=c(5,15,4,6),xpd=TRUE)
bp = barplot(data$V2, names.arg=data$V1, cex.axis=2, 
             cex.names=1.5, cex.lab=2, xlab="odds ratio",
             horiz=T, las=1)
abline(v=1, col="red", lty="dotted")
for (i in 1:length(bp))
{
  if(data$V3[i] < 0.05)
  {
    text(data$V2[i], bp[i], labels = "*", pos=4, cex=2) 
  }
  text(data$V2[i]+0.2, bp[i], labels = round(data$V3[i],2), pos=4, cex=2)
}

## large plot 122 cat
# data = read.table('fishersresults_954-asb.txt',stringsAsFactors = F)
## p val < 0.05 (for horiz barplot)
# par(mar=c(5,8,4,6),oma=c(0,0,-3,1),xpd=TRUE)
# bp = barplot(data$V2, names.arg=data$V1, cex.axis=1, 
#              cex.names=1, cex.lab=1, xlab="odds ratio",
#              horiz=T, las=1)
# abline(v=1, col="red", lty="dotted")
# for (i in 1:length(bp))
# {
#   if(data$V3[i] < 0.05)
#   {
#     text(data$V2[i], bp[i], labels = "*", pos=4, cex=1) 
#   }
# #   text(data$V2[i]+0.5, bp[i], labels = round(data$V3[i],2), pos=4, cex=1)
# }

## LOF
setwd("C:/Users/JM/thesis/mark_work/allele_specificity/enrichment/lof")

data = read.table('fishersresults_lof.txt',stringsAsFactors = F)
x11()

# p val < 0.05 (for vertical barplot)
par(mar=c(5,8,4,4),xpd=TRUE)
bp = barplot(data$V2, width=c(0.2,0.2),xlim=c(0,2),names.arg=data$V1, cex.axis=2, 
             cex.names=1.5, cex.lab=2,ylab="odds ratio", ylim=c(0,2))
abline(h=1, col="red", lty="dotted")
for (i in 1:length(bp))
{
  if(data$V3[i] < 0.05)
  {
    text(bp[i], data$V2[i], labels = "*", pos=3, cex=2) 
  }
    text(bp[i], data$V2[i]+0.1, labels = round(data$V3[i],2), pos=3, cex=2)
}

## enhancers
# setwd("C:/Users/JM/thesis/mark_work/allele_specificity/enrichment/enhancers")
setwd("C:/Users/JM/thesis/mark_work/allele_specificity/enrichment/2-others-gencode")
# setwd("C:/Users/JM/thesis/mark_work/allele_specificity/enrichment/1-954")
# data = read.table('fishersresults_enhancers.txt',stringsAsFactors = F)
# data = read.table('enhancers-asb.txt',stringsAsFactors = F)
# data = read.table('enhancers-ase.txt',stringsAsFactors = F)
# data = read.table('fig-prox-promoter.txt',stringsAsFactors = F)
data = read.table('fig-genes.txt',stringsAsFactors = F)
x11()

# p val < 0.05 (for horizontal barplot)
par(mar=c(5,10,4,6),xpd=TRUE)
bp = barplot(data$V2, names.arg=data$V1, cex.axis=2, 
             cex.names=1.5, cex.lab=2, xlab="odds ratio",
             horiz=T, las=1, xlim=c(0,3))
#              col=c("green","blue","green","blue"))
abline(v=1, col="red", lty="dotted")
for (i in 1:length(bp))
{
  if(data$V3[i] < 0.001)
  {    text(data$V2[i], bp[i], labels = "***", pos=4, cex=2)   }
  else if(data$V3[i] < 0.01)
  {    text(data$V2[i], bp[i], labels = "**", pos=4, cex=2)    }
  else if(data$V3[i] < 0.05)
  {    text(data$V2[i], bp[i], labels = "*", pos=4, cex=2)     }
  text(data$V2[i]+0.15, bp[i], labels = round(data$V3[i],2), pos=4, cex=2)
}
