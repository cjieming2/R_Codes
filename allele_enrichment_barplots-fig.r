setwd("C:/Users/Jieming/Documents/thesis/mark_work/allele_specificity/enrichment/fig")

data = read.table('expanded fig v1.txt',stringsAsFactors = F)
logdata = log10(data$V2)
bonf.p = apply(as.matrix(data$V3), 1,function(x) min(x * 20052,1))
data = cbind(data,logdata,bonf.p)
x11()

# p val < 0.05 (for horizontal barplot)
par(mar=c(5,15,4,6),xpd=TRUE)
# bp = barplot(data$V2, names.arg=data$V1, cex.axis=2, 
#              cex.names=1.5, cex.lab=2, xlab="odds ratio",
#              horiz=T, las=1, xlim=c(0,3))
#              col=c("green","blue","green","blue"))
bp = barplot(logdata, names.arg=data$V1, cex.axis=2, 
             cex.names=1.5, cex.lab=2, xlab="log odds ratio",
             horiz=T, las=1, xlim=c(-0.6,0.6),
             col=c("blue","green"))
abline(v=1, col="red", lty="dotted")
for (i in 1:length(bp))
{
  if(data$bonf.p[i] < 0.001)
  {    text(logdata[i], bp[i], labels = "***", pos=4, cex=2)   }
  else if(data$bonf.p[i] < 0.01)
  {    text(logdata[i], bp[i], labels = "**", pos=4, cex=2)    }
  else if(data$bonf.p[i] < 0.05)
  {    text(logdata[i], bp[i], labels = "*", pos=4, cex=2)     }
#   text(logdata[i]+0.15, bp[i], labels = round(data$V3[i],2), pos=4, cex=2)
}
