setwd("C:/Users/JM/thesis/mark_work/allele_specificity/enrichment/")

data = read.table('ranking-asb-954+genes.txt',stringsAsFactors = F,header=T)
logdata = log(data$OR)
data = cbind(data,logdata)
x11()

# p val < 0.05 (for horizontal barplot)
par(mar=c(5,0.00000005,4,10),xpd=TRUE)

## no border
lty.o <- par("lty") 

par(lty = 0) 
# bp = barplot(logdata, names.arg=data$category, cex.axis=2, 
#              cex.names=1.5, cex.lab=2, xlab="log odds ratio",
#              horiz=T, las=1, xlim=c(-1.5,1.5),
#              col=c("blue","green"))
bp = barplot(logdata, cex.axis=2, 
             cex.names=1.5, cex.lab=2, xlab="log odds ratio",
             horiz=T, las=1, xlim=c(-1.5,1.5),
             col=data$color)

# for (i in 1:length(bp))
# {
#   if(i == 46)
#   {    text(data$logdata[i], bp[i], labels = i, pos=4, cex=1)   }
#   else if(i == 67)
#   {    text(data$logdata[i], bp[i], labels = i, pos=4, cex=1)    }
# }