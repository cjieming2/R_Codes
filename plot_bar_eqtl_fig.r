setwd("C:/Users/JM/thesis/mark_work/FIG/eQTLsnps")

## data input
# filename1 = "fishersresults_enrichment-tfm-sf-yaoorder.txt"
filename1 = "fishersresults_enrichment-tfp-sf-yaoorder.txt"
data1 = read.table(filename1,sep="\t");

## barplots
x11()
par(mar=c(5,10,5,5))
# barplot(data1$V2,xlab="Fold change", names.arg=data1$V1, horiz=TRUE, 
#         cex.names=1.5,
#         las=1) # set hori y axis labe
# abline(v=1,col='red',lty=2)

## for peaks
barplot(data1$V2,xlab="Fold change", names.arg=data1$V1, horiz=TRUE, 
        cex.names=1.5, col=c("black","black","grey","black","black","grey","black","black","black","black","black","black","grey","grey","black"),
        las=1) # set hori y axis label
abline(v=1,col='red',lty=2)