# setwd("C:/Users/JM/thesis/mark_work/FIG/eQTLsnps")
setwd("C:/Users/JM/thesis/mark_work/FIG/allelicSNPs")

## data input
# filename1 = "fishersresults_enrichment-tfm-sf-yaoorder.txt"
# filename1 = "fishersresults_enrichment-tfp-sf-yaoorder.txt"
filename1 = "proportions_AS_rare_lt005_no0bins.txt";
# filename1 = "proportions_AS_rare_lt005.txt";

data1 = read.table(filename1,sep="\t");

## barplots
x11()
par(mar=c(5,15,5,5))
# barplot(data1$V2,xlab="Fold change", names.arg=data1$V1, horiz=TRUE, 
#         cex.names=1.5,
#         las=1) # set hori y axis labe
# abline(v=1,col='red',lty=2)

## for peaks
# colors = c("black","black","grey","black","black","grey","black","black","black","black","black","black","grey","grey","black");
colors = c("cyan","red","red","red","green","green","green") # no0bins
# colors = c("cyan","red","red","green","green","orange","orange") # 0bins
barplot(data1$V2,xlab="Proportion of rare SNPs (DAF<0.005)", names.arg=data1$V1, horiz=TRUE, 
        cex.names=1.5, col=colors,xlim=c(0.01,0.016),xpd=FALSE,
        las=1) # set hori y axis label
# abline(v=1,col='red',lty=2)
# xlim=c(0.01,0.016),