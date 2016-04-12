library(graphics)

setwd('C:/Users/JM/thesis/mark_work/allele_specificity/datasets')
data = read.table('interestingHets.min6.merged.130127snps.agg.380samples.bed.forR',stringsAsFactors = F)
dataAcc = read.table('counts.nointhets.minN.merged.341500snps.agg.380samples.bed.forR',stringsAsFactors=F)

a=table(data$V1)
b=table(dataAcc$V1)
counts = rbind(a,b)

## snp chr distribution
x11()
par(mar=c(5,5,4,2),xpd=TRUE)
par(cex.axis=1.5, cex.lab=2, cex.main=2)
barplot(counts, main='SNV distribution by chromosomes', xlab='chromosome', ylab='frequency',
        col=c('blue','red'), beside=T, 
        names.arg=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X'))

# plot(b, ylab='frequency', xlab='chromosome', xaxt="n", col="blue", cex.axis=1.5,cex.lab=1.5)
# axis(1, at=seq(1,23),
#      labels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X'))

# par(new=T)
# 
# lines(a,col="red")
# 
legend(45, 30000, c('AS SNVs','acc SNVs'),
       col=c('blue','red'), 
       text.col = "black", pch = 15, bg = 'white', cex=1.5)


