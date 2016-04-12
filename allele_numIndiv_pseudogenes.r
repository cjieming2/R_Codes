setwd("C:/Users/JM/thesis/mark_work/allele_specificity/enrichment/ase/pgenes")

data6 = read.table('interestingHets.min6.merged.snps.380samples.Unitary.Pgene.bed',stringsAsFactors = F)
data5 = read.table('interestingHets.min6.merged.snps.380samples.Polymorphic.Pgene.bed',stringsAsFactors = F)
data4 = read.table('interestingHets.min6.merged.snps.380samples.Transcribed.Duplicated.Pgene.bed',stringsAsFactors = F)
data3 = read.table('interestingHets.min6.merged.snps.380samples.Transcribed.Processed.Pgene.bed',stringsAsFactors = F)
data2 = read.table('interestingHets.min6.merged.snps.380samples.Duplicated.Pgene.bed',stringsAsFactors = F)
data1 = read.table('interestingHets.min6.merged.snps.380samples.Processed.Pgene.bed',stringsAsFactors = F)
data7 = read.table('interestingHets.min6.merged.snps.380samples.promoter.bed',stringsAsFactors=F)
x11()
par(mar=c(5,15,4,2),xpd=TRUE)
hp = boxplot(data1$V4,data2$V4, data3$V4,data4$V4,data5$V4,data6$V4,data7$V4,
             xlab='num of samples',
             horizontal=T,
             cex.lab=2, cex.axis=2)

mtext(c('processed.pg','dup.pg','trans.processed.pg','trans.dup.pg',
        'polymorphic.pg','unitary.pg','promoter'),
      side=2, line=1, las=1, at=c(1,2,3,4,5,6,7), cex=2)