setwd("C:/Users/Jieming/Documents/thesis/mark_work/allele_specificity/confident set")

filename = 'interestingHets.betabinom.min6.snps.381samples.ase.merged.auto.bed'
data = read.table(filename, header=F, stringsAsFactors = F, comment.char="")

totalLimit = 10
data$V4[data$V4 > totalLimit] = totalLimit

numbins = 20
a = hist(data$V4,numbins,plot=F)
x11()
par(cex.axis=1.5, cex.lab=1.5)
plot(a$breaks[2:length(a$breaks)], a$counts/sum(a$counts), 
     col="blue", type="b", lty=2, xlab='minor allele frequency', ylab='fraction of SNVs', 
     cex.lab=2, cex.axis=2)
