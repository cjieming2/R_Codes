setwd("C:/Users/JM/thesis/lynne_work/motifVar_pNets/num SNVs in genome_proteome")

## data
# data = read.table('hist.1KG.snps.in.gencode.17.genes.uniq.txt', header=T, stringsAsFactors = F)
data = read.table('hist.1KG.snps.in.gencode.17.cds.uniq.txt', header=T, stringsAsFactors = F)

snpdensity = data$numSnps/data$size * 1000
data = cbind(data,snpdensity)

## plot 
numbins = 100
x11()
par(mar=c(5,8,4,4),xpd=TRUE)

# par(cex.axis=2, cex.lab=1, cex.main=1.2, cex.sub=1)
par(cex.axis=2, cex.lab=2,cex.main=2)
hist.density <- hist(data$snpdensity,numbins)
hist.numSnps <- hist(data$numSnps,numbins)
hist.size    <- hist(data$size,numbins)

plot(hist.density$mids, hist.density$counts/sum(hist.density$counts), 
     col="blue", type = 'b', pch=16, main='density',
     lty=1, xlab='snp density (E-3)', ylab='fraction of genes', 
     cex.lab=2, cex.axis=2)

x11()
par(mar=c(5,8,4,4),xpd=TRUE)
par(cex.axis=2, cex.lab=2,cex.main=2)
plot(hist.numSnps$mids, hist.numSnps$counts/sum(hist.numSnps$counts), 
     col="red", type = 'b', pch=16,main='counts',
     lty=1, xlab='num of SNPs',ylab='fraction of genes')

x11()
par(mar=c(5,8,4,4),xpd=TRUE)
par(cex.axis=2, cex.lab=2,cex.main=2)
plot(hist.size$mids, hist.size$counts/sum(hist.size$counts), 
     col="green", type = 'b', pch=16,main='geneSize',
     lty=1, xlab='CDS size',ylab='fraction of genes')
#      lty=1, xlab='gene size',ylab='fraction of genes')


x11()
plot(data$numSnps, data$snpdensity, 
     col="black", pch=16,main='numSnps v snpdensity',
     xlab='numSnps',ylab='snpdensity')
cr1=cor(data$numSnps,data$snpdensity,method='spearman')
text(400,60,paste("rho=",cr1))

x11()
plot(data$numSnps, data$size, 
     col="black", pch=16,main='numSnps v Size',
     xlab='numSnps',ylab='Size')
cr2=cor(data$numSnps,data$size,method='spearman')
text(400,40000,paste("rho=",cr2))

x11()
plot(data$size, data$snpdensity, 
     col="black", pch=16,main='cdsSize v snpdensity',
     xlab='cdsSize',ylab='snpdensity')
cr3=cor(data$numSnps,data$snpdensity,method='spearman')
text(40000,0.05,paste("rho=",cr3))

## legend
# legend(30,0.6,c('density','absNumSnps'), 
#        col=c('blue','red'), 
#        text.col = "black", pch=16, bg = 'white',cex=2)