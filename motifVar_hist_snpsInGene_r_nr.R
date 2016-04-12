setwd("C:/Users/JM/thesis/lynne_work/motifVar_pNets/hist repeats nonrepeats")

## data
data.r = read.table('hist.1KG.snps.vat.18domains.repeats.uniq.txt', header=T, stringsAsFactors = F)
data.nr = read.table('hist.1KG.snps.vat.16domains.nonrepeats.uniq.txt', header=T, stringsAsFactors = F)


## plot 
numbins = 100
x11()
par(mar=c(5,8,4,4),xpd=TRUE)

# par(cex.axis=2, cex.lab=1, cex.main=1.2, cex.sub=1)
par(cex.axis=2, cex.lab=2,cex.main=2)
hist.numSnps.r <- hist(data.r$occurrences,numbins)

plot(hist.numSnps.r$mids, hist.numSnps.r$counts/sum(hist.numSnps.r$counts), 
     col="red", type = 'b', pch=16,main=paste('proteins with 18 repeat domains, binSize ',numbins),
     lty=1, xlab='num of SNPs',ylab='fraction of genes')

numbins = 50
x11()
par(mar=c(5,8,4,4),xpd=TRUE)
# par(cex.axis=2, cex.lab=1, cex.main=1.2, cex.sub=1)
par(cex.axis=2, cex.lab=2,cex.main=2)
hist.numSnps.nr <- hist(data.nr$occurrences,numbins)

plot(hist.numSnps.nr$mids, hist.numSnps.nr$counts/sum(hist.numSnps.nr$counts), 
     col="red", type = 'b', pch=16,main=paste('proteins with 16 non-repeat domains, binSize ',numbins),
     lty=1, xlab='num of SNPs',ylab='fraction of genes')
## legend
# legend(30,0.6,c('density','absNumSnps'), 
#        col=c('blue','red'), 
#        text.col = "black", pch=16, bg = 'white',cex=2)

###################################################################
setwd("C:/Users/JM/thesis/lynne_work/motifVar_pNets/hist repeats nonrepeats")
## snpdensity data
data.sd.r = read.table('hist.snpdensity.18domains.repeats.txt', header=T, stringsAsFactors = F, sep="\t")
data.sd.nr = read.table('hist.snpdensity.16domains.nonrepeats.txt', header=T, stringsAsFactors = F, sep="\t")

factor = 1000
sd.r = data.sd.r$numSNPs / data.sd.r$only.domainSize.in.protein * factor
data.sd.r = cbind(data.sd.r,sd.r)
sd.nr = data.sd.nr$numSNPs / data.sd.nr$only.domainSize.in.protein * factor
data.sd.nr = cbind(data.sd.nr,sd.nr)

## plot 
numbins = 50
hist.sd.r <- hist(data.sd.r$sd.r,numbins)
hist.sd.nr <- hist(data.sd.nr$sd.nr,numbins)

x11()
par(mar=c(5,8,4,4),xpd=TRUE)
# par(cex.axis=2, cex.lab=1, cex.main=1.2, cex.sub=1)
par(cex.axis=2, cex.lab=2,cex.main=2)
plot(hist.sd.r$mids, hist.sd.r$counts/sum(hist.sd.r$counts), 
     col="red", type = 'b', pch=16,main=paste('snpdensity proteins with domains, binSize ',numbins),
     lty=1, xlab=paste('snp density (',1/factor,')'),ylab='fraction of genes',
     xlim=c(0,50),ylim=c(0,0.15))

par(new=TRUE)
plot(hist.sd.nr$mids, hist.sd.nr$counts/sum(hist.sd.nr$counts), axes=F,
     col="blue", type = 'b',pch=16, xlab="", ylab="",
     xlim=c(0,50),ylim=c(0,0.15))

## legend
legend(30,0.05,c('18 repeat','16 nonrepeat'), 
       col=c('red','blue'), 
       text.col = "black", pch=16, bg = 'white',cex=2)