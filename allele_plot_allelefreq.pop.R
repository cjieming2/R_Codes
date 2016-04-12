setwd("C:/Users/JM/thesis/mark_work/allele_specificity/allelefreq")

## data 382
ase.cou = read.table('counts.nonintHets.minN.merged.521777snps.382samples.agg.ase.maf.auto.bed',stringsAsFactors = F)
ase.int = read.table('interestingHets.min6.merged.144083snps.382samples.maf.auto.bed',stringsAsFactors = F)

asb.cou = read.table('counts.nonint.peaks.minN.merged.239500snps.agg.382samples.asb.auto.maf.bed',stringsAsFactors = F)
asb.int = read.table('interestingHets.peaks.min6.merged.169235snps.agg.asb.382samples.auto.maf.bed',stringsAsFactors = F)

asb.maf.pop.CEU.int = read.table('asb.maf.pop.CEU.bed.382',stringsAsFactors=F)
asb.maf.pop.CHB.int = read.table('asb.maf.pop.CHB.bed.382',stringsAsFactors=F)
asb.maf.pop.JPT.int = read.table('asb.maf.pop.JPT.bed.382',stringsAsFactors=F)
asb.maf.pop.YRI.int = read.table('asb.maf.pop.YRI.bed.382',stringsAsFactors=F)

ase.maf.pop.CEU.int = read.table('ase.maf.pop.CEU.bed.382',stringsAsFactors=F)
ase.maf.pop.CHB.int = read.table('ase.maf.pop.CHB.bed.382',stringsAsFactors=F)
ase.maf.pop.FIN.int = read.table('ase.maf.pop.FIN.bed.382',stringsAsFactors=F)
ase.maf.pop.GBR.int = read.table('ase.maf.pop.GBR.bed.382',stringsAsFactors=F)
ase.maf.pop.JPT.int = read.table('ase.maf.pop.JPT.bed.382',stringsAsFactors=F)
ase.maf.pop.TSI.int = read.table('ase.maf.pop.TSI.bed.382',stringsAsFactors=F)
ase.maf.pop.YRI.int = read.table('ase.maf.pop.YRI.bed.382',stringsAsFactors=F)

asb.maf.cou.pop.CEU.int = read.table('asb.maf.cou.pop.CEU.bed',stringsAsFactors=F)
asb.maf.cou.pop.CHB.int = read.table('asb.maf.cou.pop.CHB.bed',stringsAsFactors=F)
asb.maf.cou.pop.JPT.int = read.table('asb.maf.cou.pop.JPT.bed',stringsAsFactors=F)
asb.maf.cou.pop.YRI.int = read.table('asb.maf.cou.pop.YRI.bed',stringsAsFactors=F)

ase.maf.cou.pop.CEU.int = read.table('ase.maf.cou.pop.CEU.bed',stringsAsFactors=F)
ase.maf.cou.pop.CHB.int = read.table('ase.maf.cou.pop.CHB.bed',stringsAsFactors=F)
ase.maf.cou.pop.FIN.int = read.table('ase.maf.cou.pop.FIN.bed',stringsAsFactors=F)
ase.maf.cou.pop.GBR.int = read.table('ase.maf.cou.pop.GBR.bed',stringsAsFactors=F)
ase.maf.cou.pop.JPT.int = read.table('ase.maf.cou.pop.JPT.bed',stringsAsFactors=F)
ase.maf.cou.pop.TSI.int = read.table('ase.maf.cou.pop.TSI.bed',stringsAsFactors=F)
ase.maf.cou.pop.YRI.int = read.table('ase.maf.cou.pop.YRI.bed',stringsAsFactors=F)


## hist
numbins=100
ase.hist.cou <- hist(ase.cou$V5, numbins, prob=FALSE)
ase.hist.int <- hist(ase.int$V5, numbins, prob=FALSE)
asb.hist.cou <- hist(asb.cou$V5, numbins, prob=FALSE)
asb.hist.int <- hist(asb.int$V5, numbins, prob=FALSE)

hist.asb.maf.pop.CEU.int <- hist(asb.maf.pop.CEU.int$V5, numbins, prob=FALSE)
hist.asb.maf.pop.CHB.int <- hist(asb.maf.pop.CHB.int$V5, numbins, prob=FALSE)
hist.asb.maf.pop.JPT.int <- hist(asb.maf.pop.JPT.int$V5, numbins, prob=FALSE)
hist.asb.maf.pop.YRI.int <- hist(asb.maf.pop.YRI.int$V5, numbins, prob=FALSE)

hist.ase.maf.pop.CEU.int <- hist(ase.maf.pop.CEU.int$V5, numbins, prob=FALSE)
hist.ase.maf.pop.CHB.int <- hist(ase.maf.pop.CHB.int$V5, numbins, prob=FALSE)
hist.ase.maf.pop.FIN.int <- hist(ase.maf.pop.FIN.int$V5, numbins, prob=FALSE)
hist.ase.maf.pop.GBR.int <- hist(ase.maf.pop.GBR.int$V5, numbins, prob=FALSE)
hist.ase.maf.pop.JPT.int <- hist(ase.maf.pop.JPT.int$V5, numbins, prob=FALSE)
hist.ase.maf.pop.TSI.int <- hist(ase.maf.pop.TSI.int$V5, numbins, prob=FALSE)
hist.ase.maf.pop.YRI.int <- hist(ase.maf.pop.YRI.int$V5, numbins, prob=FALSE)

hist.asb.maf.cou.pop.CEU.int <- hist(asb.maf.cou.pop.CEU.int$V5, numbins, prob=FALSE)
hist.asb.maf.cou.pop.CHB.int <- hist(asb.maf.cou.pop.CHB.int$V5, numbins, prob=FALSE)
hist.asb.maf.cou.pop.JPT.int <- hist(asb.maf.cou.pop.JPT.int$V5, numbins, prob=FALSE)
hist.asb.maf.cou.pop.YRI.int <- hist(asb.maf.cou.pop.YRI.int$V5, numbins, prob=FALSE)

hist.ase.maf.cou.pop.CEU.int <- hist(ase.maf.cou.pop.CEU.int$V5, numbins, prob=FALSE)
hist.ase.maf.cou.pop.CHB.int <- hist(ase.maf.cou.pop.CHB.int$V5, numbins, prob=FALSE)
hist.ase.maf.cou.pop.FIN.int <- hist(ase.maf.cou.pop.FIN.int$V5, numbins, prob=FALSE)
hist.ase.maf.cou.pop.GBR.int <- hist(ase.maf.cou.pop.GBR.int$V5, numbins, prob=FALSE)
hist.ase.maf.cou.pop.JPT.int <- hist(ase.maf.cou.pop.JPT.int$V5, numbins, prob=FALSE)
hist.ase.maf.cou.pop.TSI.int <- hist(ase.maf.cou.pop.TSI.int$V5, numbins, prob=FALSE)
hist.ase.maf.cou.pop.YRI.int <- hist(ase.maf.cou.pop.YRI.int$V5, numbins, prob=FALSE)




## plot col5 AF spectra ASE

xbtmlm=0
xlimit=0.5
ylimit=0.25
x11()
par(mar=c(5,8,4,4),xpd=TRUE)
# par(cex.axis=2, cex.lab=1, cex.main=1.2, cex.sub=1)
par(cex.axis=1.5, cex.lab=1.5)
plot(ase.hist.cou$breaks[2:(numbins+1)], ase.hist.cou$counts/sum(ase.hist.cou$counts), 
     col="blue", 
     type="b", lty=2, xlab='minor allele frequency', ylab='fraction of SNVs', 
     xlim=c(xbtmlm,xlimit),ylim=c(0,ylimit),
     cex.lab=2, cex.axis=2)

par(new=TRUE)
plot(ase.hist.int$breaks[2:(numbins+1)], ase.hist.int$counts/sum(ase.hist.int$counts), 
     col="blue", 
     type="b", pch=16, axes=FALSE, xlab='',ylab='',
     xlim=c(xbtmlm,xlimit),ylim=c(0,ylimit))
# plot(ase.hist.cou$breaks[2:(numbins+1)], ase.hist.cou$density, col="blue", 
#      type="b", lty=2, xlab='allele frequency', ylab='percentage', 
#      xlim=c(0,xlimit))
# par(new=TRUE)
# plot(ase.hist.int$breaks[2:(numbins+1)], ase.hist.int$density, col="blue", 
#      type="b", pch=16, axes=FALSE, xlab='',ylab='',
#      xlim=c(0,xlimit))

## ASB
par(new=TRUE)
plot(asb.hist.cou$breaks[2:(numbins+1)], asb.hist.cou$counts/sum(asb.hist.cou$counts),
     col="limegreen", 
     type="b", lty=2, axes=FALSE, xlab='',ylab='',
     xlim=c(xbtmlm,xlimit),ylim=c(0,ylimit))
par(new=TRUE)
plot(asb.hist.int$breaks[2:(numbins+1)], asb.hist.int$counts/sum(asb.hist.int$counts), 
     col="limegreen", 
     type="b", pch=16, axes=FALSE, xlab='',ylab='',
     xlim=c(xbtmlm,xlimit),ylim=c(0,ylimit))
# par(new=TRUE)
# plot(asb.hist.cou$breaks[2:(numbins+1)], asb.hist.cou$density, col="red", 
#      type="b", lty=2, axes=FALSE, xlab='',ylab='',
#      xlim=c(0,xlimit))
# par(new=TRUE)
# plot(asb.hist.int$breaks[2:(numbins+1)], asb.hist.int$density, col="red", 
#      type="b", pch=16, axes=FALSE, xlab='',ylab='',
#      xlim=c(0,xlimit))

# par(cex.axis=2, cex.lab=1, cex.main=1.2, cex.sub=1)
## pops
par(new=TRUE)
plot(hist.asb.maf.pop.CEU.int$breaks[2:(numbins+1)],
     hist.asb.maf.pop.CEU.int$counts/sum(hist.asb.maf.pop.CEU.int$counts), 
     col="red", 
     type="b", pch=16, axes=FALSE, xlab='',ylab='',
     xlim=c(xbtmlm,xlimit),ylim=c(0,ylimit))
par(new=TRUE)
plot(hist.asb.maf.cou.pop.CEU.int$breaks[2:(numbins+1)],
     hist.asb.maf.cou.pop.CEU.int$counts/sum(hist.ase.maf.cou.pop.CEU.int$counts), 
     col="red", 
     type="b", axes=FALSE, xlab='',ylab='',
     xlim=c(xbtmlm,xlimit),ylim=c(0,ylimit))

par(new=TRUE)
plot(hist.asb.maf.pop.CHB.int$breaks[2:(numbins+1)],
     hist.asb.maf.pop.CHB.int$counts/sum(hist.asb.maf.pop.CHB.int$counts), 
     col="magenta", 
     type="b", pch=16, axes=FALSE, xlab='',ylab='',
     xlim=c(xbtmlm,xlimit),ylim=c(0,ylimit))
par(new=TRUE)
plot(hist.asb.maf.cou.pop.CHB.int$breaks[2:(numbins+1)],
     hist.asb.maf.cou.pop.CHB.int$counts/sum(hist.ase.maf.cou.pop.CHB.int$counts), 
     col="magenta", 
     type="b", axes=FALSE, xlab='',ylab='',
     xlim=c(xbtmlm,xlimit),ylim=c(0,ylimit))

par(new=TRUE)
plot(hist.asb.maf.pop.JPT.int$breaks[2:(numbins+1)],
     hist.asb.maf.pop.JPT.int$counts/sum(hist.asb.maf.pop.JPT.int$counts), 
     col="magenta4", 
     type="b", pch=16, axes=FALSE, xlab='',ylab='',
     xlim=c(xbtmlm,xlimit),ylim=c(0,ylimit))
par(new=TRUE)
plot(hist.asb.maf.cou.pop.JPT.int$breaks[2:(numbins+1)],
     hist.asb.maf.cou.pop.JPT.int$counts/sum(hist.ase.maf.cou.pop.JPT.int$counts), 
     col="magenta4", 
     type="b", axes=FALSE, xlab='',ylab='',
     xlim=c(xbtmlm,xlimit),ylim=c(0,ylimit))

par(new=TRUE)
plot(hist.asb.maf.pop.YRI.int$breaks[2:(numbins+1)],
     hist.asb.maf.pop.YRI.int$counts/sum(hist.asb.maf.pop.YRI.int$counts), 
     col="black", 
     type="b", pch=16, axes=FALSE, xlab='',ylab='',
     xlim=c(xbtmlm,xlimit),ylim=c(0,ylimit))
par(new=TRUE)
plot(hist.asb.maf.cou.pop.YRI.int$breaks[2:(numbins+1)],
     hist.asb.maf.cou.pop.YRI.int$counts/sum(hist.ase.maf.cou.pop.YRI.int$counts), 
     col="black", 
     type="b", axes=FALSE, xlab='',ylab='',
     xlim=c(xbtmlm,xlimit),ylim=c(0,ylimit))

par(new=TRUE)
plot(hist.ase.maf.pop.CEU.int$breaks[2:(numbins+1)],
     hist.ase.maf.pop.CEU.int$counts/sum(hist.ase.maf.pop.CEU.int$counts), 
     col="orange", 
     type="b", pch=16, axes=FALSE, xlab='',ylab='',
     xlim=c(xbtmlm,xlimit),ylim=c(0,ylimit))
par(new=TRUE)
plot(hist.ase.maf.cou.pop.CEU.int$breaks[2:(numbins+1)],
     hist.ase.maf.cou.pop.CEU.int$counts/sum(hist.ase.maf.cou.pop.CEU.int$counts), 
     col="orange", 
     type="b", axes=FALSE, xlab='',ylab='',
     xlim=c(xbtmlm,xlimit),ylim=c(0,ylimit))

par(new=TRUE)
plot(hist.ase.maf.pop.CHB.int$breaks[2:(numbins+1)],
     hist.ase.maf.pop.CHB.int$counts/sum(hist.ase.maf.pop.CHB.int$counts), 
     col="orange3", 
     type="b", pch=16, axes=FALSE, xlab='',ylab='',
     xlim=c(xbtmlm,xlimit),ylim=c(0,ylimit))
par(new=TRUE)
plot(hist.ase.maf.cou.pop.CHB.int$breaks[2:(numbins+1)],
     hist.ase.maf.cou.pop.CHB.int$counts/sum(hist.ase.maf.cou.pop.CHB.int$counts), 
     col="orange3", 
     type="b", axes=FALSE, xlab='',ylab='',
     xlim=c(xbtmlm,xlimit),ylim=c(0,ylimit))

par(new=TRUE)
plot(hist.ase.maf.pop.JPT.int$breaks[2:(numbins+1)],
     hist.ase.maf.pop.JPT.int$counts/sum(hist.ase.maf.pop.JPT.int$counts), 
     col="orange4", 
     type="b", pch=16, axes=FALSE, xlab='',ylab='',
     xlim=c(xbtmlm,xlimit),ylim=c(0,ylimit))
par(new=TRUE)
plot(hist.ase.maf.cou.pop.JPT.int$breaks[2:(numbins+1)],
     hist.ase.maf.cou.pop.JPT.int$counts/sum(hist.ase.maf.cou.pop.JPT.int$counts), 
     col="orange4", 
     type="b", axes=FALSE, xlab='',ylab='',
     xlim=c(xbtmlm,xlimit),ylim=c(0,ylimit))

par(new=TRUE)
plot(hist.ase.maf.pop.FIN.int$breaks[2:(numbins+1)],
     hist.ase.maf.pop.FIN.int$counts/sum(hist.ase.maf.pop.FIN.int$counts), 
     col="lightsalmon", 
     type="b", pch=16, axes=FALSE, xlab='',ylab='',
     xlim=c(xbtmlm,xlimit),ylim=c(0,ylimit))
par(new=TRUE)
plot(hist.ase.maf.cou.pop.FIN.int$breaks[2:(numbins+1)],
     hist.ase.maf.cou.pop.FIN.int$counts/sum(hist.ase.maf.cou.pop.FIN.int$counts), 
     col="lightsalmon", 
     type="b", axes=FALSE, xlab='',ylab='',
     xlim=c(xbtmlm,xlimit),ylim=c(0,ylimit))


par(new=TRUE)
plot(hist.ase.maf.pop.GBR.int$breaks[2:(numbins+1)],
     hist.ase.maf.pop.GBR.int$counts/sum(hist.ase.maf.pop.GBR.int$counts), 
     col="lightsalmon3", 
     type="b", pch=16, axes=FALSE, xlab='',ylab='',
     xlim=c(xbtmlm,xlimit),ylim=c(0,ylimit))
par(new=TRUE)
plot(hist.ase.maf.cou.pop.GBR.int$breaks[2:(numbins+1)],
     hist.ase.maf.cou.pop.GBR.int$counts/sum(hist.ase.maf.cou.pop.GBR.int$counts), 
     col="lightsalmon3", 
     type="b", axes=FALSE, xlab='',ylab='',
     xlim=c(xbtmlm,xlimit),ylim=c(0,ylimit))

par(new=TRUE)
plot(hist.ase.maf.pop.TSI.int$breaks[2:(numbins+1)],
     hist.ase.maf.pop.TSI.int$counts/sum(hist.ase.maf.pop.TSI.int$counts), 
     col="lightsalmon4", 
     type="b", pch=16, axes=FALSE, xlab='',ylab='',
     xlim=c(xbtmlm,xlimit),ylim=c(0,ylimit))
par(new=TRUE)
plot(hist.ase.maf.cou.pop.TSI.int$breaks[2:(numbins+1)],
     hist.ase.maf.cou.pop.TSI.int$counts/sum(hist.ase.maf.cou.pop.TSI.int$counts), 
     col="lightsalmon4", 
     type="b", axes=FALSE, xlab='',ylab='',
     xlim=c(xbtmlm,xlimit),ylim=c(0,ylimit))


par(new=TRUE)
plot(hist.ase.maf.pop.YRI.int$breaks[2:(numbins+1)],
     hist.ase.maf.pop.YRI.int$counts/sum(hist.ase.maf.pop.YRI.int$counts), 
     col="navajowhite4", 
     type="b", pch=16, axes=FALSE, xlab='',ylab='',
     xlim=c(xbtmlm,xlimit),ylim=c(0,ylimit))
par(new=TRUE)
plot(hist.ase.maf.cou.pop.YRI.int$breaks[2:(numbins+1)],
     hist.ase.maf.cou.pop.YRI.int$counts/sum(hist.ase.maf.cou.pop.YRI.int$counts), 
     col="navajowhite4", 
     type="b", axes=FALSE, xlab='',ylab='',
     xlim=c(xbtmlm,xlimit),ylim=c(0,ylimit))



## legend
legend(0.015,0.2,
       c('ASE-','ASE+','ASB-','ASB+','ASB+_CEU','ASB+_CHB','ASB+_JPT','ASB+_YRI',
         'ASE+_CEU','ASE+_CHB','ASE+_JPT','ASE+_FIN','ASE+_GBR','ASE+_TSI','ASE+_YRI'), 
       col=c('blue','blue','limegreen','limegreen','red','magenta','magenta4','black',
          'orange','orange3','orange4','lightsalmon','lightsalmon3','lightsalmon4','navajowhite4'), 
       text.col = "black", pch = c(21,16,21,16,21,21,21,21,21,21), bg = 'white',cex=2)

