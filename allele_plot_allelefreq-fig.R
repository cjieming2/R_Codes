setwd("C:/Users/Jieming/Documents/thesis/mark_work/allele_specificity/allelefreq")

## data
# ase.cou = read.table('pop.maf.noninterestingHets.minN.snps.381samples.ase.merged.auto.bed',stringsAsFactors = F)
# ase.int = read.table('pop.maf.interestingHets.betabinom.min6.snps.381samples.ase.merged.auto.bed',stringsAsFactors = F)
# 
# asb.cou = read.table('pop.maf.noninterestingHets.peaks.minN.snps.381samples.asb.merged.auto.bed',stringsAsFactors = F)
# asb.int = read.table('pop.maf.interestingHets.betabinom.peaks.min6.snps.381samples.asb.merged.auto.bed',stringsAsFactors = F)

ase.cou = read.table('pop.maf.cds.spec.znonintHets.ase.acc.381samples.merged.uniq.auto.bed',stringsAsFactors = F)
ase.int = read.table('pop.maf.cds.spec.zintHets.ase.acc.381samples.merged.uniq.auto.bed',stringsAsFactors = F)

# asb.cou = read.table('tfm.pop.maf.noninterestingHets.peaks.minN.snps.381samples.asb.merged.auto.bed',stringsAsFactors = F)
# asb.int = read.table('tfm.pop.maf.interestingHets.betabinom.peaks.min6.snps.381samples.asb.merged.auto.bed',stringsAsFactors = F)
# asb.cou = read.table('tfm.spec.pop.maf.noninterestingHets.peaks.minN.snps.381samples.asb.merged.auto.bed',stringsAsFactors = F)
# asb.int = read.table('tfm.spec.pop.maf.interestingHets.betabinom.peaks.min6.snps.381samples.asb.merged.auto.bed',stringsAsFactors = F)

asb.cou = read.table('pop.maf.tfm.spec.znonintHets.asb.acc.13samples.merged.uniq.auto.bed',stringsAsFactors = F)
asb.int = read.table('pop.maf.tfm.spec.zintHets.asb.acc.13samples.merged.uniq.auto.bed',stringsAsFactors = F)
# asb.cou = read.table('pop.maf.tfm.znonintHets.asb.acc.13samples.merged.uniq.auto.bed',stringsAsFactors = F)
# asb.int = read.table('pop.maf.tfm.zintHets.asb.acc.13samples.merged.uniq.auto.bed',stringsAsFactors = F)

# asb.cou = read.table('bound.motifs.spec.pop.maf.noninterestingHets.peaks.minN.snps.381samples.asb.merged.auto.bed',stringsAsFactors = F)
# asb.int = read.table('bound.motifs.spec.pop.maf.interestingHets.betabinom.peaks.min6.snps.381samples.asb.merged.auto.bed',stringsAsFactors = F)


## plot col5 AF spectra ASE
x11()
numbins=100
xbtmlm=0
xlimit=0.5
ylimit=0.6
par(mar=c(5,8,4,4),xpd=TRUE)
ase.hist.cou <- hist(ase.cou$V5, numbins, prob=FALSE)
ase.hist.int <- hist(ase.int$V5, numbins, prob=FALSE)
asb.hist.cou <- hist(asb.cou$V5, numbins, prob=FALSE)
asb.hist.int <- hist(asb.int$V5, numbins, prob=FALSE)

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

## legend
legend(0.015,0.5,c('ASE-','ASE+','ASB-','ASB+'), 
       col=c('blue','blue','limegreen','limegreen'), 
       text.col = "black", pch = c(21,16,21,16), bg = 'white',cex=2)


## p values
## use /scratch/fas/gerstein/jc2296/alleledb/alleleseq-runs/aggregate_bySamples/allelefreq_sharing_pop
## 13120 cds.pop.maf.interestingHets.betabinom.min6.snps.381samples.ase.merged.auto.bed
## 87954 cds.pop.maf.noninterestingHets.minN.snps.381samples.ase.merged.auto.bed
## 782 bound.motifs.pop.maf.interestingHets.betabinom.peaks.min6.snps.381samples.asb.merged.auto.bed
## 5929 bound.motifs.pop.maf.noninterestingHets.peaks.minN.snps.381samples.asb.merged.auto.bed
afdata.asb <- matrix(c(sum(asb.int$V5 <=0.005),sum(asb.cou$V5 <=0.005),
                       sum(asb.int$V5 >0.005),sum(asb.cou$V5 >0.005)), 2,2,
               dimnames = list(Pathways = c("ASB+","ASB-"),SNPs = c("MAF<=0.005","MAF>0.005")))
x = fisher.test(afdata.asb,alternative="two.sided")
y = c("ASB+.rare",paste(x$estimate),(x$p.value)*2)
## multiple hypothesis *2

afdata.ase <- matrix(c(sum(ase.int$V5 <=0.005),sum(ase.cou$V5 <=0.005),
                       sum(ase.int$V5 >0.005),sum(ase.cou$V5 >0.005)), 2,2,
                 dimnames = list(Pathways = c("ASE+","ASE-"),SNPs = c("MAF<=0.005","MAF>0.005")))
x = fisher.test(afdata.ase,alternative="two.sided")
y = c("ASB+.rare",paste(x$estimate),(x$p.value)*2)
