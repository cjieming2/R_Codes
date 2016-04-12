setwd("C:/Users/JM/thesis/mark_work/FIG/allelicSNPs")

## wilcoxon rank sum test - compare 2 independent samples
chipseqas     = read.table("noCDS.no0bin-AAedited.uniq.interestingHets-chipseq-peaks-joels-1KGp1-hg19-auto.all.daf.bed",sep="\t")
chipseqnonas  = read.table("noCDS.no0bin-AAedited.uniq.noninterestingHets-chipseq-peaks-joels-1KGp1-hg19-auto.all.daf.bed",sep="\t")
chipseqnonpeak= read.table("noCDS.no0bin-AAedited.uniq.nonpeaks-chipseq-joels-1KGp1-hg19-auto.all.daf.bed", sep="\t")

wilcox.test(chipseqas$V4,chipseqnonas$V4,alternative="greater",exact=TRUE,paired=FALSE,conf.int=TRUE)
x11()
a <- hist(chipseqas$V4,freq=FALSE); # instead of hist use kernel density plot
b <- hist(chipseqnonas$V4,freq=FALSE);

## use density to show that the curves are similar except at the DAF < 0.05
plot(b$breaks[2:21],b$density, col=rgb(0,0,1,1/4), 
     xlim=c(0,1),ylim=c(0,2),
     type="l",
     xaxt="n")   #blue nonAS
axis(1, at=seq(0,1,0.05))
par(new = TRUE)
plot(a$breaks[2:21],a$density, col=rgb(1,0,0,1/4), 
     axes=FALSE, xlab='', ylab='',
     xlim=c(0,1),ylim=c(0,2),
     type="l") # red AS
abline(v=(seq(0,1,0.05)), col="lightgray", lty="dotted")
abline(h=(seq(0,3,0.1)), col="lightgray", lty="dotted")

rnaseqas      = read.table("noCDS.no0bin-AAedited-uniq.interestingHets-rnaseq-peaks-arifs-1KGp1-hg19.all.daf.bed",sep="\t")
rnaseqnonas   = read.table("noCDS.no0bin-AAedited-uniq.noninterestingHets-rnaseq-peaks-arifs-1KGp1-hg19.all.daf.bed",sep="\t")
rnaseqnonpeak = read.table("noCDS.no0bin-AAedited-uniq.nonpeaks-rnaseq-peaks-arifs-1KGp1-hg19.all.daf.bed",sep="\t")

wilcox.test(rnaseqas$V4,rnaseqnonas$V4,alternative="greater",exact=TRUE,paired=FALSE,conf.int=TRUE)
x11()
c <- hist(rnaseqas$V4,prob=FALSE); # instead of hist use kernel density plot
d <- hist(rnaseqnonas$V4,prob=FALSE);

plot(d$breaks[2:21], d$density, col=rgb(0,0,1,1/4),
     xlim=c(0,1), ylim=c(0,2),
     type="l") # blue nonAS
axis(1, at=seq(0,1,0.05))
par(new=TRUE)
plot(c$breaks[2:21], c$density, col=rgb(1,0,0,1/4),
     axes=FALSE, xlab='', ylab='',
     xlim=c(0,1), ylim=c(0,2),
     type="l") # red

median(chipseqas$V4)
median(chipseqnonas$V4)
median(rnaseqas$V4)
median(rnaseqnonas$V4)