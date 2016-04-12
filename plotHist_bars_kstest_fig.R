setwd("C:/Users/JM/thesis/mark_work/FIG/allelicSNPs")
library(sm)

## ks test
chipseqas     = read.table("noCDS.no0bin-AAedited.uniq.interestingHets-chipseq-peaks-joels-1KGp1-hg19-auto.all.daf.bed",sep="\t")
chipseqnonas  = read.table("noCDS.no0bin-AAedited.uniq.noninterestingHets-chipseq-peaks-joels-1KGp1-hg19-auto.all.daf.bed",sep="\t")
chipseqnonpeak= read.table("noCDS.no0bin-AAedited.uniq.nonpeaks-chipseq-joels-1KGp1-hg19-auto.all.daf.bed", sep="\t")

ks.test(chipseqas$V4,chipseqnonas$V4,alternative="less",exact=TRUE)
x11()
a <- hist(chipseqas$V4,freq=FALSE); # instead of hist use kernel density plot
b <- hist(chipseqnonas$V4,freq=FALSE);

plot(b, col=rgb(0,0,1,1/4)) # blue nonAS
plot(a, col=rgb(1,0,0,1/4), add=T) # red

rnaseqas      = read.table("noCDS.no0bin-AAedited-uniq.interestingHets-rnaseq-peaks-arifs-1KGp1-hg19.all.daf.bed",sep="\t")
rnaseqnonas   = read.table("noCDS.no0bin-AAedited-uniq.noninterestingHets-rnaseq-peaks-arifs-1KGp1-hg19.all.daf.bed",sep="\t")
rnaseqnonpeak = read.table("noCDS.no0bin-AAedited-uniq.nonpeaks-rnaseq-peaks-arifs-1KGp1-hg19.all.daf.bed",sep="\t")

ks.test(rnaseqas$V4,rnaseqnonas$V4,alternative="less",exact=TRUE)
x11()
c <- hist(rnaseqas$V4,prob=FALSE); # instead of hist use kernel density plot
d <- hist(rnaseqnonas$V4,prob=FALSE);

plot(d, col=rgb(0,0,1,1/4)) # blue nonAS
plot(c, col=rgb(1,0,0,1/4), add=T) # red
