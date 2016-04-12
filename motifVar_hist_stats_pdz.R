setwd('C:/Users/Jieming/Documents/thesis/lynne_work/motifVar/pdz')

domain = 'PDZ'
## data
# num of TPR proteins in each chromosome (not num of motifs) for 30aa
chr      = read.table('PDZ.chr.byprot.txt', header=T, stringsAsFactors = F)
motifsize= read.table('PDZ.motifsize.txt', header=T, stringsAsFactors = F)
motifnum = read.table('PDZ.motifnum.txt', header=T, stringsAsFactors = F)

# motifnum.34aa = read.table('ank.motifnum.30aa.txt', header=T, stringsAsFactors = F)

## plot
x11()
par(cex.axis=1.5, cex.lab=2, cex.main=2)
ff = barplot(chr$protein,names.arg=chr$chr,main=paste(domain,' proteins by chromosomes'))
lines(x=ff,y=chr$gene,type="b",pch=16,col='red')
legend(14,40,c("by protein", "by gene"), cex=2, pch=c(15,16), col=c("grey50","red"))

x11()
par(cex.axis=2, cex.lab=1.5, cex.main=2)
barplot(motifsize$occurrences,names.arg=motifsize$value,main=paste(domain,' motifsize'), las=2)

x11()
par(cex.axis=2, cex.lab=2, cex.main=2)
barplot(motifnum$occurrences,names.arg=motifnum$value,main=paste(domain,' motifnum'))

# x11()
# par(cex.axis=2, cex.lab=2, cex.main=2)
# barplot(motifnum.34aa$occurrences,names.arg=motifnum.34aa$value,main='ank motifnum 30aa',ylim=c(0,150))