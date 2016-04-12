setwd('C:/Users/JM/thesis/lynne_work/motifVar_pNets/ank/align')

## data
# num of TPR proteins in each chromosome (not num of motifs) for 30aa
chr      = read.table('ank.chr.byprot.txt', header=T, stringsAsFactors = F)
motifsize= read.table('ank.motifsize.txt', header=T, stringsAsFactors = F)
motifnum = read.table('ank.motifnum.txt', header=T, stringsAsFactors = F)

motifnum.34aa = read.table('ank.motifnum.30aa.txt', header=T, stringsAsFactors = F)

## plot
x11()
par(cex.axis=1.5, cex.lab=2, cex.main=2)
ff = barplot(chr$protein,names.arg=chr$chr,main='ank proteins by chromosomes 30aa')
lines(x=ff,y=chr$gene,type="b",pch=16,col='red')

x11()
par(cex.axis=2, cex.lab=2, cex.main=2)
barplot(motifsize$occurrences,names.arg=motifsize$value,main='ank motifsize')

x11()
par(cex.axis=2, cex.lab=2, cex.main=2)
barplot(motifnum$occurrences,names.arg=motifnum$value,main='ank motifnum')

x11()
par(cex.axis=2, cex.lab=2, cex.main=2)
barplot(motifnum.34aa$occurrences,names.arg=motifnum.34aa$value,main='ank motifnum 30aa',ylim=c(0,150))