setwd('C:/Users/JM/thesis/lynne_work/motifVar_pNets/bromo/stats')

## data
# num of TPR proteins in each chromosome (not num of motifs)
chr      = read.table('bromo.chr.byprot.txt', header=T, stringsAsFactors = F)
motifsize= read.table('bromo.motifsize.txt', header=T, stringsAsFactors = F)
motifnum = read.table('bromo.motifnum.txt', header=T, stringsAsFactors = F)

motifnum.34aa = read.table('bromo.109aa.motifnum.txt', header=T, stringsAsFactors = F)
chr.34aa = read.table('bromo.109aa.chr.byprot.txt', header=T, stringsAsFactors = F)

## plot
x11()
par(cex.axis=1.5, cex.lab=2, cex.main=2)
ff = barplot(chr$protein,names.arg=chr$chr,main='bromo proteins by chromosomes')
lines(x=ff,y=chr$gene,type="b",pch=16,col='red')

x11()
par(cex.axis=1.5, cex.lab=2, cex.main=2)
ff = barplot(chr.34aa$protein,names.arg=chr.34aa$chr,main='bromo proteins by chromosomes 109aa')
lines(x=ff,y=chr.34aa$gene,type="b",pch=16,col='red')


x11()
par(cex.axis=1, cex.lab=2, cex.main=2)
barplot(motifsize$occurrences,names.arg=motifsize$value,main='bromo motifsize')

x11()
par(cex.axis=2, cex.lab=2, cex.main=2)
barplot(motifnum$occurrences,names.arg=motifnum$value,main='bromo motifnum')

x11()
par(cex.axis=2, cex.lab=2, cex.main=2)
barplot(motifnum.34aa$occurrences,names.arg=motifnum.34aa$value,main='bromo motifnum 109aa',ylim=c(0,150))
