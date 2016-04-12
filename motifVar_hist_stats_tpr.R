setwd('C:/Users/JM/thesis/lynne_work/motifVar_pNets/tpr/align')

## data
# num of TPR proteins in each chromosome (not num of motifs)
chr      = read.table('tpr.chr.byprot.txt', header=T, stringsAsFactors = F)
motifsize= read.table('tpr.motifsize.txt', header=T, stringsAsFactors = F)
motifnum = read.table('tpr.motifnum.txt', header=T, stringsAsFactors = F)

motifnum.34aa = read.table('tpr.motifnum.34aa.txt', header=T, stringsAsFactors = F)

## plot
x11()
par(cex.axis=1.5, cex.lab=2, cex.main=2)
ff = barplot(chr$protein,names.arg=chr$chr,main='tpr proteins by chromosomes')
lines(x=ff,y=chr$gene,type="b",pch=16,col='red')

x11()
par(cex.axis=2, cex.lab=2, cex.main=2)
barplot(motifsize$occurrences,names.arg=motifsize$value,main='tpr motifsize')

x11()
par(cex.axis=2, cex.lab=2, cex.main=2)
barplot(motifnum$occurrences,names.arg=motifnum$value,main='tpr motifnum')

x11()
par(cex.axis=2, cex.lab=2, cex.main=2)
barplot(motifnum.34aa$all,names.arg=motifnum.34aa$value,main='tpr motifnum 34aa',ylim=c(0,100))