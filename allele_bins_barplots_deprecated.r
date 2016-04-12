# setwd("C:/Users/JM/thesis/mark_work/allele_specificity/enrichment/bins/asb")
setwd("C:/Users/JM/thesis/mark_work/allele_specificity/enrichment/bins/ase")

int.cds  = read.table('uniq.int.cds.ase.bin',stringsAsFactors=F,sep='\t')
int.5utr = read.table('uniq.int.5utr.ase.bin',stringsAsFactors=F,sep='\t')
int.3utr = read.table('uniq.int.3utr.ase.bin',stringsAsFactors=F,sep='\t')

counts.cds  = read.table('uniq.counts.cds.ase.bin',stringsAsFactors=F,sep='\t')
counts.5utr = read.table('uniq.counts.5utr.ase.bin',stringsAsFactors=F,sep='\t')
counts.3utr = read.table('uniq.counts.3utr.ase.bin',stringsAsFactors=F,sep='\t')

## make those 11s 10s
grp11 = which(int.5utr$V5 == 11)
for (i in 1:length(grp11))
{
  int.5utr$V5[grp11[i]] = 10
}

grp11 = which(int.cds$V5 == 11)
for (i in 1:length(grp11))
{
  int.cds$V5[grp11[i]] = 10
}

grp11 = which(int.3utr$V5 == 11)
for (i in 1:length(grp11))
{
  int.3utr$V5[grp11[i]] = 10
}

## make those 11s 10s
grp11 = which(counts.5utr$V5 == 11)
for (i in 1:length(grp11))
{
  counts.5utr$V5[grp11[i]] = 10
}

grp11 = which(counts.cds$V5 == 11)
for (i in 1:length(grp11))
{
  counts.cds$V5[grp11[i]] = 10
}

grp11 = which(counts.3utr$V5 == 11)
for (i in 1:length(grp11))
{
  counts.3utr$V5[grp11[i]] = 10
}

## PLOT
table.int.cds  = table(int.cds$V5)
table.int.5utr = table(int.5utr$V5)
table.int.3utr = table(int.3utr$V5)

table.counts.cds  = table(counts.cds$V5)
table.counts.5utr = table(counts.5utr$V5)
table.counts.3utr = table(counts.3utr$V5)

x11()
plot(table.int.5utr,type='b',ylab='',main='5UTR')
par(new=TRUE)
plot(table.counts.5utr,type='b',col='red',yaxt='n',ylab='')
legend(2,500,c('int','counts'), 
       col=c('black','red'), 
       text.col = "black", pch = 20, bg = 'white')

x11()
plot(table.int.cds,type='b',ylab='',main='CDS')
par(new=TRUE)
plot(table.counts.cds,type='b',col='red',yaxt='n',ylab='')
legend(2,1000,c('int','counts'), 
       col=c('black','red'), 
       text.col = "black", pch = 20, bg = 'white')

x11()
plot(table.int.3utr,type='b',ylab='',main='3UTR')
par(new=TRUE)
plot(table.counts.3utr,type='b',col='red',yaxt='n',ylab='')
legend(2,1000,c('int','counts'), 
       col=c('black','red'), 
       text.col = "black", pch = 20, bg = 'white')
