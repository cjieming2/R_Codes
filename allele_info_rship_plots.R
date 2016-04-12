library(graphics)
source("C:/Users/Jieming/Documents/Shared/scripts-R_perl_shell_macros/R codes/jmRlib.R")
setwd('C:/Users/Jieming/Documents/thesis/mark_work/allele_specificity/datasets_accN/z_statistics')
rnaseq = read.table('info.indiv.383samples.rnaseq.txt',header=T,stringsAsFactors = F)
chipseq  = read.table('info.indiv.383samples.chipseq.txt',header=T,stringsAsFactors=F)
rnaseq.db = read.table('info.indiv.383samples.rnaseq.alleleDB.txt',header=T,stringsAsFactors = F)
chipseq.db  = read.table('info.indiv.383samples.chipseq.alleleDB.txt',header=T,stringsAsFactors=F)

rnaseq.pooled = read.table('info.pooled.383samples.rnaseq.txt',header=T,stringsAsFactors=F)
rnaseq.pooled2 = rnaseq.pooled[(!rnaseq.pooled$sample == "NA12878" & !rnaseq.pooled$sample == "NA12891" & !rnaseq.pooled$sample == "NA12892"),]
chipseq.pooled = read.table('info.pooled.383samples.chipseq.txt',header=T,stringsAsFactors=F)

## individual
x11()
xlimit = c(0,1)
ylimit = c(0,max(rnaseq$size_numReads)+1000)
par(xpd=TRUE,cex.axis=1, cex.lab=1.5, cex.main=1.5)
plot(rnaseq$b.grad,rnaseq$size_numReads,col='red',pch=3,xlim=xlimit,ylim=ylimit,xlab='overdispersion',ylab='number of reads')
par(new=TRUE)
plot(chipseq$b.grad,chipseq$size_numReads,col='blue',pch=3,xlim=xlimit,ylim=ylimit,xlab='',ylab='',xaxt="n",yaxt="n")
par(new=TRUE)
plot(rnaseq.db$b.grad,rnaseq.db$size_numReads,col='red',pch=16,xlim=xlimit,ylim=ylimit,xlab='',ylab='',xaxt="n",yaxt="n")
par(new=TRUE)
plot(chipseq.db$b.grad,chipseq.db$size_numReads,col='blue',pch=16,xlim=xlimit,ylim=ylimit,xlab='',ylab='',xaxt="n",yaxt="n")
legend(0.7,5.4e8,c("rnaseq 987","chipseq 274","rnaseq.db 955", "chipseq.db 186"),
       col=c("red","blue","red","blue"), cex=2, pt.cex=2,
       text.col = "black", pch = c(3,3,16,16), bg = 'white')

## pooled
x11()
xlimit = c(0,max(rnaseq.pooled$intHets.betabin)+1000)
ylimit = c(0,max(rnaseq.pooled$numReads)+1000)
par(xpd=TRUE,cex.axis=1, cex.lab=1.5, cex.main=1.5)
plot(rnaseq.pooled$intHets.betabin,rnaseq.pooled$numReads,col='red',pch=3,xlim=xlimit,ylim=ylimit,xlab='num_intHets',ylab='number of reads')
par(new=TRUE)
plot(chipseq.pooled$intHets.betabin,chipseq.pooled$numReads,col='blue',pch=3,xlim=xlimit,ylim=ylimit,xlab='',ylab='',xaxt="n",yaxt="n")
legend(0,4e9,c("rnaseq 382","chipseq 83"),
       col=c("red","blue"), cex=2, pt.cex=2,
       text.col = "black", pch = c(3,3,16,16), bg = 'white')
identifyPch(rnaseq.pooled$intHets.betabin, rnaseq.pooled$numReads, n=length(rnaseq.pooled), 
            pch=3, col="red", tag=rnaseq.pooled$sample, sizeofpoint=2)
identifyPch(chipseq.pooled$intHets.betabin, chipseq.pooled$numReads, n=length(chipseq.pooled), 
            pch=3, col="blue", tag=chipseq.pooled$indiv_fastq, sizeofpoint=2)
a=cor(rnaseq.pooled$intHets.betabin,rnaseq.pooled$numReads,method="pearson")
b=cor(chipseq.pooled$intHets.betabin,chipseq.pooled$numReads,method="spearman")
c=cor(rnaseq.pooled2$intHets.betabin,rnaseq.pooled2$numReads,method="spearman")
text(4000,5e9,paste("rnaseq.spearman=",round(a,2),"\nrnaseq.nontrio.spearman=",round(c,2),
                    "\nchipseq.spearman=",round(b,2)), col='black', adj=c(0,0))


## pooled zoomed
x11()
xlimit = c(0,4500)
ylimit = c(0,2e9)
par(xpd=TRUE,cex.axis=1, cex.lab=1.5, cex.main=1.5)
plot(rnaseq.pooled$intHets.betabin,rnaseq.pooled$numReads,col='red',pch=3,xlim=xlimit,ylim=ylimit,xlab='num_intHets',ylab='number of reads')
par(new=TRUE)
plot(chipseq.pooled$intHets.betabin,chipseq.pooled$numReads,col='blue',pch=3,xlim=xlimit,ylim=ylimit,xlab='',ylab='',xaxt="n",yaxt="n")
legend(0,2e9,c("rnaseq 382","chipseq 83"),
       col=c("red","blue"), cex=2, pt.cex=2,
       text.col = "black", pch = c(3,3,16,16), bg = 'white')
identifyPch(rnaseq.pooled$intHets.betabin, rnaseq.pooled$numReads, n=length(rnaseq.pooled), 
            pch=3, col="red", tag=rnaseq.pooled$sample, sizeofpoint=2)
identifyPch(chipseq.pooled$intHets.betabin, chipseq.pooled$numReads, n=length(chipseq.pooled), 
            pch=3, col="blue", tag=chipseq.pooled$indiv_fastq, sizeofpoint=2)
