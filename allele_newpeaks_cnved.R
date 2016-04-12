setwd("C:/Users/Jieming/Documents/thesis/mark_work/allele_specificity/newpeaks_cnved")
library(VGAM)
### data
filename1 = "asb.interestingHets.betabinom.min6.stats.txt"
data1 = read.table(filename1, header=T, stringsAsFactors=F)

x11()
par(cex.axis=1,cex.lab=1.4,mar=c(8,5,1,1))
barplot(rbind(data1$bbinom.min6.oldpeaks,data1$bbinom.min6.newpeaks,data1$bbinom.min6.newpeaks.cnv),
        names.arg=data1$trio_asb,beside=T,col=c("blue","lightblue","navy"), las=2,
        xlab="", ylab="number of SNVs")

legend(100,2500, c("bbinom.min6.oldpeaks","bbinom.min6.newpeaks","bbinom.min6.newpeaks.cnv"),
       col=c("blue","lightblue","navy"), cex=1.5, pch=15)


###############################################################

filename2 = "asb.accN.min6.peaks.cnv.stats.txt"
data2 = read.table(filename2, header=T, stringsAsFactors=F)

x11()
par(cex.axis=1,cex.lab=1.4,mar=c(8,5,1,1))
barplot(rbind(data2$counts.min101.oldpeaks,data2$counts.min101.newpeaks.cnved),
        names.arg=data2$trio_asb,beside=T,col=c("blue","lightblue"), las=2,
        xlab="", ylab="number of SNVs")

legend(50,50000, c("counts.minN.oldpeaks","counts.minN.newpeaks.cnved"),
       col=c("blue","lightblue"), cex=1.5, pch=15)

###############################################################

filename3 = "ase.interestingHets.betabinom.min6.stats.txt"
data3 = read.table(filename3, header=T, stringsAsFactors=F)

x11()
par(cex.axis=1,cex.lab=1.4,mar=c(8,5,1,1))
barplot(rbind(data3$bbinom.min6.old,data3$bbinom.min6.new.cnv),border=NA,
        names.arg=data3$trio_ase,beside=T,col=c("blue","lightblue"), las=2,
        xlab="", ylab="number of SNVs")

legend(50,12000, c("bbinom.min6.old","bbinom.min6.new.cnv"),
       col=c("blue","lightblue"), cex=1.5, pch=15)

###############################################################

filename4 = "ase.accN.txt"
data4 = read.table(filename4, header=T, stringsAsFactors=F)

x11()
par(cex.axis=1,cex.lab=1.4,mar=c(8,5,1,1))
barplot(rbind(data4$counts.min10.old,data4$counts.min10.cnved),border=NA,
        names.arg=data3$trio_ase,beside=T,col=c("blue","lightblue"), las=2,
        xlab="", ylab="number of SNVs")

legend(50,250000, c("counts.minN.old","counts.minN.cnved"),
       col=c("blue","lightblue"), cex=1.5, pch=15)
