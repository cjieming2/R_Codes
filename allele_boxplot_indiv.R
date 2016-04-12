setwd("C:/Users/JM/thesis/mark_work/allele_specificity/boxplot_indiv")

## data
# data.indiv = read.table('boxplot.txt',stringsAsFactors = F, header=T)
data.indiv = read.table('boxplot_minus2.txt',stringsAsFactors = F, header=T)



## boxplot for indiv hets
x11()
par(mar=c(5,8,4,4),xpd=TRUE)
par(cex.axis=1.5, cex.lab=1.5)
boxplot(het~popnum,data.indiv, 
        names=c('CEU','FIN','GBR','TSI','CHB','JPT','YRI'),
        col  =c('blue','blue','blue','blue','green','green','red'),
        xlab='population',ylab='number of heterozygous SNVs')

## boxplot for indiv intHets ASE
x11()
par(mar=c(5,8,4,4),xpd=TRUE)
par(cex.axis=1.5, cex.lab=1.5)
boxplot(ASE~popnum,data.indiv, 
        names=c('CEU','FIN','GBR','TSI','CHB','JPT','YRI'),
        col  =c('blue','blue','blue','blue','green','green','red'),
        xlab='population',ylab='number of ASE SNVs')

## boxplot for indiv intHets ASB
x11()
par(mar=c(5,8,4,4),xpd=TRUE)
par(cex.axis=1.5, cex.lab=1.5)
asb.data = cbind(ASB=data.indiv$ASB[which(data.indiv$ASB != 'NA')], 
                 popnum=data.indiv$popnum[which(data.indiv$ASB != 'NA')])
boxplot(ASB~popnum,asb.data, 
        names=c('CEU','CHB','JPT','YRI'),
        col  =c('blue','green','green','red'),
        xlab='population',ylab='number of ASB SNVs')
