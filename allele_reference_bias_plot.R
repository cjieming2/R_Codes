setwd("C:/Users/Jieming/Documents/thesis/mark_work/allele_specificity/overdispersion_read depth/bias_plot")

### set parameters
# filename1 = "HG00096.rnaseq.counts.min6.allelicRatio.mod.auto.txt"
# filename2 = "HG00096.rnaseq.interestingHets.min6.allelicRatio.mod.auto.txt"
# filename3 = "HG00096.rnaseq.interestingHets.betabinom.min6.allelicRatio.mod.auto.txt"
filename1 = "NA12891.MYC.counts.min6.allelicRatio.mod.auto.txt"
filename2 = "NA12891.MYC.interestingHets.min6.allelicRatio.mod.auto.txt"
filename3 = "NA12891.MYC.interestingHets.betabinom.min6.allelicRatio.mod.auto.txt"
data1 = read.table(filename1, header=T, stringsAsFactors=F)
data2 = read.table(filename2, header=T, stringsAsFactors=F)
data3 = read.table(filename3, header=T, stringsAsFactors=F)

## binomial parameters
p=0.5               #null probability
minN=6             #min total num of reads (since it's left open right closed)
# maxN=max(data$total) #max total num of reads
if(max(data1$total) < 2500){ maxN=max(data1$total) }else { maxN=2500 }
apropor = length(data1$total[data1$total <= 2500]) / nrow(data1)
yuplimit=0.15
binSize=40
bins=pretty(0:1,binSize)

## empirical allelic Ratio
data.match=data1[data1$total <= maxN & data1$total >= minN, ]
h = hist(data1$allelicRatio, breaks=bins, plot=FALSE)
# plot right closed interval left open (range] for x axis
# note that you have pseudozeroes as counts in your data so thats fine
empirical = h$counts / sum(h$counts)
i = hist(data2$allelicRatio, breaks=bins, plot=FALSE)
allelic = i$counts / sum(h$counts)
j = hist(data3$allelicRatio, breaks=bins, plot=FALSE)
allelic2 = j$counts / sum(h$counts)

################################################
# pdf(paste(filename1,"-biasplot.pdf", sep=""), width=10, height=7)
par(cex.axis=1.5, cex.lab=2, cex.main=2, mar=c(5,5,5,5))
barplot(empirical, ylab='freq', xlab='allelicRatio', 
        names.arg=h$mids, ylim=c(0,0.15), main=paste("n=",minN,'-',maxN), 
        space=0, col='blue')
par(new=TRUE)
barplot((empirical-allelic2),ylim=c(0,0.15),col='red',
     bty='n',ylab='',xlab='',yaxt='n',xaxt='n',
     space=0)
par(new=TRUE)
barplot((empirical-allelic),ylim=c(0,0.15),col='grey',
        bty='n',ylab='',xlab='',yaxt='n',xaxt='n',
        space=0)

# barplot((empirical-allelic),ylim=c(0,max(empirical)),col='grey',
#         bty='n',ylab='',xlab='',yaxt='n',xaxt='n',
#         space=0)

# abline(h=(max(allelic) + 150), lty=2, col="red")
# dev.off()
# text(2,2750,"pcut_bin=0.009286")

##########################################################3
## counts

empirical = h$counts/ sum(h$counts)
i = hist(data2$allelicRatio, breaks=bins, plot=FALSE)
allelic = i$counts / sum(h$counts)
j = hist(data3$allelicRatio, breaks=bins, plot=FALSE)
allelic2 = j$counts/ sum(h$counts)
# pdf(paste(filename1,"-biasplot.pdf", sep=""), width=10, height=7)
par(cex.axis=1.5, cex.lab=2, cex.main=2, mar=c(5,5,5,5))
barplot(empirical, ylab='freq', xlab='allelicRatio', 
        names.arg=h$mids, ylim=c(0,0.15), main=paste("n=",minN,'-',maxN), 
        space=0, col='blue')
par(new=TRUE)
barplot((empirical-allelic),ylim=c(0,0.15),col='grey',
        bty='n',ylab='',xlab='',yaxt='n',xaxt='n',
        space=0)
par(new=TRUE)
barplot((empirical-allelic2),ylim=c(0,0.15),col='grey',
        bty='n',ylab='',xlab='',yaxt='n',xaxt='n',
        space=0)
