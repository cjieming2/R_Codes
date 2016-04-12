setwd("C:/Users/Jieming/Documents/thesis/mark_work/allele_specificity/manuscript/figures/fig2-overdispersion")
library(VGAM)
### set parameters
filename1 = "HG00096.1.M_111124_6_1.fastq.gz.dir.counts.min6.allelicRatio.mod.auto.txt"
filename2 = "HG00096.1.M_111124_6_1.fastq.gz.dir.interestingHets.betabinom.min6.allelicRatio.mod.auto.txt"
filename3 = "HG00096.1.M_111124_6_1.fastq.gz.dir.interestingHets.min6.allelicRatio.mod.auto.txt"
filename4 = "HG00096.1.M_111124_6_1.fastq.gz.dir.rnaseq.b_chosen.grad.txt"
# filename1 = "NA11894_ERR009150_1.fastq.gz.dir.counts.min6.allelicRatio.mod.auto.txt"
# filename2 = "NA11894_ERR009150_1.fastq.gz.interestingHets.betabinom.min6.allelicRatio.mod.auto.txt"
# filename3 = "NA11894_ERR009150_1.fastq.gz.dir.interestingHets.min6.allelicRatio.mod.auto.txt"
# filename4 = "NA11894_ERR009150_1.fastq.gz.dir.rnaseq.b_chosen.grad.txt"
data1 = read.table(filename1, header=T, stringsAsFactors=F)
data2 = read.table(filename2, header=T, stringsAsFactors=F)
data3 = read.table(filename3, header=T, stringsAsFactors=F)
data4 = read.table(filename4, header=T, stringsAsFactors=F)


## binomial parameters
p=0.5               #null probability
minN=6             #min total num of reads (since it's left open right closed)
maxNactual=max(data1$total) #max total num of reads
if(max(data1$total) < 2500){ maxN=max(data1$total) }else { maxN=2500 }
# apropor = length(data1$total[data1$total <= 2500]) / nrow(data1)
yuplimit=0.25
binSize=20
bins=pretty(0:1,binSize)

## empirical allelic Ratio
# data.match=data1[data1$total <= maxN & data1$total >= minN, ]
h = hist(data1$allelicRatio, breaks=bins, plot=FALSE)
# plot right closed interval left open (range] for x axis
# note that you have pseudozeroes as counts in your data so thats fine
empirical = h$counts / sum(h$counts)
i = hist(data2$allelicRatio, breaks=bins, plot=FALSE)
allelic = i$counts / sum(h$counts) ## betabin
j = hist(data3$allelicRatio, breaks=bins, plot=FALSE)
allelic2 = j$counts / sum(h$counts) ## bin

################################################
## empirical barplots only
x11()
#pdf(paste(filename1,"-biasplot.pdf", sep=""), width=10, height=7)
par(cex.axis=1.5, cex.lab=2, cex.main=2, mar=c(5,5,5,5))
barplot(empirical, ylab='freq', xlab='allelicRatio', 
        names.arg=h$mids, ylim=c(0,yuplimit), main=paste("n=",minN,'-',maxNactual), 
        space=0, col='grey')
par(new=TRUE)
barplot((allelic2),ylim=c(0,yuplimit),col='red',
     bty='n',ylab='',xlab='',yaxt='n',xaxt='n',
     space=0)
par(new=TRUE)
barplot((allelic),ylim=c(0,yuplimit),col='blue',
        bty='n',ylab='',xlab='',yaxt='n',xaxt='n',
        space=0)

# barplot((empirical-allelic),ylim=c(0,max(empirical)),col='grey',
#         bty='n',ylab='',xlab='',yaxt='n',xaxt='n',
#         space=0)

# abline(h=(max(allelic) + 150), lty=2, col="red")
# dev.off()
# text(2,2750,"pcut_bin=0.009286")

##########################################################3
## bar + null line plots
nulldistrib <- function(minN,maxN,p,w,binSize,yuplimit,distrib="binomial",b=0)
{
  d.combined = matrix(0,sum(seq(minN+1,maxN+1)),2)
  ptr = 1
  
  for (i in minN:maxN)
  {  
    ## doing the distribution
    k=seq(0,i)
    
    if(distrib == "binomial")
    {
      d = dbinom(k,i,p) ## binomial
    }
    else if(distrib == "betabinomial")
    {
      d = dbetabinom(k,i,p,b)
    }
    
    ## weight each with actual counts in empirical
    d.w = d*w[i,1]
    
    if(i == minN)
    {
      d.combined[ptr:length(k),1] = k/i
      d.combined[ptr:length(k),2] = d.w
      colnames(d.combined) = c('allelicRatio','wBinDist')
    }
    else
    {
      d.combined[ptr:(ptr+length(k)-1),1] = k/i
      d.combined[ptr:(ptr+length(k)-1),2] = d.w
    }
    
    ptr = ptr + length(k)
  }
  
  ## sort the d.combined distribution of all the n's
  d.combined.sorted = d.combined[ order(d.combined[,1],d.combined[,2]), ]
  
  ## bin it according to empirical distribution
  bins=pretty(0:1,binSize)
  start=0
  end=0
  
  d.combined.sorted.binned = matrix(0,length(bins)-1,2)
  
  for (z in 2:length(bins)) ##skip 0
  {
    start=bins[z-1]
    end=bins[z]
    
    row=z-1
    d.combined.sorted.binned[row,1] = (end-start)/2 +start ## equi of a $mid in hist
    
    d.combined.sorted.binned[row,2] = sum(d.combined.sorted[(d.combined.sorted[,1]<=end & 
                                                               d.combined.sorted[,1]>start),2])
    ## empirical right closed, left open ?hist; right=TRUE
    ## (range] so no double counts
    ## but zero gets excluded!!
    if(row==1)
    {
      d.combined.sorted.binned[row,2] = sum(d.combined.sorted[(d.combined.sorted[,1]<=end & 
                                                                 d.combined.sorted[,1]>=start),2])
    }
    
    ## empirical right closed, left open ?hist; right=TRUE
    ## (range] so no double counts
    ## but zero gets excluded!!
    if(row==1)
    {
      d.combined.sorted.binned[row,2] = sum(d.combined.sorted[(d.combined.sorted[,1]<=end & 
                                                                 d.combined.sorted[,1]>=start),2])
    }
    
  }
  
  ## change "counts" into density
  d.combined.sorted.binned[,2] = d.combined.sorted.binned[,2]/sum(d.combined.sorted.binned[,2])
  
  return(d.combined.sorted.binned)
}

t = as.data.frame(table(data1$total), stringsAsFactors=F)
w = matrix(0,max(data1$total),1)

for (jj in 1:nrow(t))
{
  w[as.integer(t[jj,1]),1] = t[jj,2]
}

overdispersion = data4$b.choice

x11()
par(cex.axis=1.5, cex.lab=2, cex.main=2, mar=c(5,5,5,5))
barplot(empirical, ylab='freq', xlab='allelicRatio', 
        names.arg=h$mids, ylim=c(0,yuplimit), main=paste("n=",minN,'-',maxN),
        space=0, col='grey')
d.combined.sorted.binned = nulldistrib(minN,maxN,p,w,binSize,yuplimit,distrib="binomial")
par(new=TRUE)
plot(d.combined.sorted.binned,ylim=c(0,yuplimit),pch=15,type='b',col='red',
     bty='n',ylab='',xlab='',yaxt='n',xaxt='n',yaxs="i")

e.combined.sorted.binned = nulldistrib(minN,maxN,p,w,binSize,yuplimit,distrib="betabinomial",b=overdispersion)
par(new=TRUE)
plot(e.combined.sorted.binned,ylim=c(0,yuplimit),pch=15,type='b',col='blue',
     bty='n',ylab='',xlab='',yaxt='n',xaxt='n',yaxs="i")
legend(0.53,0.25,c("empirical","binomial",
                   paste("betabinomial=",overdispersion)),
       col=c("grey","red","blue"), cex=2, pt.cex=2,
       text.col = "black", pch = 15, bg = 'white')

