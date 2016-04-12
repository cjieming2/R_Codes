library(VGAM)

# setwd("C:/Users/JM/thesis/mark_work/allele_specificity/read depth")
setwd("/Users/Jieming/Documents/mark_work/allele_specificity/read depth")

### set parameters
title = 'NA12878_ase_polyA'
data = read.table("counts.min6.allelicRatio.mod.auto.txt", header=T, stringsAsFactors=F)
# data = read.table("test.txt", header=T, stringsAsFactors=F)

colors <- c("green","red","blue","orange","cyan","pink","purple",
            "brown","black","slategray1","violetred","tan","deeppink","darkgreen", 
            "orchid","darksalmon","antiquewhite3","magenta","darkblue","peru","slateblue",
            "thistle","tomato","rosybrown1","royalblue","olivedrab") ##Set of 27 colours to use for all plots

## binomial parameters
p=0.5               #null probability
minN=12             #min total num of reads (since it's left open right closed)
# maxN=max(data$total) #max total num of reads
maxN=12 #max total num of reads
yuplimit=0.25
binSize=40

## betabinomial parameters
#overdispersion

## empirical allelic Ratio
# jpeg(filename=paste("allele_readdepth_table_binomial_distribution_compare",minN,"-",maxN,".jpeg", 
#                                           sep=""), quality=90)
x11()
par(cex.axis=1.5, cex.lab=2, cex.main=2, mar=c(5,5,5,5))

data.match=data[data$total <= maxN & data$total >=minN, ]
h = hist(data.match$allelicRatio, xlim=range(0,1),breaks=binSize,right=TRUE)
# plot right closed interval left open (range] for x axis
# note that you have pseudozeroes as counts in your data so thats fine
barplot(h$counts/sum(h$counts), ylab='density', xlab='allelicRatio', 
        names.arg=h$mids, ylim=c(0,yuplimit), main=paste("n=",minN,'-',maxN))

### expected binomial 
## plot the binomial distribution for each n
## dbinom(seq(0,500),n=500,p=0.5) gives the pdf of 0-500, at n=500, p=0.5 in binomial distrib
## weight each probability by the number of SNPs with n reads

# weight by empirical counts
t = as.data.frame(table(data$total), stringsAsFactors=F)
w = matrix(0,max(data$total),1)

for (jj in 1:nrow(t))
{
  w[as.integer(t[jj,1]),1] = t[jj,2]
}


## weighted binomial distribution
# d.combined collects all results and correspond it to an allelic ratio:
# col1=allelicRatio (based on binomial n=6, ar=0,1/6,2/6...)
# col2=corresponding weighted value in binomial distribution, i.e. pdf(n,k,p)*(num of empirical SNPs at n counts)
gc()

d.combined = matrix(0,sum(seq(minN+1,maxN+1)),2)
ptr = 1

system.time({
for (i in minN:maxN)
{  
  ## doing the distribution
  k=seq(0,i)
  d = dbinom(k,i,p) ## binomial
  
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
  ## debug
#   plot(k,d,type="p",xlab='',ylab='',cex.lab=2, 
#        cex.axis=2, pch=16 , bty='n', 
#        xaxt='n', yaxt='n',col='red', ylim=c(0,yuplimit))
#   par(new=TRUE)
#   
#   plot(k,d.w,type="p",xlab='',ylab='',cex.lab=2, 
#        cex.axis=2, pch=16 , bty='n', 
#        xaxt='n', yaxt='n',col='blue', ylim=c(0,yuplimit))
#   par(new=TRUE)
#   
#   Sys.sleep(0.5)
}
})

gc()

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

par(cex.axis=1.5, cex.lab=2, cex.main=2, mar=c(5,5,5,5))
par(new=TRUE)
plot(d.combined.sorted.binned,ylim=c(0,yuplimit),pch=16,type='b',col='red',
     bty='n',ylab='',xlab='',yaxt='n',xaxt='n',yaxs="i")

# dev.off()

gc()

##################################################################################
## total
# max = 100
# data2 = data
# data2$total[data2$total>max] = max
# x11()
# par(cex.axis=1.5, cex.lab=2, cex.main=2, mar=c(5,5,5,5))
# h = hist(data2$total, breaks=40)
# barplot(h$counts/sum(h$counts), ylab='density', xlab='total',
#         names.arg=h$mids)
