setwd("C:/Users/JM/thesis/mark_work/allele_specificity/read depth")

### set parameters
title = 'NA12878_ase_polyA'
data = read.table("counts.min6.allelicRatio.mod.auto.txt", header=T, stringsAsFactors=F)
# data = read.table("test.txt", header=T, stringsAsFactors=F)

colors <- c("green","red","blue","orange","cyan","pink","purple",
            "brown","black","slategray1","violetred","tan","deeppink","darkgreen", 
            "orchid","darksalmon","antiquewhite3","magenta","darkblue","peru","slateblue",
            "thistle","tomato","rosybrown1","royalblue","olivedrab") ##Set of 27 colours to use for all plots

p=0.5
minN=10             #min total num of reads (since it's left open right closed)
maxN=12 #max total num of reads
# maxN=12 #max total num of reads
yuplimit=0.17
binSize=40

## empirical allelic Ratio
x11()
# jpeg(filename=paste("allele_readdepth_table_binomial_distribution_compare",minN,"-",maxN,".jpeg", 
#                                                           sep=""), quality=90)
par(cex.axis=1.5, cex.lab=2, cex.main=2, mar=c(5,5,5,5))
h = hist(data$allelicRatio, xlim=range(0,1),breaks=binSize,right=TRUE)
# plot right closed interval left open (range] for x axis
# note that you have pseudozeroes as counts in your data so thats fine
barplot(h$counts/sum(h$counts), ylab='density', xlab='allelicRatio', 
        names.arg=h$mids, ylim=c(0,yuplimit))

### expected binomial 
## plot the binomial distribution for each n
## dbinom(seq(0,500),n=500,p=0.5) gives the pdf of 0-500, at n=500, p=0.5 in binomial distrib
## weight each probability by the number of SNPs with n reads

# weight by empirical counts
# sort the combined distribution of all the n's
combined.sorted = read.table("counts.min6.allelicRatio.mod.auto.nonzero.bin", header=T, stringsAsFactors=F)

## bin it according to empirical distribution
bins=pretty(0:1,binSize)
start=0
end=0

combined.sorted.binned = matrix(0,length(bins)-1,2)
for (z in 2:length(bins)) ##skip 0
{
    start=bins[z-1]
    end=bins[z]
    
    row=z-1
    combined.sorted.binned[row,1] = (end-start)/2 +start ## equi of a $mid in hist
    
    combined.sorted.binned[row,2] = sum(combined.sorted[(combined.sorted[,1]<=end & 
                                                             combined.sorted[,1]>start),2])
                  ## empirical right closed, left open ?hist; right=TRUE
                  ## (range] so no double counts
                  ## but zero gets excluded!!
    if(row==1)
    {
      combined.sorted.binned[row,2] = sum(combined.sorted[(combined.sorted[,1]<=end & 
                                                             combined.sorted[,1]>=start),2])
    }
                
    
}

## change "counts" into density
combined.sorted.binned[,2] = combined.sorted.binned[,2]/sum(combined.sorted.binned[,2])

par(new=TRUE)
plot(combined.sorted.binned,ylim=c(0,yuplimit),pch=16,type='b',col='red',
     bty='n',ylab='',xlab='',yaxt='n',xaxt='n',yaxs="i")

par(new=TRUE)
abline(0,0,col='red')
# dev.off()

gc()
# ##################################################################################
# ## total
# xuplimit = 200
# max_a = max(data$total)
# ave = mean(data$total)
# data2 = data
# data2$total[data2$total>xuplimit] = xuplimit
# 
# x11()
# par(cex.axis=1.5, cex.lab=2, cex.main=2, mar=c(5,5,5,5))
# h = hist(data2$total, breaks=40)
# barplot(h$counts/sum(h$counts), ylab='density', xlab='total',
#         names.arg=h$mids)
# hp = dpois(seq(1,max_a),ave)
# par(new=TRUE)
# plot(seq(1,max_a),hp,ylim=c(0,yuplimit),xlim=c(0,xuplimit),pch=16,type='b',col='red',
#      bty='n',ylab='',xlab='',yaxt='n',xaxt='n')
# par(new=TRUE)
# abline(0,0,col='red')