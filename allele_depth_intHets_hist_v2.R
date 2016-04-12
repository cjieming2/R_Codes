library(graphics)
source("C:/Users/Jieming/Documents/Shared/scripts-R_perl_shell_macros/R codes/jmRlib.R")
setwd('C:/Users/Jieming/Documents/thesis/mark_work/allele_specificity/datasets_accN/z_statistics')
all_intHets_asb = read.table('sorted.11samples.interestingHets.betabinom.peaks.chipseq.full.auto.bed',header=F,stringsAsFactors = F)
all_intHets_ase = read.table('sorted.379samples.interestingHets.betabinom.rnaseq.full.auto.bed',header=F,stringsAsFactors = F)
all_accHets_asb = read.table('sorted.11samples.acc.betabinom.peaks.chipseq.full.auto.bed',header=F,stringsAsFactors=F)
all_accHets_ase = read.table('sorted.379samples.acc.betabinom.rnaseq.full.auto.bed',header=F,stringsAsFactors=F)

asb_samples=unique(all_accHets_asb$V20)
ase_samples=unique(all_accHets_ase$V20)

## common parameters ###
numBins = 200000
xlimit.small = c(0,100)
xlimit.large = c(0,500)

#### ASE accHets ####
## calc total number of ASE acc reads between max and second highest
## table()[1:5] gives  3     4     5     6     7 
##                     5   470 61910 41994 31077 
## the min should be 8; which means that there are some sites that are tri/quad allelic
## but a small number of sites compared to the total - 135456/5714785=2%
# maxAlleleNum.acc.ase = apply(all_accHets_ase[,11:14],1,max)
# secondmaxAlleleNum.acc.ase = apply(all_accHets_ase[,11:14],1,max2)
# totalAlleleCount.acc.ase = maxAlleleNum.acc.ase + secondmaxAlleleNum.acc.ase

## use the sum of the 4 instead then since that's you used before
totalAlleleCount.acc.ase = apply(all_accHets_ase[,11:14],1,sum)

## plot hist: read depth and %significant sites
x11()
par(mar=c(5,5,1,1),cex.axis=1, cex.lab=2, cex.main=2)
hh.acc.ase = hist(totalAlleleCount.acc.ase,numBins)

## smalllimit
# density
# barplot(hh.acc.ase$counts/sum(hh.acc.ase$counts), names.arg=hh.acc.ase$mids, xlim=xlimit.small, las=2, xlab="read depth", ylab="percentage ASE sites")
# counts
ylimit=c(0,max(hh.acc.ase$counts))
barplot(hh.acc.ase$counts, names.arg=hh.acc.ase$mids, xlim=xlimit.small, ylim=ylimit,las=2, 
        col=add.alpha("cadetblue",alpha=0.5), border=add.alpha("cadetblue",alpha=0.5),
        xlab="read depth", ylab="counts ASE sites")
## largelimit
# barplot(hh.acc.ase$counts/sum(hh.acc.ase$counts), xlim=xlimit.large, las=2, xlab="read depth", ylab="percentage ASE sites")


#### ASE intHets ####
## calc total number of ASE reads between max and second highest
## table()[1:5] gives  8     9    10    11    12 
##                   186 22870 20262 16496 13364
## but to be consistent with above we use sum
# maxAlleleNum.ase = apply(all_intHets_ase[,11:14],1,max)
# secondmaxAlleleNum.ase = apply(all_intHets_ase[,11:14],1,max2)
# totalAlleleCount.ase = maxAlleleNum.ase + secondmaxAlleleNum.ase

## use the sum of the 4 instead then since that's you used before
totalAlleleCount.ase = apply(all_intHets_ase[,11:14],1,sum)

## plot hist: read depth and %significant sites
# x11()
# par(mar=c(5,5,1,1),cex.axis=1, cex.lab=2, cex.main=2)
hh.ase = hist(totalAlleleCount.ase,numBins,plot=F)

## smalllimit
# density
# barplot(hh.ase$counts/sum(hh.ase$counts), names.arg=hh.ase$mids, xlim=xlimit.small, las=2, xlab="read depth", ylab="percentage ASE sites")
# counts
barplot(hh.ase$counts, names.arg=hh.ase$mids, xlim=xlimit.small, ylim=ylimit, las=2, 
        col=add.alpha("cadetblue4",alpha=1) , border="black", xlab="",ylab="", xaxt="n", yaxt="n", 
        add=T)
## largelimit
# barplot(hh.ase$counts/sum(hh.ase$counts), xlim=xlimit.large, las=2, xlab="read depth", ylab="percentage ASE sites")

legend(50,250000,c("accHets","intHets"),col=c(add.alpha("cadetblue",alpha=0.5),"cadetblue4"),pch=15,cex=2)

# counts zoomed in
ylimit.zoomed=c(0,25000)
x11()
par(mar=c(5,5,1,1),cex.axis=1, cex.lab=2, cex.main=2)
barplot(hh.acc.ase$counts, names.arg=hh.acc.ase$mids, xlim=xlimit.small, ylim=ylimit.zoomed,las=2, 
        col=add.alpha("cadetblue",alpha=0.5), border=add.alpha("grey50",alpha=0.5),
        xlab="read depth", ylab="counts ASE sites")
barplot(hh.ase$counts, names.arg=hh.ase$mids, xlim=xlimit.small, ylim=ylimit.zoomed, las=2, 
        col=add.alpha("cadetblue4",alpha=1) , border="black", xlab="",ylab="", xaxt="n", yaxt="n", 
        add=T)