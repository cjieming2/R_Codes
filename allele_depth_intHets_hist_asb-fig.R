library(graphics)
source("C:/Users/Jieming/Documents/Shared/scripts-R_perl_shell_macros/R codes/jmRlib.R")
# setwd('C:/Users/Jieming/Documents/thesis/mark_work/allele_specificity/datasets_accN/z_statistics')
setwd('C:/Users/Jieming/Documents/thesis/mark_work/allele_specificity/datasets_accN/datasets_pooled')
# all_intHets_asb = read.table('sorted.11samples.interestingHets.betabinom.peaks.chipseq.full.auto.bed',header=F,stringsAsFactors = F)
# all_conHets_asb = read.table('sorted.11samples.acc.betabinom.peaks.chipseq.full.auto.bed',header=F,stringsAsFactors=F)
all_intHets_asb = read.table('zintHets.asb.acc.13samples.merged.redun.auto.bed',header=F,stringsAsFactors = F)
all_conHets_asb = read.table('znonintHets.asb.acc.13samples.merged.redun.auto.bed',header=F,stringsAsFactors=F)

## since acc here means control SNVs you add them up to form acc
## note there are casbs where control in one sample is found as intHets in another and vice versa
## so this is a redundant set in terms of SNVs but non-redun in terms of read depth even for same SNV
## the only reason we are simply concat them, is we are assuming we are concat for each indiv samples

## intHets has one more column?!
# all_intHets_asb.mod = all_intHets_asb[,-20,drop=FALSE]
# names(all_intHets_asb.mod)[20]="V20";names(all_intHets_asb.mod)[21]="V21"


all_accHets_asb = rbind(all_intHets_asb,all_conHets_asb)

asb_samples=unique(all_accHets_asb$V21)

################################################################
###### the following section plots a histogram
###### how many intHets and accHets for each read depth
## common parameters ###
numBins = 2000
xlimit.small = c(0,100)
xlimit.large = c(0,500)

#### asb accHets ####
## calc total number of asb acc reads between max and second highest
## table()[1:5] gives  3     4     5     6     7 
##                     5   470 61910 41994 31077 
## the min should be 8; which means that there are some sites that are tri/quad allelic
## but a small number of sites compared to the total - 135456/5714785=2%
# maxAlleleNum.acc.asb = apply(all_accHets_asb[,11:14],1,max)
# secondmaxAlleleNum.acc.asb = apply(all_accHets_asb[,11:14],1,max2)
# totalAlleleCount.acc.asb = maxAlleleNum.acc.asb + secondmaxAlleleNum.acc.asb

## use the sum of the 4 instead then since that's you used before
totalAlleleCount.acc.asb = apply(all_accHets_asb[,11:14],1,sum)

## plot hist: read depth and %significant sites
x11()
par(mar=c(5,5,1,1),cex.axis=1, cex.lab=2, cex.main=2)
hh.acc.asb = hist(totalAlleleCount.acc.asb,numBins)

## smalllimit
# density
# barplot(hh.acc.asb$counts/sum(hh.acc.asb$counts), names.arg=hh.acc.asb$mids, xlim=xlimit.small, las=2, xlab="read depth", ylab="percentage asb sites")
# counts
ylimit=c(0,max(hh.acc.asb$counts))
mybars = barplot(hh.acc.asb$counts[1:200], names.arg=hh.acc.asb$mids[1:200], xlim=xlimit.small, ylim=ylimit,las=2, 
        col=add.alpha("aquamarine3",alpha=0.5), border=add.alpha("aquamarine3",alpha=0.5),
        xlab="read depth", ylab="number of SNVs")
## largelimit
# barplot(hh.acc.asb$counts/sum(hh.acc.asb$counts), xlim=xlimit.large, las=2, xlab="read depth", ylab="percentage asb sites")


#### asb intHets ####
## calc total number of asb reads between max and second highest
## table()[1:5] gives  8     9    10    11    12 
##                   186 22870 20262 16496 13364
## but to be consistent with above we use sum
# maxAlleleNum.asb = apply(all_intHets_asb[,11:14],1,max)
# secondmaxAlleleNum.asb = apply(all_intHets_asb[,11:14],1,max2)
# totalAlleleCount.asb = maxAlleleNum.asb + secondmaxAlleleNum.asb

## use the sum of the 4 instead then since that's you used before
totalAlleleCount.asb = apply(all_intHets_asb[,11:14],1,sum)

## plot hist: read depth and %significant sites
# x11()
# par(mar=c(5,5,1,1),cex.axis=1, cex.lab=2, cex.main=2)
hh.asb = hist(totalAlleleCount.asb,numBins/3,plot=F)

## smalllimit
# density
# barplot(hh.asb$counts/sum(hh.asb$counts), names.arg=hh.asb$mids, xlim=xlimit.small, las=2, xlab="read depth", ylab="percentage asb sites")
# counts
barplot(hh.asb$counts, names.arg=hh.asb$mids, xlim=xlimit.small, ylim=ylimit, las=2, 
        col=add.alpha("aquamarine4",alpha=1) , border="black", xlab="",ylab="", xaxt="n", yaxt="n", 
        add=T)
## largelimit
# barplot(hh.asb$counts/sum(hh.asb$counts), xlim=xlimit.large, las=2, xlab="read depth", ylab="percentage asb sites")

### write things to files debug
### comparison shows that there are no zeroes within 677
### cut fat 677
# write.table(cbind(hh.asb$mids,hh.asb$counts),file="hh.asb.txt", sep="\t",
#             row.names=FALSE,quote=FALSE)
# write.table(cbind(hh.acc.asb$mids,hh.acc.asb$counts),file="hh.acc.asb.txt", sep="\t",
#             row.names=FALSE,quote=FALSE)

## % intHets/acc
par(new=TRUE)
plot(mybars,hh.asb$counts[1:200]/hh.acc.asb$counts[1:200], type="b", col="black", 
     xlab="",ylab="",bty="n", xaxt="n", yaxt="n",
     xlim=c(0,100), ylim=c(0,0.5), yaxs="i")
axis(side=4, pos=90)
mtext("% ASBsites/AccSites",side=4,cex=2,line=-5)

legend(50,0.4,c("ACC","ASB","%ASB"),col=c(add.alpha("aquamarine3",alpha=0.5),"aquamarine4","black"),
       pch=c(15,15,1),cex=2,pt.lwd=2)


#############################################
# counts zoomed in
ylimit.zoomed=c(0,1000)
x11()
par(mar=c(5,5,1,1),cex.axis=1, cex.lab=2, cex.main=2)
barplot(hh.acc.asb$counts, names.arg=hh.acc.asb$mids, xlim=xlimit.small, ylim=ylimit.zoomed,las=2, 
        col=add.alpha("aquamarine3",alpha=0.5), border=add.alpha("aquamarine3",alpha=0.5),
        xlab="read depth", ylab="counts ASB sites")
barplot(hh.asb$counts, names.arg=hh.asb$mids, xlim=xlimit.small, ylim=ylimit.zoomed, las=2, 
        col=add.alpha("aquamarine4",alpha=1) , border="black", xlab="",ylab="", xaxt="n", yaxt="n", 
        add=T)
axis(side=4, pos=90)
mtext("% asbsites/accSites",side=4,cex=2,line=-5)


########################################################################
## plot read depth for allelic ratio ####
library(gplots)
library(spatstat)

## get only the first few bars 8:160
## massage data
roundtoDecPlaces = 1
all_intHets_asb_total = cbind(all_intHets_asb,totalAlleleCount.asb)
mat_data = all_intHets_asb_total[(all_intHets_asb_total[,23]>=8) & (all_intHets_asb_total[,23]<=160),]
## col23 - read depth ; col22 - allelic ratio
mat_table_unready = as.matrix(as.data.frame.matrix(table(mat_data[,23], round(mat_data[,22],roundtoDecPlaces))))

### when roundtoDecPlaces=1 ####
mm11 = matrix(0,nrow=length(rownames(mat_table_unready)),ncol=5) ## only 0.3 0.4 0.5 0.6 0.7 missing
mat_table = cbind(mat_table_unready[,1:3], mm11, mat_table_unready[,4:ncol(mat_table_unready)])
colnames(mat_table) = seq(0,1,0.1)
x11()
par(mar=c(5,5,1,1),cex.axis=1, cex.lab=2, cex.main=2)
cols = c(colorRampPalette(c("white","cornflowerblue"))(8), colorRampPalette(c("yellow","red"))(19))
### note that the filled contour plot.axes do not work nor using axes help
## so my solution is to not customize the axes then use illustrator to manusally change!!
filled.contour(mat_table, col=cols,
               xlab="read depth", ylab="allelic ratio",
               key.title = title(main="num asb",cex.main=1))

### when roundtoDecPlaces=2 ####
# mm1 = matrix(0,nrow=length(rownames(mat_table_unready)),ncol=39) ## missing data 0.28 to 0.66
# mm2 = matrix(0,nrow=length(rownames(mat_table_unready)),ncol=1) ## 0.68,0.70
# mat_table = cbind(mat_table_unready[,1:28], mm1, mat_table_unready[,29], mm2, mat_table_unready[,30], mm2, mat_table_unready[,31:60])
# colnames(mat_table) = seq(0,1,0.01)
# x11()
# par(mar=c(5,5,1,1),cex.axis=1, cex.lab=2, cex.main=2)
# cols = c(colorRampPalette(c("white","cornflowerblue"))(8), colorRampPalette(c("yellow","red"))(19))
# filled.contour(mat_table, col=cols,
#                plot.axes={ axis(side=1, at=seq(8,80), labels=rownames(mat_table))
#                            axis(side=2, at=seq(0,1,0.01), labels=colnames(mat_table), las=1) })

## normal plot
# x11()
# par(mar=c(5,5,1,1),cex.axis=1, cex.lab=2, cex.main=2)
# plot(mat_data[,23], mat_data[,22],xlab="read depth", ylab="allelic ratio")
# 
# ## zoomed in
# x11()
# xlimit.rd.ar = c(8,80)
# par(mar=c(5,5,1,1),cex.axis=1, cex.lab=2, cex.main=2)
# plot(mat_data[,23], mat_data[,22],
#      xlab="read depth", ylab="allelic ratio", xlim=xlimit.rd.ar)


