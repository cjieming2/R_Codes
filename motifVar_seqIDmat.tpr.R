source("C:/Users/Jieming/Documents/Shared/scripts-R_perl_shell_macros/R codes/jmRlib.R")

setwd('C:/Users/Jieming/Documents/thesis/lynne_work/motifVar/tpr/seqID/all')
domain = 'tprs.34aa'
snvs = read.table('TPRs.canonical.auto.122.txt', header=T, stringsAsFactors = F)
aa=hist(snvs[,3], breaks=100,plot=FALSE)
x11()
par(mar=c(5,6,4,4),xpd=TRUE,cex.axis=1, cex.lab=2, cex.main=2)
barplot(aa$counts/sum(aa$counts),names.arg=aa$mids,main=paste('distrib num SNVs ', domain, ' motif domains'),
        xlab='num SNVs', ylab='density',col='blue', las=2)

##################################################################################################################
setwd('C:/Users/Jieming/Documents/thesis/lynne_work/motifVar/tpr/blosum_sift_polyphen2_gerp/blosum_seqID')

#######################################################################
## process a nxn symmetric matrix of sequence identities
## and returns the lower triangular as a vector n*(n+1) / 2 rows * 1
processMatrix <- function(data.mod)
{
  ## remove diagonals by first setting them as NA then remove NA
  ## also transforms them into a vector
  diag(data.mod) = NA
  data.ltri=vech(data.mod)
  d <- data.ltri[!is.na(data.ltri)]
  
  return(d)
}

library(matrixcalc)
library(plotrix)
library(vioplot)
domain = 'tprs.34aa'

## data
filename = 'final.merged.ExAC.ens73.TPR.34aa.vepoutput.coding.canonical.sorted.auto.574.reduseq.seqidmat'
data = read.table(filename, sep="\t", header=T, stringsAsFactors = F)
data.mod = as.matrix(data[,2:ncol(data)])

d = processMatrix(data.mod)

## find max, min and histogram
norm.interval = function(data, variance = var(data), conf.level = 0.95) {
  z = qnorm((1 - conf.level)/2, lower.tail = FALSE)
  xbar = mean(data)
  sdx = sqrt(variance/length(t(data)))
  c(xbar - z * sdx, xbar + z * sdx)
}

## summary stats
summary.stats(d,filename)

## doing a histogram of sequences identities
a=hist(d, breaks=40,plot=FALSE)
x11()
par(mar=c(5,6,4,4),xpd=TRUE,cex.axis=1, cex.lab=2, cex.main=2)
barplot(a$counts/sum(a$counts),names.arg=a$mids,main=paste('distrib % sequence identity ', domain, ' motif pos'),
        xlab='% sequence identity', ylab='probability',col='blue', las=2)
quantiles=as.data.frame(quantile(d, probs=c(10,25,50,75)/100, type=8, names = T))
df<-data.frame(as.character(rownames(quantiles)),quantiles)
text(39,0.15,paste("max=",max(d),";min=",min(d),"\nmean=",round(mean(d),2),sep=""),cex=2)
names(df) = c("quantiles","seqIdentity")
addtable2plot(31,0.06,df, cex=2)


#########################################################################################3
## compare seq ids between each protein
setwd("C:/Users/Jieming/Documents/thesis/lynne_work/motifVar/tpr/seqID")
## process a nxn symmetric matrix of sequence identities
## and returns the lower triangular as a vector n*(n+1) / 2 rows * 1
readAndprocessMatrix <- function(filename)
{
  data = read.table(filename, sep="\t", header=T, stringsAsFactors = F)
  data.mod = as.matrix(data[,2:ncol(data)])
  
  ## remove diagonals by first setting them as NA then remove NA
  ## also transforms them into a vector
  diag(data.mod) = NA
  data.ltri=vech(data.mod)
  d <- data.ltri[!is.na(data.ltri)]
  
  return(d)
}

readAndnrow <- function(filename)
{
  data = read.table(filename, sep="\t", header=T, stringsAsFactors = F)
  return(nrow(data))
}

readAndmean <- function(filename)
{
  data = readAndprocessMatrix(filename)
  return(mean(data))
}

fnames = list.files(".", pattern="*.seqIDmat", full.names=TRUE)
frows = unlist(lapply(fnames,readAndnrow))
flist = lapply(fnames,readAndprocessMatrix)
flist.unlist = unlist(lapply(fnames,readAndprocessMatrix)) 

## plot boxplots of seq id for all ENSPs
x11()
par(cex.axis=1,cex.lab=1)
do.call(boxplot,c(flist,col="yellow",xaxt="n"))
axis(1, at=seq(1:length(fnames)),labels=frows)
mtext("TPR proteins",side=1, line=2, cex=1.5)
mtext("seq id",side=2, line=2, cex=1.5)

## plot distribution of seq id of pairwise comparisons of TPR within each ENSP for all ENSPs (seq ids within each ENSP)
hf=hist(flist.unlist, breaks=40,plot=FALSE)
x11()
par(mar=c(5,6,4,4),xpd=TRUE,cex.axis=1, cex.lab=2, cex.main=2)

## put up empty bars
newhfcounts = c(hf$counts,rep(0,length(a$mids)-length(hf$mids)))
newhfmids = c(hf$mids,seq(max(hf$mids)+0.02,1,0.02))
barplot(newhfcounts/sum(hf$counts),names.arg=newhfmids,main=paste('distrib mean % sequence identity within TPR proteins'),
        xlab='(mean) % sequence identity', ylab='probability',col=rgb(1,0,0,0.3), las=2, ylim=c(0,0.2))

## compare that with the distribution of seq id of ALL TPRs in all ENSPs
## do a comparison for seq id between and within ENSPs
barplot(a$counts/sum(a$counts),names.arg=a$mids,col=rgb(0,0,1,0.3), las=2,  ylim=c(0,0.2), xlim=c(0,1),add=T)

legend(33,0.2,c("within tpr proteins","between all tpr proteins"),
       col=c(rgb(1,0,0,0.3),rgb(0,0,1,0.3)), cex=2, pt.cex=2,
       text.col = "black", pch = 15, bg = 'white')
