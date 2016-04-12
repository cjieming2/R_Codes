setwd('C:/Users/JM/thesis/lynne_work/motifVar_pNets/tpr/')

## relative entropy data
sorted.rel.ent = read.table('1KG.snps.nonmono.smartDomain2gPos.TPR.34aa.sorted.sorting', header=T, stringsAsFactors = F)

######ratio.comm2rare.noS
column = sorted.rel.ent$ratio.comm2rare.noS
colname = "ratio.comm2rare.noS"
sorted.rel.ent <- sorted.rel.ent[order(column,decreasing=TRUE), ]
column = sorted.rel.ent$ratio.comm2rare.noS

## calc total entropy for each position/row
x11()
par(cex.axis=0.9, cex.lab=2, cex.main=2,mar=c(5,5,4,3),xpd=TRUE)

# plot bar for info content
ff=barplot(sorted.rel.ent$relativeEntropy,names.arg=sorted.rel.ent$resNum,main=colname,
           ylab="relative entropy (bits)", xlab="TPR motif position", 
           col=sorted.rel.ent$color, beside=T)
par(new=TRUE)

# plot line for ratio comm2rare; denominator is totvar.noSingle

plot(ff+0.5, column, type="b", pch=16,col='green',
     xaxt="n", xlab="", yaxt="n", ylab="", bty="n", ylim=c(0,3))
axis(side=4, at = seq(0,3,0.1))
mtext(colname, side=4, line=3,cex=1.5)
c = cor(sorted.rel.ent$relativeEntropy, column, method="pearson")
text(10,1.5,paste("r = ", round(c,3)),cex=1.5,col=c('green'))

######ratio.comm2rare
column = sorted.rel.ent$ratio.comm2rare
colname = "ratio.comm2rare"
sorted.rel.ent <- sorted.rel.ent[order(column,decreasing=TRUE), ]
column = sorted.rel.ent$ratio.comm2rare

## calc total entropy for each position/row
x11()
par(cex.axis=0.9, cex.lab=2, cex.main=2,mar=c(5,5,4,3),xpd=TRUE)

# plot bar for info content
ff=barplot(sorted.rel.ent$relativeEntropy,names.arg=sorted.rel.ent$resNum,main=colname,
           ylab="relative entropy (bits)", xlab="TPR motif position", 
           col=sorted.rel.ent$color, beside=T)
par(new=TRUE)

# plot line for ratio comm2rare; denominator is totvar.noSingle

plot(ff+0.5, column, type="b", pch=16,col='green',
     xaxt="n", xlab="", yaxt="n", ylab="", bty="n", ylim=c(0,1))
axis(side=4, at = seq(0,1,0.1))
mtext(colname, side=4, line=3,cex=1.5)
c = cor(sorted.rel.ent$relativeEntropy, column, method="pearson")
text(10,0.9,paste("r = ", round(c,3)),cex=1.5,col=c('green'))

######P.NS.noS
column = sorted.rel.ent$P.NS.noS
colname = "P.NS.noS"
sorted.rel.ent <- sorted.rel.ent[order(column,decreasing=TRUE), ]
column = sorted.rel.ent$P.NS.noS

## calc total entropy for each position/row
x11()
par(cex.axis=0.9, cex.lab=2, cex.main=2,mar=c(5,5,4,3),xpd=TRUE)

# plot bar for info content
ff=barplot(sorted.rel.ent$relativeEntropy,names.arg=sorted.rel.ent$resNum,main=colname,
           ylab="relative entropy (bits)", xlab="TPR motif position", 
           col=sorted.rel.ent$color, beside=T)
par(new=TRUE)

# plot line for ratio comm2rare; denominator is totvar.noSingle

plot(ff+0.5, column, type="b", pch=16,col='green',
     xaxt="n", xlab="", yaxt="n", ylab="", bty="n", ylim=c(0,1))
axis(side=4, at = seq(0,1,0.1))
mtext(colname, side=4, line=3,cex=1.5)
c = cor(sorted.rel.ent$relativeEntropy, column, method="pearson")
text(10,0.9,paste("r = ", round(c,3)),cex=1.5,col=c('green'))

######P.NS
column = sorted.rel.ent$P.NS
colname = "P.NS"
sorted.rel.ent <- sorted.rel.ent[order(column,decreasing=TRUE), ]
column = sorted.rel.ent$P.NS

## calc total entropy for each position/row
x11()
par(cex.axis=0.9, cex.lab=2, cex.main=2,mar=c(5,5,4,3),xpd=TRUE)

# plot bar for info content
ff=barplot(sorted.rel.ent$relativeEntropy,names.arg=sorted.rel.ent$resNum,main=colname,
           ylab="relative entropy (bits)", xlab="TPR motif position", 
           col=sorted.rel.ent$color, beside=T)
par(new=TRUE)

# plot line for ratio comm2rare; denominator is totvar.noSingle

plot(ff+0.5, column, type="b", pch=16,col='green',
     xaxt="n", xlab="", yaxt="n", ylab="", bty="n", ylim=c(0,1))
axis(side=4, at = seq(0,1,0.1))
mtext(colname, side=4, line=3,cex=1.5)
c = cor(sorted.rel.ent$relativeEntropy, column, method="pearson")
text(10,0.9,paste("r = ", round(c,3)),cex=1.5,col=c('green'))

######P.rare.noS
column = sorted.rel.ent$P.rare.noS
colname = "P.rare.noS"
sorted.rel.ent <- sorted.rel.ent[order(column,decreasing=TRUE), ]
column = sorted.rel.ent$P.rare.noS

## calc total entropy for each position/row
x11()
par(cex.axis=0.9, cex.lab=2, cex.main=2,mar=c(5,5,4,3),xpd=TRUE)

# plot bar for info content
ff=barplot(sorted.rel.ent$relativeEntropy,names.arg=sorted.rel.ent$resNum,main=colname,
           ylab="relative entropy (bits)", xlab="TPR motif position", 
           col=sorted.rel.ent$color, beside=T)
par(new=TRUE)

# plot line for ratio comm2rare; denominator is totvar.noSingle

plot(ff+0.5, column, type="b", pch=16,col='green',
     xaxt="n", xlab="", yaxt="n", ylab="", bty="n", ylim=c(0,1))
axis(side=4, at = seq(0,1,0.1))
mtext(colname, side=4, line=3,cex=1.5)
c = cor(sorted.rel.ent$relativeEntropy, column, method="pearson")
text(10,20,paste("r = ", round(c,3)),cex=1.5,col=c('green'))

######P.rare
column = sorted.rel.ent$P.rare
colname = "P.rare"
sorted.rel.ent <- sorted.rel.ent[order(column,decreasing=TRUE), ]
column = sorted.rel.ent$P.rare

## calc total entropy for each position/row
x11()
par(cex.axis=0.9, cex.lab=2, cex.main=2,mar=c(5,5,4,3),xpd=TRUE)

# plot bar for info content
ff=barplot(sorted.rel.ent$relativeEntropy,names.arg=sorted.rel.ent$resNum,main=colname,
           ylab="relative entropy (bits)", xlab="TPR motif position", 
           col=sorted.rel.ent$color, beside=T)
par(new=TRUE)

# plot line for ratio comm2rare; denominator is totvar.noSingle

plot(ff+0.5, column, type="b", pch=16,col='green',
     xaxt="n", xlab="", yaxt="n", ylab="", bty="n", ylim=c(0,1))
axis(side=4, at = seq(0,1,0.1))
mtext(colname, side=4, line=3,cex=1.5)
c = cor(sorted.rel.ent$relativeEntropy, column, method="pearson")
text(10,9,paste("r = ", round(c,3)),cex=1.5,col=c('green'))

######P.rare.noS
column = sorted.rel.ent$P.rare.noS
colname = "P.rare.noS"
sorted.rel.ent <- sorted.rel.ent[order(column,decreasing=TRUE), ]
column = sorted.rel.ent$P.rare.noS

## calc total entropy for each position/row
x11()
par(cex.axis=0.9, cex.lab=2, cex.main=2,mar=c(5,5,4,3),xpd=TRUE)

# plot bar for info content
ff=barplot(sorted.rel.ent$relativeEntropy,names.arg=sorted.rel.ent$resNum,main=colname,
           ylab="relative entropy (bits)", xlab="TPR motif position", 
           col=sorted.rel.ent$color, beside=T)
par(new=TRUE)

# plot line for ratio comm2rare; denominator is totvar.noSingle

plot(ff+0.5, column, type="b", pch=16,col='green',
     xaxt="n", xlab="", yaxt="n", ylab="", bty="n", ylim=c(0,1))
axis(side=4, at = seq(0,1,0.1))
mtext(colname, side=4, line=3,cex=1.5)
c = cor(sorted.rel.ent$relativeEntropy, column, method="pearson")
text(10,9,paste("r = ", round(c,3)),cex=1.5,col=c('green'))


######P.RareNS.noS
column = sorted.rel.ent$P.RareNS.noS
colname = "P.RareNS.noS"
sorted.rel.ent <- sorted.rel.ent[order(column,decreasing=TRUE), ]
column = sorted.rel.ent$P.RareNS.noS

## calc total entropy for each position/row
x11()
par(cex.axis=0.9, cex.lab=2, cex.main=2,mar=c(5,5,4,3),xpd=TRUE)

# plot bar for info content
ff=barplot(sorted.rel.ent$relativeEntropy,names.arg=sorted.rel.ent$resNum,main=colname,
           ylab="relative entropy (bits)", xlab="TPR motif position", 
           col=sorted.rel.ent$color, beside=T)
par(new=TRUE)

# plot line for ratio comm2rare; denominator is totvar.noSingle

plot(ff+0.5, column, type="b", pch=16,col='green',
     xaxt="n", xlab="", yaxt="n", ylab="", bty="n", ylim=c(0,1))
axis(side=4, at = seq(0,1,0.1))
mtext(colname, side=4, line=3,cex=1.5)
c = cor(sorted.rel.ent$relativeEntropy, column, method="pearson")
text(10,0.9,paste("r = ", round(c,3)),cex=1.5,col=c('green'))

######P.RareNS
column = sorted.rel.ent$P.RareNS
colname = "P.RareNS"
sorted.rel.ent <- sorted.rel.ent[order(column,decreasing=TRUE), ]
column = sorted.rel.ent$P.RareNS

## calc total entropy for each position/row
x11()
par(cex.axis=0.9, cex.lab=2, cex.main=2,mar=c(5,5,4,3),xpd=TRUE)

# plot bar for info content
ff=barplot(sorted.rel.ent$relativeEntropy,names.arg=sorted.rel.ent$resNum,main=colname,
           ylab="relative entropy (bits)", xlab="TPR motif position", 
           col=sorted.rel.ent$color, beside=T)
par(new=TRUE)

# plot line for ratio comm2rare; denominator is totvar.noSingle

plot(ff+0.5, column, type="b", pch=16,col='green',
     xaxt="n", xlab="", yaxt="n", ylab="", bty="n", ylim=c(0,1))
axis(side=4, at = seq(0,1,0.1))
mtext(colname, side=4, line=3,cex=1.5)
c = cor(sorted.rel.ent$relativeEntropy, column, method="pearson")
text(10,0.9,paste("r = ", round(c,3)),cex=1.5,col=c('green'))