setwd('C:/Users/JM/thesis/lynne_work/motifVar_pNets/ank/')

## relative entropy data
sorted.rel.ent = read.table('1KG.snps.nonmono.smartDomain2gPos.ANK.30aa.mod.enrich', header=T, stringsAsFactors = F)

## calc total entropy for each position/row
x11()
par(cex.axis=1, cex.lab=2, cex.main=2,mar=c(5,5,4,1),xpd=TRUE)

# plot bar for info content
ff=barplot(sorted.rel.ent$relativeEntropy,names.arg=sorted.rel.ent$resNum,main='tpr relative entropy',
            ylab="relative entropy (bits)", xlab="TPR motif position", 
            col=c("red","red","red","red","red","red","red","red","red","red","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue"), 
            beside=T)
par(new=TRUE)

legend(17,1.5,c('conserved','none','hypervar'), 
                 col=c('red','grey','blue'), 
                 text.col = "black", pch=16, bg = 'white',cex=1.5)

# plot line for ratio comm2rare; denominator is totvar.noSingle
column = sorted.rel.ent$ratio.comm2rare
plot(ff+0.5, column, type="b", pch=16,col='red',
     xaxt="n", xlab="", yaxt="n", ylab="", bty="n", ylim=c(0,1))
axis(side=4, at = seq(0,1,0.1))
mtext("%SNVs", side=4, line=3,cex=1.5)
c = cor(sorted.rel.ent$relativeEntropy, column, method="pearson")
text(10,0.65,paste("r = ", round(c,3)),cex=1.5,col=c('red'))

# plot line for ratio comm2rare.noS; denominator is totvar.noSingle
par(new=TRUE)
column = sorted.rel.ent$ratio.comm2rare.noS
plot(ff+0.5, column, type="b", pch=16,col='blue',
     xaxt="n", xlab="", yaxt="n", ylab="", bty="n", ylim=c(0,1))
axis(side=4, at = seq(0,1,0.1))
mtext("%SNVs", side=4, line=3,cex=1.5)
c = cor(sorted.rel.ent$relativeEntropy, column, method="pearson")
text(10,0.68,paste("r = ", round(c,3)),cex=1.5,col=c('blue'))

# plot line for rare.ns.noSingle; denominator is totvar.noSingle
par(new=TRUE)
column = sorted.rel.ent$P.RareNS.noS
plot(ff+0.5, column, type="b", pch=16,col='red',
     xaxt="n", xlab="", yaxt="n", ylab="", bty="n", ylim=c(0,1))
axis(side=4, at = seq(0,1,0.1))
mtext("%SNVs", side=4, line=3,cex=1.5)
c = cor(sorted.rel.ent$relativeEntropy, column, method="pearson")
text(10,0.65,paste("r = ", round(c,3)),cex=1.5,col=c('red'))

# plot line for ns
par(new=TRUE)
column = sorted.rel.ent$numCommNS.noS / sorted.rel.ent$totVar.noS
plot(ff+0.5, column, type="b", pch=16,col='blue',
     xaxt="n", xlab="", yaxt="n", ylab="", bty="n", ylim=c(0,1))
axis(side=4, at = seq(0,1,0.1))
mtext("%SNVs", side=4, line=3,cex=1.5)
c = cor(sorted.rel.ent$relativeEntropy, column, method="pearson")
text(10,0.68,paste("r = ", round(c,3)),cex=1.5,col=c('blue'))

# plot line for ns
par(new=TRUE)
column = sorted.rel.ent$P.NS.noS
plot(ff+0.5, column, type="b", pch=16,col='darkgreen',
     xaxt="n", xlab="", yaxt="n", ylab="", bty="n", ylim=c(0,1))
axis(side=4, at = seq(0,1,0.1))
mtext("%SNVs", side=4, line=3,cex=1.5)
c = cor(sorted.rel.ent$relativeEntropy, column, method="pearson")
text(10,0.71,paste("r = ", round(c,3)),cex=1.5,col=c('darkgreen'))


# plot line for rare
par(new=TRUE)
column = sorted.rel.ent$P.rare.noS
plot(ff+0.5, column, type="b", pch=16,col='green',
     xaxt="n", xlab="", yaxt="n", ylab="", bty="n", ylim=c(0,1))
axis(side=4, at = seq(0,1,0.1))
mtext("%SNVs", side=4, line=3,cex=1.5)
c = cor(sorted.rel.ent$relativeEntropy, column, method="pearson")
text(10,0.74,paste("r = ", round(c,3)),cex=1.5,col=c('green'))

# plot line for rare.S
par(new=TRUE)
column = sorted.rel.ent$numRareS.noS / sorted.rel.ent$totVar.noS
plot(ff+0.5, column, type="b", pch=16,col='magenta',
     xaxt="n", xlab="", yaxt="n", ylab="", bty="n", ylim=c(0,1))
axis(side=4, at = seq(0,1,0.1))
mtext("%SNVs", side=4, line=3,cex=1.5)
c = cor(sorted.rel.ent$relativeEntropy, column, method="pearson")
text(10,0.77,paste("r = ", round(c,3)),cex=1.5,col=c('magenta'))

# plot line for comm.S
par(new=TRUE)
column = sorted.rel.ent$numCommS.noS / sorted.rel.ent$totVar.noS
plot(ff+0.5, column, type="b", pch=16,col='cyan',
     xaxt="n", xlab="", yaxt="n", ylab="", bty="n", ylim=c(0,1))
axis(side=4, at = seq(0,1,0.1))
mtext("%SNVs", side=4, line=3,cex=1.5)
c = cor(sorted.rel.ent$relativeEntropy, column, method="pearson")
text(10,0.80,paste("r = ", round(c,3)),cex=1.5,col=c('cyan'))

## legend
legend(8,1.07,c('information content','%rare.ns.noS','%comm.ns.noS','%ns.noS','%rare.noS', 
                                      '%rare.s.noS','%comm.s.noS'), 
       col=c('gray','red','blue','darkgreen','green','magenta','cyan'), 
       text.col = "black", pch=16, bg = 'white',cex=1.5)

# legend(8,1.05,c('information content','%rare.ns.noSingle',), 
#        col=c('gray','red','blue'), 
#        text.col = "black", pch=16, bg = 'white',cex=1.5)

### pca
