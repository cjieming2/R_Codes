setwd('C:/Users/JM/thesis/lynne_work/motifVar_pNets/tpr/')

## data
data = read.table('hist.1KG.snps.nonmono.smartDomain2gPos.TPR.34aa.bed', header=T, stringsAsFactors = F)

tab.var                = table(data$resNum)
tab.var.rare           = table(subset(data,data$maf<0.005)$resNum)
tab.var.rare.noSingle  = table(subset(data,data$maf<0.005 & data$maf>0.0005)$resNum)

tab.var.noSingle       = table(subset(data,data$maf>0.0005)$resNum) ## maf < 0.0005 ==> 1.092 singleton

## NS includes deletionframeshift, prematreStop, splieceoverlap
tab.var.NS             = table(subset(data,data$NS=="nonsynonymous" | data$NS=="deletionFS" | 
                                        data$NS=="prematureStop" | data$NS=="spliceOverlap")$resNum)
## only 1 in TPR and not a singleton
tab.var.delFS          = table(subset(data,data$NS=="deletionFS")$resNum)
## 5, and all singletons
## to accurately make the plot add 1 to noSingleton for that position manually in the file
tab.var.preStop        = table(subset(data,data$NS=="prematureStop")$resNum)
## 1, singleton
## to accurately make the plot add 1 to noSingleton for that position manually in the file
tab.var.splice         = table(subset(data,data$NS=="spliceOverlap")$resNum)
tab.var.NS.noSingle    = table(subset(data,(data$NS=="nonsynonymous" | data$NS=="deletionFS" | 
                                        data$NS=="prematureStop" | data$NS=="spliceOverlap") & data$maf>0.0005)$resNum)

## histograms
hist.var      = read.table('hist.1KG.tpr.variants.snps.txt',header=F,stringsAsFactors=F)
hist.var.rare = read.table('hist.1KG.tpr.variants.rare.txt',header=F,stringsAsFactors=F)
hist.var.ns   = read.table('hist.1KG.tpr.variants.ns.txt',header=F,stringsAsFactors=F)

## plot
# x11()
# par(cex.axis=1.5, cex.lab=2, cex.main=2)
# plot(tab.var, ylab='counts', xlab='TPR motif position', xaxt="n", col="blue", 
#      cex.axis=1.5,cex.lab=1.5, type="h")
# axis(1, at=seq(1,34), labels=seq(1,34))
# 
# par(new=T)
# 
# lines(tab.var.rare,col="red")

# legend(15, 25000, c('accHets','intHets'), 
#        col=c('blue','red'), 
#        text.col = "black", pch = 20, bg = 'white', horiz = 1, cex=1.5)
x11()
par(cex.axis=1.5, cex.lab=2, cex.main=2)
barplot(as.matrix(hist.var[2,]),names.arg=hist.var[1,],main='variants by tpr motif pos')
barplot(as.matrix(hist.var[3,]),add=T,xaxt="n",col='red') # tab.var.noSingle
barplot(as.matrix(hist.var[4,]),add=T,xaxt="n",col='blue') # tab.var.rare.noSingle

# legend(3,40,c("totvar","totvar.noSingle","totvar.noSingle.rare"),
#        col=c('darkgrey','red','blue'), cex=2, pt.cex=2,
#        text.col = "black", pch = 15, bg = 'white')
legend(3,40,c("Singletons","notSingletons","rare.notSingletons"),
              col=c('darkgrey','red','blue'), cex=2, pt.cex=2,
              text.col = "black", pch = 15, bg = 'white')

x11()
par(cex.axis=1.5, cex.lab=2, cex.main=2)
barplot(as.matrix(hist.var.rare[2,]),names.arg=hist.var.rare[1,],main='variants by tpr motif pos')
barplot(as.matrix(hist.var.rare[3,]),add=T,xaxt="n",col='green') # tab.var.rare
barplot(as.matrix(hist.var.rare[4,]),add=T,xaxt="n",col='blue')  # tab.var.rare.noSingle

# legend(3,40,c("totvar","totvar.rare","totvar.rare.noSingle"),
#        col=c('darkgrey','green','blue'), cex=2, pt.cex=2,
#        text.col = "black", pch = 15, bg = 'white')
legend(3,40,c("common","rare.singletons","rare.notSingletons"),
       col=c('darkgrey','green','blue'), cex=2, pt.cex=2,
       text.col = "black", pch = 15, bg = 'white')


x11()
par(cex.axis=1.5, cex.lab=2, cex.main=2)
barplot(as.matrix(hist.var.ns[2,]),names.arg=hist.var.ns[1,],main='variants by tpr motif pos')
barplot(as.matrix(hist.var.ns[3,]),add=T,xaxt="n",col='magenta') # ALL NS
barplot(as.matrix(hist.var.ns[4,]),add=T,xaxt="n",col='orange') # prematurestop
barplot(as.matrix(hist.var.ns[5,]),add=T,xaxt="n",col='deepskyblue') # spliceOverlap
barplot(as.matrix(hist.var.ns[6,]),add=T,xaxt="n",col='pink') # nosingle
barplot(as.matrix(hist.var.ns[7,]),add=T,xaxt="n",col='brown') # deletionFS (included as nt a singleton)

# legend(3,45,c("tot.var","totvar.NS","totvar.NS.noSingle","prematureStop.sing","spliceOverlap.sing","deletionFS.nsing"),
#        col=c('gray22','magenta','pink',"deepskyblue","orange","brown"), cex=2, pt.cex=2,
#        text.col = "black", pch = 15, bg = 'white')

legend(3,40,c("S","NS.singletons","NS.notSingletons","prematureStop.sing","spliceOverlap.sing","deletionFS.nsing"),
       col=c('gray22','magenta','pink',"deepskyblue","orange","brown"), cex=2, pt.cex=2,
       text.col = "black", pch = 15, bg = 'white')

########################
## combinations
### with singletons
########################
rare.ns = table(subset(data,data$maf<0.005 & (data$NS=="nonsynonymous" | data$NS=="deletionFS" | 
                                                data$NS=="prematureStop" | data$NS=="spliceOverlap"))$resNum)
rare.s  = table(subset(data,data$maf<0.005 & data$NS=="synonymous")$resNum)
comm.ns = table(subset(data,data$maf>=0.005 & (data$NS=="nonsynonymous" | data$NS=="deletionFS" | 
                                                 data$NS=="prematureStop" | data$NS=="spliceOverlap"))$resNum)
comm.s  = table(subset(data,data$maf>=0.005 & data$NS=="synonymous")$resNum)

hist.combo = read.table('hist.1KG.tpr.variants.combo.txt',header=F,stringsAsFactors=F)

x11()
par(cex.axis=1.5, cex.lab=2, cex.main=2)
# rare.ns, rare.s, comm.ns, comm.s
barplot(as.matrix(hist.combo[2:5,]),names.arg=hist.var.ns[1,],main='combo variants by tpr motif pos/withSingletons',
        col=c("red","pink","blue","cyan"))

# legend(3,45,c("tot.var","totvar.NS","totvar.NS.noSingle","prematureStop.sing","spliceOverlap.sing","deletionFS.nsing"),
#        col=c('gray22','magenta','pink',"deepskyblue","orange","brown"), cex=2, pt.cex=2,
#        text.col = "black", pch = 15, bg = 'white')

legend(3,40,c("rare.ns","rare.s","comm.ns","comm.s"),
       col=c("red","pink","blue","cyan"), cex=2, pt.cex=2,
       text.col = "black", pch = 15, bg = 'white')

## ns only
x11()
par(cex.axis=1.5, cex.lab=2, cex.main=2)
hist.combo.ns = rbind(hist.combo[2,],hist.combo[4,])
# rare.ns, rare.s, comm.ns, comm.s
barplot(as.matrix(hist.combo.ns),names.arg=hist.var.ns[1,],main='combo variants by tpr motif pos/withSingletons',
        col=c("red","blue"), ylim=c(0,40))

# legend(3,45,c("tot.var","totvar.NS","totvar.NS.noSingle","prematureStop.sing","spliceOverlap.sing","deletionFS.nsing"),
#        col=c('gray22','magenta','pink',"deepskyblue","orange","brown"), cex=2, pt.cex=2,
#        text.col = "black", pch = 15, bg = 'white')
legend(3,40,c("rare.ns","comm.ns"),
       col=c("red","blue"), cex=2, pt.cex=2,
       text.col = "black", pch = 15, bg = 'white')





### without Singletons
rare.ns = table(subset(data,data$maf<0.005 & data$maf>0.0005 & (data$NS=="nonsynonymous" | data$NS=="deletionFS" | 
                         data$NS=="prematureStop" | data$NS=="spliceOverlap"))$resNum)
rare.s  = table(subset(data,data$maf<0.005 & data$maf>0.0005 & data$NS=="synonymous")$resNum)
comm.ns = table(subset(data,data$maf>=0.005 & data$maf>0.0005 & (data$NS=="nonsynonymous" | data$NS=="deletionFS" | 
                         data$NS=="prematureStop" | data$NS=="spliceOverlap"))$resNum)
comm.s  = table(subset(data,data$maf>=0.005 & data$maf>0.0005 & data$NS=="synonymous")$resNum)

hist.combo = read.table('hist.1KG.tpr.variants.combo.noSingle.txt',header=F,stringsAsFactors=F)

x11()
par(cex.axis=1.5, cex.lab=2, cex.main=2)
# rare.ns, rare.s, comm.ns, comm.s
barplot(as.matrix(hist.combo[2:5,]),names.arg=hist.var.ns[1,],main='combo variants by tpr motif pos/w/oSingletons',
        col=c("red","pink","blue","cyan"), ylim=c(0,40))

# legend(3,45,c("tot.var","totvar.NS","totvar.NS.noSingle","prematureStop.sing","spliceOverlap.sing","deletionFS.nsing"),
#        col=c('gray22','magenta','pink',"deepskyblue","orange","brown"), cex=2, pt.cex=2,
#        text.col = "black", pch = 15, bg = 'white')

legend(3,30,c("rare.ns","rare.s","comm.ns","comm.s"),
       col=c("red","pink","blue","cyan"), cex=2, pt.cex=2,
       text.col = "black", pch = 15, bg = 'white')

## ns only
x11()
par(cex.axis=1.5, cex.lab=2, cex.main=2)
hist.combo.ns = rbind(hist.combo[2,],hist.combo[4,])
# rare.ns, rare.s, comm.ns, comm.s
barplot(as.matrix(hist.combo.ns),names.arg=hist.var.ns[1,],main='combo variants by tpr motif pos/w/oSingletons',
        col=c("red","blue"), ylim=c(0,40))



# legend(3,45,c("tot.var","totvar.NS","totvar.NS.noSingle","prematureStop.sing","spliceOverlap.sing","deletionFS.nsing"),
#        col=c('gray22','magenta','pink',"deepskyblue","orange","brown"), cex=2, pt.cex=2,
#        text.col = "black", pch = 15, bg = 'white')

legend(3,30,c("rare.ns","comm.ns"),
       col=c("red","blue"), cex=2, pt.cex=2,
       text.col = "black", pch = 15, bg = 'white')


## rare ns only
x11()
par(cex.axis=1.5, cex.lab=2, cex.main=2)
par(mar=c(5,5,5,0))
hist.combo.rare.ns = rbind(hist.combo[2,])
# rare.ns, rare.s, comm.ns, comm.s
barplot(as.matrix(hist.combo.rare.ns),names.arg=hist.var.ns[1,],main='combo variants by tpr motif pos/w/oSingletons rare NS',
        col=c("black"), ylim=c(0,40), xlab="TPR motif position", ylab="number of variants")



