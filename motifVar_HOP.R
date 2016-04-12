setwd('C:/Users/JM/thesis/lynne_work/motifVar_pNets/tpr/HOP')

## data
data = read.table('variations in HOP.txt', header=T, stringsAsFactors = F)


### plot bar plots
# rare.ns vs comm.ns
rare.comm.ns = data$numRareNS.noS + data$numCommNS.noS
x11()
par(cex.axis=1.5, cex.lab=2, cex.main=2)
barplot(rare.comm.ns,names.arg=data$resNum,main='variants by tpr motif pos',col='blue')
barplot(data$numRareNS.noS,add=T,xaxt="n",col='red') # tab.var.noSingle


legend(3,12,c("rare.ns","comm.ns"),
       col=c("red","blue"), cex=2, pt.cex=2,
       text.col = "black", pch = 15, bg = 'white')

# rare.s vs comm.s
rare.comm.s = data$numRareS.noS + data$numCommS.noS
x11()
par(cex.axis=1.5, cex.lab=2, cex.main=2)
barplot(rare.comm.s,names.arg=data$resNum,main='variants by tpr motif pos',col='cyan')
barplot(data$numRareNS.noS,add=T,xaxt="n",col='magenta') # tab.var.noSingle


legend(1,10,c("rare.s","comm.s"),
       col=c("magenta","cyan"), cex=2, pt.cex=2,
       text.col = "black", pch = 15, bg = 'white')