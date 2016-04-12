setwd('C:/Users/Jieming/Documents/thesis/lynne_work/motifVar/manuscript/figures/figure1-protConserv')

file = 'figure1c.txt'
data = read.table(file,header=T,stringsAsFactors=F)
datat = setNames(data.frame(t(data[,-1])), data[,1])

x11()
par(cex.axis=1, cex.lab=2, cex.main=2)
barplot(as.matrix(datat), xlab="residue number", ylab ="number of variants",
        col=c("red","blue","orange","black"),beside=TRUE,ylim=c(0,200))

legend(40,200,c("TTC21B_1","TTC21B_2","TTC21B_3","1TPR"),col=c("red","blue","orange","black"),
       cex=1.5, pt.cex=1.5,
       text.col = "black", pch = 15, bg = 'white')

## zoomed
x11()
par(cex.axis=1, cex.lab=2, cex.main=2)
barplot(as.matrix(datat), xlab="residue number", ylab ="number of variants",
        col=c("red","blue","orange","black"),beside=TRUE,ylim=c(0,3))

## just the 1TPR
x11()
par(mar=c(5,5,1,1),cex.axis=1, cex.lab=2, cex.main=2)
barplot(as.matrix(datat)[4,], xlab="residue number", ylab ="number of variants",
        col=c("black"))
