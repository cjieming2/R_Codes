setwd("C:/Users/Jieming/Documents/thesis/mark_work/allele_specificity/overdispersion_read depth")

################################## MAIN #######################################
### set parameters
filename = "zb_chosen_combined.sizeDataset.txt"
data = read.table(filename, header=T, stringsAsFactors=F)

## run parameters
colors <- c("red","green","blue","orange","cyan","pink","purple",
            "brown","black","slategray1","violetred","tan","deeppink","darkgreen", 
            "orchid","darksalmon","antiquewhite3","magenta","darkblue","peru","slateblue",
            "thistle","tomato","rosybrown1","royalblue","olivedrab") ##Set of 26 colours (no red) to use for all plots

binSize = 40

## histogram of b parameter
# pdf(paste(filename,"-check-",minN,"-",maxN,".pdf", sep=""),width=17, height=9)
h = hist(data$b.compare, xlim=range(0,1),breaks=pretty(0:1,binSize),right=TRUE)
i = hist(data$b.grad, xlim=range(0,1),breaks=pretty(0:1,binSize),right=TRUE)
empirical.h = h$counts/sum(h$counts)
empirical.i = i$counts/sum(i$counts)

x11(width=15,height=9)
par(cex.axis=1.5, cex.lab=2, cex.main=2, mar=c(5,5,5,5))
barplot(empirical.h, ylab='density', xlab='b parameter', 
        names.arg=h$mids, main=filename, col=rgb(0,0,1,1/4))
barplot(empirical.i, 
        names.arg=i$mids, col=rgb(1,0,0,1/4), add=T)
legend(40,0.5,c("b.compare","b.grad"),col=c("blue","red"),cex=2, pt.cex=2,
       text.col = "black", pch = 15, bg = 'white')

## distribution of b parameters
b.compare.mean.all = signif(mean(data$b.compare),3)
b.compare.mean.c = signif(mean(data$b.compare[data$assay == "chipseq"]),3)
b.compare.mean.r = signif(mean(data$b.compare[data$assay == "rnaseq"]),3)
b.compare.median.all = signif(median(data$b.compare),3)
b.compare.median.c = signif(median(data$b.compare[data$assay == "chipseq"]),3)
b.compare.median.r = signif(median(data$b.compare[data$assay == "rnaseq"]),3)
b.compare.max.all = signif(max(data$b.compare),3)
b.compare.max.c = signif(max(data$b.compare[data$assay == "chipseq"]),3)
b.compare.max.r = signif(max(data$b.compare[data$assay == "rnaseq"]),3)
b.compare.min.all = signif(min(data$b.compare),3)
b.compare.min.c = signif(min(data$b.compare[data$assay == "chipseq"]),3)
b.compare.min.r = signif(min(data$b.compare[data$assay == "rnaseq"]),3)

b.grad.mean.all = signif(mean(data$b.grad),3)
b.grad.mean.c = signif(mean(data$b.grad[data$assay == "chipseq"]),3)
b.grad.mean.r = signif(mean(data$b.grad[data$assay == "rnaseq"]),3)
b.grad.median.all = signif(median(data$b.grad),3)
b.grad.median.c = signif(median(data$b.grad[data$assay == "chipseq"]),3)
b.grad.median.r = signif(median(data$b.grad[data$assay == "rnaseq"]),3)
b.grad.max.all = signif(max(data$b.grad),3)
b.grad.max.c = signif(max(data$b.grad[data$assay == "chipseq"]),3)
b.grad.max.r = signif(max(data$b.grad[data$assay == "rnaseq"]),3)
b.grad.min.all = signif(min(data$b.grad),3)
b.grad.min.c = signif(min(data$b.grad[data$assay == "chipseq"]),3)
b.grad.min.r = signif(min(data$b.grad[data$assay == "rnaseq"]),3)


x11(width=15,height=9)
par(cex.axis=1.5, cex.lab=2, cex.main=2, mar=c(5,5,5,5))
boxplot(b.compare~assay,data,
        boxwex=0.25, at=1:2 - 0.2,
        col  =c('blue','red'),
        xlab='assay',ylab='b parameter',main=filename)
boxplot(b.grad~assay,data,
        boxwex=0.25, at=1:2 + 0.2,
        col  =c('lightblue','magenta'), add=TRUE)
legend(1.5,0.8,c("b.compare.chipseq","b.grad.chipseq",
                 "b.compare.rnaseq","b.grad.rnaseq"),
       col=c("blue","lightblue","red","magenta"),
       cex=2, pt.cex=2,
       text.col = "black", pch = 15, bg = 'white')

## relationship between b.choice and sse
x11(width=15,height=9)
par(cex.axis=1.5, cex.lab=2, cex.main=2, mar=c(5,5,5,5))
plot(data$b.compare[data$assay == "chipseq"],data$sse.compare[data$assay == "chipseq"],col="blue",
     xlab="b parameter", ylab="SSE", ylim=c(0,max(data$sse.compare)), xlim=c(0,max(data$b.compare)))
par(new=TRUE)
plot(data$b.compare[data$assay == "rnaseq"],data$sse.compare[data$assay == "rnaseq"],col="red",
     bty='n',ylab='',xlab='',yaxt='n',xaxt='n',
     ylim=c(0,max(data$sse.compare)), xlim=c(0,max(data$b.compare)))

legend(0.7,0.02,c("chipseq","rnaseq"),
         col=c("blue","red"), cex=2, pt.cex=2,
         text.col = "black", pch = 15, bg = 'white')

# dev.copy2pdf(file = paste(filename,"-check-",minN,"-",maxN,".pdf", sep=""))
# dev.off()
identify(data$b.compare,data$sse.compare,labels=data$sample)

####################################################################################
## correlation of size of dataset (# reads) vs b.choice
x11(width=15,height=9)
par(cex.axis=1.5, cex.lab=2, cex.main=2, mar=c(5,5,5,5))
plot(data$b.compare[data$assay == "chipseq"],data$datasetSize[data$assay == "chipseq"],col="blue",
     xlab="b parameter", ylab="datasetSize", ylim=c(0,max(data$datasetSize)), xlim=c(0,max(data$b.compare)))
par(new=TRUE)
plot(data$b.compare[data$assay == "rnaseq"],data$datasetSize[data$assay == "rnaseq"],col="red",
     bty='n',ylab='',xlab='',yaxt='n',xaxt='n',
     ylim=c(0,max(data$datasetSize)), xlim=c(0,max(data$b.compare)))

c = cor(data$b.compare,data$datasetSize, method="spearman")
c.c = cor(data$b.compare[data$assay == "chipseq"],data$datasetSize[data$assay == "chipseq"], method="spearman")
c.r = cor(data$b.compare[data$assay == "rnaseq"],data$datasetSize[data$assay == "rnaseq"], method="spearman")
text(0.6,5e09,bquote("Spearman.all," * rho * "=" * .(signif(c,3))))
text(0.6,4.7e09,bquote("Spearman.chipseq," * rho * "=" * .(signif(c.c,3))))
text(0.6,4.4e09,bquote("Spearman.rnaseq," * rho * "=" * .(signif(c.r,3))))
identify(data$b.compare,data$datasetSize,labels=data$sample)

x11(width=15,height=9)
par(cex.axis=1.5, cex.lab=2, cex.main=2, mar=c(5,5,5,5))
plot(data$b.grad[data$assay == "chipseq"],data$datasetSize[data$assay == "chipseq"],col="blue",
     xlab="b parameter", ylab="datasetSize", ylim=c(0,max(data$datasetSize)), xlim=c(0,max(data$b.grad)))
par(new=TRUE)
plot(data$b.grad[data$assay == "rnaseq"],data$datasetSize[data$assay == "rnaseq"],col="red",
     bty='n',ylab='',xlab='',yaxt='n',xaxt='n',
     ylim=c(0,max(data$datasetSize)), xlim=c(0,max(data$b.grad)))

c = cor(data$b.grad,data$datasetSize, method="spearman")
c.c = cor(data$b.grad[data$assay == "chipseq"],data$datasetSize[data$assay == "chipseq"], method="spearman")
c.r = cor(data$b.grad[data$assay == "rnaseq"],data$datasetSize[data$assay == "rnaseq"], method="spearman")
text(0.6,5e09,bquote("Spearman.all," * rho * "=" * .(signif(c,3))))
text(0.6,4.7e09,bquote("Spearman.chipseq," * rho * "=" * .(signif(c.c,3))))
text(0.6,4.4e09,bquote("Spearman.rnaseq," * rho * "=" * .(signif(c.r,3))))
