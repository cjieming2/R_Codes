setwd("C:/Users/Jieming/Documents/thesis/mark_work/allele_specificity/overdispersion_read depth")

################################## MAIN #######################################
### set parameters
filename = "zb_chosen_combined.indivDataset.txt"
data = read.table(filename, header=T, stringsAsFactors=F)

## run parameters
colors <- c("red","green","blue","orange","cyan","pink","purple",
            "brown","black","slategray1","violetred","tan","deeppink","darkgreen", 
            "orchid","darksalmon","antiquewhite3","magenta","darkblue","peru","slateblue",
            "thistle","tomato","rosybrown1","royalblue","olivedrab") ##Set of 26 colours (no red) to use for all plots

binSize = 40

chipseq = as.data.frame(data[grep("chipseq",data$sample),])
rnaseq = as.data.frame(data[grep("rnaseq",data$sample),])

## histogram of b parameter
# pdf(paste(filename,"-check-",minN,"-",maxN,".pdf", sep=""),width=17, height=9)
h = hist(data$b.choice, xlim=range(0,1),breaks=pretty(0:1,binSize),right=TRUE)
empirical.h = h$counts/sum(h$counts)
i = hist(rnaseq$b.choice, xlim=range(0,1),breaks=pretty(0:1,binSize),right=TRUE)
empirical.i = i$counts/sum(i$counts)
j = hist(chipseq$b.choice, xlim=range(0,1),breaks=pretty(0:1,binSize),right=TRUE)
empirical.j = j$counts/sum(j$counts)

x11(width=15,height=9)
par(cex.axis=1.5, cex.lab=2, cex.main=2, mar=c(5,5,5,5),mfrow=c(1,3))
boxplot(data$b.choice, col  ='grey',xlab='all',ylab='b parameter',ylim=c(0,1))
boxplot(rnaseq$b.choice, col  ='red',xlab='rnaseq',ylab='',ylim=c(0,1))
boxplot(chipseq$b.choice, col  ='blue',xlab='chipseq',ylab='',ylim=c(0,1))

##barplot not useful
# barplot(empirical.h, ylab='density', xlab='b parameter', 
#         names.arg=h$mids, main=filename, col=rgb(0,0,1,1/4))
# barplot(empirical.i, 
#         names.arg=i$mids, col=rgb(0,1,0,1/4), add=T)
# barplot(empirical.j, 
#         names.arg=j$mids, col=rgb(1,0,0,1/4), add=T)

# legend(40,0.5,c("b.choice","b.choice"),col=c("blue","red"),cex=2, pt.cex=2,
#        text.col = "black", pch = 15, bg = 'white')

## distribution of b parameters
b.choice.mean.all = signif(mean(data$b.choice),3)
b.choice.mean.c = signif(mean(data$b.choice[grep("chipseq",data$sample)]),3)
b.choice.mean.r = signif(mean(data$b.choice[grep("rnaseq",data$sample)]),3)
b.choice.median.all = signif(median(data$b.choice),3)
b.choice.median.c = signif(median(data$b.choice[grep("chipseq",data$sample)]),3)
b.choice.median.r = signif(median(data$b.choice[grep("rnaseq",data$sample)]),3)
b.choice.max.all = signif(max(data$b.choice),3)
b.choice.max.c = signif(max(data$b.choice[grep("chipseq",data$sample)]),3)
b.choice.max.r = signif(max(data$b.choice[grep("rnaseq",data$sample)]),3)
b.choice.min.all = signif(min(data$b.choice),3)
b.choice.min.c = signif(min(data$b.choice[grep("chipseq",data$sample)]),3)
b.choice.min.r = signif(min(data$b.choice[grep("rnaseq",data$sample)]),3)

# ## relationship between b.choice and sse
x11(width=15,height=9)
par(cex.axis=1.5, cex.lab=2, cex.main=2, mar=c(5,5,5,5))
plot(chipseq$b.choice,chipseq$sse,col="blue",
     xlab="b parameter", ylab="SSE", ylim=c(0,max(chipseq$sse)), xlim=c(0,max(chipseq$b.choice)))
identify(chipseq$b.choice,chipseq$sse,labels=chipseq$sample)
par(new=TRUE)
plot(rnaseq$b.choice,rnaseq$sse,col="red",
     xlab="b parameter", ylab="SSE", ylim=c(0,max(chipseq$sse)), xlim=c(0,max(chipseq$b.choice)))
identify(rnaseq$b.choice,rnaseq$sse,labels=rnaseq$sample)

legend(0.7,0.08,c("chipseq","rnaseq"),
         col=c("blue","red"), cex=2, pt.cex=2,
         text.col = "black", pch = 15, bg = 'white')

# dev.copy2pdf(file = paste(filename,"-check-",minN,"-",maxN,".pdf", sep=""))
# dev.off()
# 
# ####################################################################################
# ## correlation of size of dataset (# reads) vs b.choice
# x11(width=15,height=9)
# par(cex.axis=1.5, cex.lab=2, cex.main=2, mar=c(5,5,5,5))
# plot(data$b.compare[data$assay == "chipseq"],data$datasetSize[data$assay == "chipseq"],col="blue",
#      xlab="b parameter", ylab="datasetSize", ylim=c(0,max(data$datasetSize)), xlim=c(0,max(data$b.compare)))
# par(new=TRUE)
# plot(data$b.compare[data$assay == "rnaseq"],data$datasetSize[data$assay == "rnaseq"],col="red",
#      bty='n',ylab='',xlab='',yaxt='n',xaxt='n',
#      ylim=c(0,max(data$datasetSize)), xlim=c(0,max(data$b.compare)))
# 
# c = cor(data$b.compare,data$datasetSize, method="spearman")
# c.c = cor(data$b.compare[data$assay == "chipseq"],data$datasetSize[data$assay == "chipseq"], method="spearman")
# c.r = cor(data$b.compare[data$assay == "rnaseq"],data$datasetSize[data$assay == "rnaseq"], method="spearman")
# text(0.6,5e09,bquote("Spearman.all," * rho * "=" * .(signif(c,3))))
# text(0.6,4.7e09,bquote("Spearman.chipseq," * rho * "=" * .(signif(c.c,3))))
# text(0.6,4.4e09,bquote("Spearman.rnaseq," * rho * "=" * .(signif(c.r,3))))
# identify(data$b.compare,data$datasetSize,labels=data$sample)
# 
# x11(width=15,height=9)
# par(cex.axis=1.5, cex.lab=2, cex.main=2, mar=c(5,5,5,5))
# plot(data$b.grad[data$assay == "chipseq"],data$datasetSize[data$assay == "chipseq"],col="blue",
#      xlab="b parameter", ylab="datasetSize", ylim=c(0,max(data$datasetSize)), xlim=c(0,max(data$b.grad)))
# par(new=TRUE)
# plot(data$b.grad[data$assay == "rnaseq"],data$datasetSize[data$assay == "rnaseq"],col="red",
#      bty='n',ylab='',xlab='',yaxt='n',xaxt='n',
#      ylim=c(0,max(data$datasetSize)), xlim=c(0,max(data$b.grad)))
# 
# c = cor(data$b.grad,data$datasetSize, method="spearman")
# c.c = cor(data$b.grad[data$assay == "chipseq"],data$datasetSize[data$assay == "chipseq"], method="spearman")
# c.r = cor(data$b.grad[data$assay == "rnaseq"],data$datasetSize[data$assay == "rnaseq"], method="spearman")
# text(0.6,5e09,bquote("Spearman.all," * rho * "=" * .(signif(c,3))))
# text(0.6,4.7e09,bquote("Spearman.chipseq," * rho * "=" * .(signif(c.c,3))))
# text(0.6,4.4e09,bquote("Spearman.rnaseq," * rho * "=" * .(signif(c.r,3))))
