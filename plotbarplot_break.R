library(plotrix)
setwd("C:/Users/JM/thesis/lynne_work/3TPR predicting function and affinity with binders/motifNum");

filename1 = '10599proteins-34aa-smart7_111128_removedByProts_distrib_forR.txt';

data <- read.table(filename1,sep="\t");

x11()
gap.barplot(data[,2],gap=c(2001,3500),ytics=c(seq(0,2000,500),4000,4100),
            xlab="Number of TPR repeats",ylab="Number of proteins",
            main="Number of TPR repeats in SMART database",
            cex.main=1.5,cex.axis=1.5,cex.lab=2, xlim=c(1,30));