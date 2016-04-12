setwd('C:/Users/Jieming/Documents/thesis/lynne_work/motifVar/tpr/exac')

## relative entropy data
data = read.table('final.merged.ExAC.ens73.TPR.34aa.vepoutput.coding.canonical.sorted.auto.txt', sep="\t", header=T, stringsAsFactors = F)

## distribution of read depth V7
bins=50
a = hist(data$DP, breaks=bins, plot=F)

x11()
barplot(a$counts/sum(a$counts), main="RD dist for SNVs in TPR.34aa", names.arg=a$mid, ylim=c(0,0.4), las=2)
