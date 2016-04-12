setwd("C:/Users/Jieming/Documents/thesis/lynne_work/motifVar/tpr/blosum_sift_polyphen2_gerp")

library(gplots)

##################################################################################################
#### setup
filename = 'final.merged.ExAC.ens73.TPR.34aa.vepoutput.coding.canonical.sorted.auto.blosum62-30-35-40-45-50.sift.gerp.bed'
data = read.table(filename, header=F, stringsAsFactors = F, comment.char="", sep="\t")
data.ns.only <- data[!(data$V56 == "SYNONYMOUS_CODING" | data$V56 == "SYNONYMOUS_CODING,SPLICE_SITE"),]
resNumCol=45
aachangeCol=64
aachange = data.frame(data.ns.only[,resNumCol],data.ns.only[,aachangeCol])
names(aachange) = c("resNum","aachange")

## split -> so that you get aa1 in col3 and aa2 in col4
## then col bind them back
df = data.frame(do.call('rbind', strsplit(as.character(aachange$aachange),"->",fixed=T)))
aachange = cbind(aachange,df)
names(aachange) = c("resNum","aachange","aa1","aa2")
mymat = as.matrix(table(aachange$aa1,aachange$aa2))

x11()
colors = heat.colors(100)
heatmap.2(mymat,col=rev(colors))

