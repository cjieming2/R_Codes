setwd("/home/dr395/STAT660")
yeastExp <- read.delim("E-GEOD-8536-processed-data-1631482821.txt", header=T)
yeastData <- as.data.frame(yeastExp[,-1]) 
setwd("/home/dr395/STAT660/hw3")
rownames(yeastData) <- yeastExp[,1]

# sort genes by expression variance across samples (genes are rows)
vars <- apply(yeastData, 1, function(x) var(x))
# add fold change to dataframe
yeastData.sorted = cbind(yeastData, vars)
# sort dataframe by fold change
yeastData.sorted <- yeastData.sorted[order(yeastData.sorted[,ncol(yeastData.sorted)]),]
# remove variance measure from dataframe and select first 5000 genes for analysis
yeastData0 <- yeastData.sorted[1:4000,-ncol(yeastData.sorted)]
yeastData <- yeastData0


pdf("euclidean_complete.dendrogram.pdf")
dist <- dist(yeastData, method="euclidean")
clust <- hclust(dist, method="complete")
plot(clust, labels=rownames(yeastData),cex=0.5,xlab="",ylab="Distance" ,main="Clustering for Genes - Euclidean/Complete")
#rect.clust(clust1, k=2)
dev.off()
pdf("euclidean_average.dendrogram.pdf")
dist <- dist(yeastData, method="euclidean")
clust <- hclust(dist, method="average")
plot(clust, labels=rownames(yeastData),cex=0.5,xlab="",ylab="Distance",main="Clustering for Genes - Euclidean/Average")
#rect.hclust(clust2, k=2)
dev.off()
pdf("euclidean_ward.dendrogram.pdf")
dist <- dist(yeastData, method="euclidean")
clust <- hclust(dist, method="ward")
plot(clust, labels=rownames(yeastData),cex=0.5,xlab="",ylab="Distance",main="Clustering for Genes - Euclidean/Ward")
rect.hclust(clust,k=20)
cuts = cutree(clust,k=20)
dev.off()

save(yeastData, dist, clust, cuts, file="Preserved_Data.RData")

# Read in probe annotation
annot = read.delim(file="A-AFFY-27.noheader.adf.txt", header=FALSE)

# Match probes in the data set to the probe IDs
probes = rownames(yeastData)
probes2annot = match(probes, annot$V1)

# Get corresponding Gene IDs
allGeneIDs = annot$V4[probes2annot]


#List genes
for (i in 1:20) 
{
        # select cluster probes
        modGenes = (cuts == i)
        # get their entrez ID codes
        modGeneIDs = allGeneIDs[modGenes]
        
        # write them into a file
        fileName = paste("cluster_", i, ".txt", sep="")
        write.table(as.data.frame(modGeneIDs[modGeneIDs!='']), file=fileName, row.names=FALSE, col.names=FALSE)
}


for (i in 1:6) {
  write(paste("Genes in cluster ", i), file="clusters.txt", append=T)
  write(rownames(yeastData)[cuts==i], file="clusters.txt", append=T)
  print(" ", file="clusters.txt", append=T)
}

pdf("manhattan_complete.dendrogram.pdf")
dist <- dist(yeastData, method="manhattan")
clust <- hclust(dist, method="complete")
plot(clust, labels=rownames(yeastData),cex=0.5,xlab="",ylab="Distance",main="Clustering for Genes - Manhattan/Complete")
#rect.hclust(clust4, k=2)
dev.off()
pdf("manhattan_average.dendrogram.pdf")
dist <- dist(yeastData, method="manhattan")
clust <- hclust(dist, method="average")
plot(clust, labels=rownames(yeastData),cex=0.5,xlab="",ylab="Distance",main="Clustering for Genes - Manhattan/Average")
#rect.hclust(clust4, k=2)
dev.off()
pdf("manhattan_ward.dendrogram.pdf")
dist <- dist(yeastData, method="manhattan")
clust <- hclust(dist, method="ward")
plot(clust, labels=rownames(yeastData),cex=0.5,xlab="",ylab="Distance",main="Clustering for Genes - Manhattan/Ward")
#rect.hclust(clust4, k=2)
dev.off()
calcHartiganROT <- function(ssK, ssKplus1, n, k){
  ((ssK/ssKplus1)-1)*(n-k-1)
}

numGenes <- nrow(yeastData)
hartiganROT <- c()
comparisons <- c()
ss <- c()
for(k in 2:25) {
  clusters <- kmeans(yeastData, centers=k)
  #pdf(paste("k", k, sep=''))
  #plot(clusters, labels=votes[,1],cex=0.5,xlab="",ylab="Distance",main=paste("Clustering for k =", k))
  ssK <- clusters$tot.withinss
  ss <- c(ss, ssK)
  if(k > 2) {
    rot <- calcHartiganROT(ssKMinus1, ssK, numGenes, k)
    hartiganROT <- c(hartiganROT, rot)
    comparisons <- c(comparisons, paste(k-1, 'v', k))
  }
  ssKMinus1 <- ssK
}
hartiganRatios <- data.frame(hartiganROT=hartiganROT, row.names=comparisons)
        hartiganROT
2 v 3   6637.817658
3 v 4   3490.238024
4 v 5   2033.065699
5 v 6   1184.394890
6 v 7    877.119234
7 v 8    597.956125
8 v 9    541.852593
9 v 10   521.766662
10 v 11  445.193349
11 v 12  345.357924
12 v 13  400.101803
13 v 14  148.439676
14 v 15  398.637469
15 v 16  256.643159
16 v 17  172.724718
17 v 18  324.969107
18 v 19  238.076586
19 v 20  117.994151
20 v 21  225.246011
21 v 22  176.318896
22 v 23  214.810349
23 v 24  196.496407
24 v 25    5.819798




