setwd('C:/Users/JM/yale/courses/courses spring 2011/stats 660 multivariate stats in social and env sciences/our dataset')
#source("C:/Shared/scripts-R_perl_shell_macros/R codes/hartiganROT.R")
# source("C:/Shared/scripts-R_perl_shell_macros/R codes/HartiganROT1.R")

## input data
file = "E-GEOD-8536-processed-data-1631482821.txt"
data <- read.table(file, header=TRUE, row.names=1, sep="\t")
dataname = "Yeast expression data"

## transform data
#data <- sqrt(data)

## doing faces
library(TeachingDemos)
x11()
faces(data,byrow=TRUE,main=dataname)

## star plot - to many observations
#stars(t(data),scale=TRUE,draw.segments=TRUE,col.segments=TRUE, main=dataname)

## hierarchical clustering on variables, euclidean/manhattan via single/ward/complete (dist)
require(graphics)
distmethod = "euclidean"
clusmethod = "single"
hc1 <- hclust(dist(t(data),method=distmethod),method=clusmethod)
x11();plot(hc1, main=paste(dataname,"Clustering method=",clusmethod,"via Distance measure=",distmethod))
x11();plot(hc1, hang=-1,main=paste(dataname,"Clustering method=",clusmethod,"via Distance measure=",distmethod))

distmethod = "manhattan"
clusmethod = "single"
hc2 <- hclust(dist(t(data),method=distmethod),method=clusmethod)
x11();plot(hc2,main=paste(dataname,"Clustering method=",clusmethod,"via Distance measure=",distmethod))
x11();plot(hc2, hang=-1,main=paste(dataname,"Clustering method=",clusmethod,"via Distance measure=",distmethod))

distmethod = "euclidean"
clusmethod = "ward"
hc3 <- hclust(dist(t(data),method=distmethod),method=clusmethod)
x11();plot(hc3,main=paste(dataname,"Clustering method=",clusmethod,"via Distance measure=",distmethod))
x11();plot(hc3, hang=-1,main=paste(dataname,"Clustering method=",clusmethod,"via Distance measure=",distmethod))

distmethod = "manhattan"
clusmethod = "ward"
hc4 <- hclust(dist(t(data),method=distmethod),method=clusmethod)
x11();plot(hc4,main=paste(dataname,"Clustering method=",clusmethod,"via Distance measure=",distmethod))
x11();plot(hc4, hang=-1,main=paste(dataname,"Clustering method=",clusmethod,"via Distance measure=",distmethod))

distmethod = "euclidean"
clusmethod = "complete"
hc5 <- hclust(dist(t(data),method=distmethod),method=clusmethod)
x11();plot(hc5,main=paste(dataname,"Clustering method=",clusmethod,"via Distance measure=",distmethod))
x11();plot(hc5, hang=-1,main=paste(dataname,"Clustering method=",clusmethod,"via Distance measure=",distmethod))

distmethod = "manhattan"
clusmethod = "complete"
hc6 <- hclust(dist(t(data),method=distmethod),method=clusmethod)
x11();plot(hc6,main=paste(dataname,"Clustering method=",clusmethod,"via Distance measure=",distmethod))
x11();plot(hc6, hang=-1,main=paste(dataname,"Clustering method=",clusmethod,"via Distance measure=",distmethod))

## kmeans
kclust3 <- kmeans(t(data),centers=3)
kclust5 <- kmeans(t(data),centers=5)
# # x11();plot(data, col = kclust$cluster)
# # hartROT <- harti.rule(t(data),have.rownames=1)
# tdata = t(data)
# tdatanew <- data.frame(rownames(tdata),scale(tdata[,2:14]))
# hartROT <- harti.rule(tdatanew)
#x11()
#plot(hartROT["SSwithin",],type="b",col="red",lwd=2,pch=19,cex=1.2,
#          xlab="# clusters", ylab="SSwithin",
#          main=paste("Scree Plot of SSwithin results against the # clusters in the",dataname))

################DARRYL#####################################
calcHartiganROT <- function(ssK, ssKplus1, n, k){
  ((ssK/ssKplus1)-1)*(n-k-1)
}

rowCount <- nrow(tdata)
hartiganROT <- c()
comparisons <- c()
ss <- c()
kkk <- c()
for(k in 2:rowCount-1) {
  clusters <- kmeans(tdata, centers=k) 
  ssK <- clusters$tot.withinss
  ss <- c(ss, ssK)
  kkk <- c(kkk, k)
  if(k > 2) {
    rot <- calcHartiganROT(ssKMinus1, ssK, dim(tdata)[1], k)
    hartiganROT <- c(hartiganROT, rot)
    comparisons <- c(comparisons, paste(k-1, 'v', k))
  }
  ssKMinus1 <- ssK
}
hartiganRatios <- data.frame(hartiganROT=hartiganROT, row.names=comparisons)
x11(); plot(kkk,ss,type="b",col="red",lwd=2,pch=19,cex=1.2,
         xlab="# clusters", ylab="SSwithin",
         main=paste("Scree Plot of SSwithin results against the # clusters in the",dataname))