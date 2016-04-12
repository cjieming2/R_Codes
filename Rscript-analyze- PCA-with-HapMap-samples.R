
###############################################
###############################################
##                                           ##                                    
## Rscript-analyze- PCA-with-HapMap-samples  ##
## Date: March 11, 2008                      ##
## Version 1                                 ##   
##                                           ##   
###############################################
###############################################

library(GDD)

##Read the pca matrix
setwd("C:/Documents and Settings/linx/Desktop/BCH-826samples-550K/PCA/PCA-with-HapMap")
evector=read.table("BCH-HapMap-1014samples-492914snps.pca",header=T)
dim(evector)
#First column of evector gives "sample_id"
evector[1:10,1:10]


##Read Phenotype data 
#HapMap phenotype file has 208 samples, but PCA is only on 206 samples)
#BCH phenotype file has 826 samples, but PCA is only on 808 samples

phenotype1=read.table("208samples-hapmap.sa",header=TRUE)
dim(phenotype1)
head(phenotype1)

setwd("C:/Documents and Settings/linx/Desktop/BCH-826samples-550K/pheno")
phenotype2=read.table("BCH-826samples.sa",header=TRUE)
dim(phenotype2)
head(phenotype2)

phenotype3=cbind(phenotype2,rep("FIN",nrow(phenotype2)))[,c(1,4)]
colnames(phenotype3)[2]="population.id"
dim(phenotype3)
head(phenotype3)
rm(phenotype2)

phenotype=rbind(phenotype3,phenotype1)
dim(phenotype)
head(phenotype)
rm(phenotype3,phenotype1)

setwd("C:/Documents and Settings/linx/Desktop/BCH-826samples-550K/PCA/PCA-with-HapMap")

##Combine phenotype file with evector, create a new outlier column
data_temp=merge(phenotype, evector, by="sample.id",sort=F)
dim(data_temp)
outlier=rep(21, times=nrow(evector))
length(outlier)
head(outlier)
data=cbind(data_temp,outlier)
dim(data)
data[1:5,1:10]
rm(data_temp,outlier)



##Plots
col = c(  "red2","royalblue","purple", "turquoise","springgreen3")


#PC1 vs PC3
pca13.id <- which(data$population.id=="FIN"&data$PC3 > 0.065)
data[pca13.id, "outlier"] <- 3
data[pca13.id,1]

x11(10,7, pointsize = 18)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
layout.show(2)
plot(data[,ncol(phenotype)+1], data[,ncol(phenotype)+3], xlab="First principal component",
     ylab="Third principal component", pch=data$outlier, bg=col[unclass(data$population.id)],main="PCA-BCH-808samples-HapMap-206samples-492914snps")
text(x=data[pca13.id[1],ncol(phenotype)+1]-0.003, y=data[pca13.id[1],ncol(phenotype)+3], labels=data[pca13.id[1],1], cex=0.6, col=col[unclass(data[pca13.id[1],"population.id"])])
text(x=data[pca13.id[2],ncol(phenotype)+1]-0.003, y=data[pca13.id[2],ncol(phenotype)+3], labels=data[pca13.id[2],1], cex=0.6, col=col[unclass(data[pca13.id[2],"population.id"])])
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(data$population.id)), pt.bg = col,text.col = "black", pch = 21,bg = 'white')


#PC1 vs PC2

x11(10,7, pointsize = 18)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
layout.show(2)
plot(data[,ncol(phenotype)+1], data[,ncol(phenotype)+2], xlab="First principal component",
     ylab="Second principal component", pch=data$outlier, bg=col[unclass(data$population.id)],main="PCA-BCH-808samples-HapMap-206samples-492914snps")
text(x=data[pca13.id[1],ncol(phenotype)+1]-0.003, y=data[pca13.id[1],ncol(phenotype)+2], labels=data[pca13.id[1],1], cex=0.6, col=col[unclass(data[pca13.id[1],"population.id"])])
text(x=data[pca13.id[2],ncol(phenotype)+1]-0.003, y=data[pca13.id[2],ncol(phenotype)+2], labels=data[pca13.id[2],1], cex=0.6, col=col[unclass(data[pca13.id[2],"population.id"])])
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(data$population.id)), pt.bg = col,text.col = "black", pch = 21,bg = 'white')


#PC1 vs PC4

x11(10,7, pointsize = 18)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
layout.show(2)
plot(data[,ncol(phenotype)+1], data[,ncol(phenotype)+4], xlab="First principal component",
     ylab="Fourth principal component", pch=data$outlier, bg=col[unclass(data$population.id)],main="PCA-BCH-808samples-HapMap-206samples-492914snps")
text(x=data[pca13.id[1],ncol(phenotype)+1]-0.003, y=data[pca13.id[1],ncol(phenotype)+4], labels=data[pca13.id[1],1], cex=0.6, col=col[unclass(data[pca13.id[1],"population.id"])])
text(x=data[pca13.id[2],ncol(phenotype)+1]-0.003, y=data[pca13.id[2],ncol(phenotype)+4], labels=data[pca13.id[2],1], cex=0.6, col=col[unclass(data[pca13.id[2],"population.id"])])
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(data$population.id)), pt.bg = col,text.col = "black", pch = 21,bg = 'white')


#PC1 vs PC5

x11(10,7, pointsize = 18)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
layout.show(2)
plot(data[,ncol(phenotype)+1], data[,ncol(phenotype)+5], xlab="First principal component",
     ylab="Fifth principal component", pch=data$outlier, bg=col[unclass(data$population.id)],main="PCA-BCH-808samples-HapMap-206samples-492914snps")
text(x=data[pca13.id[1],ncol(phenotype)+1]-0.003, y=data[pca13.id[1],ncol(phenotype)+5], labels=data[pca13.id[1],1], cex=0.6, col=col[unclass(data[pca13.id[1],"population.id"])])
text(x=data[pca13.id[2],ncol(phenotype)+1]-0.003, y=data[pca13.id[2],ncol(phenotype)+5], labels=data[pca13.id[2],1], cex=0.6, col=col[unclass(data[pca13.id[2],"population.id"])])
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(data$population.id)), pt.bg = col,text.col = "black", pch = 21,bg = 'white')


#PC1 vs PC6

x11(10,7, pointsize = 18)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
layout.show(2)
plot(data[,ncol(phenotype)+1], data[,ncol(phenotype)+6], xlab="First principal component",
     ylab="Sixth principal component", pch=data$outlier, bg=col[unclass(data$population.id)],main="PCA-BCH-808samples-HapMap-206samples-492914snps")
text(x=data[pca13.id[1],ncol(phenotype)+1]-0.003, y=data[pca13.id[1],ncol(phenotype)+6], labels=data[pca13.id[1],1], cex=0.6, col=col[unclass(data[pca13.id[1],"population.id"])])
text(x=data[pca13.id[2],ncol(phenotype)+1]-0.003, y=data[pca13.id[2],ncol(phenotype)+6], labels=data[pca13.id[2],1], cex=0.6, col=col[unclass(data[pca13.id[2],"population.id"])])
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(data$population.id)), pt.bg = col,text.col = "black", pch = 21,bg = 'white')



#PC1 vs PC7

pca17.id <- which(data$population.id=="FIN"&data$PC7 > 0.2)
data[pca17.id, "outlier"] <- 3
data[pca17.id,1]

x11(10,7, pointsize = 18)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
layout.show(2)
plot(data[,ncol(phenotype)+1], data[,ncol(phenotype)+7], xlab="First principal component",
     ylab="Seventh principal component", pch=data$outlier, bg=col[unclass(data$population.id)],main="PCA-BCH-808samples-HapMap-206samples-492914snps")
text(x=data[pca17.id,ncol(phenotype)+1]-0.003, y=data[pca17.id,ncol(phenotype)+7], labels=data[pca17.id,1], cex=0.6, col=col[unclass(data[pca17.id,"population.id"])])
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(data$population.id)), pt.bg = col,text.col = "black", pch = 21,bg = 'white')


#PC1 vs PC8

pca18.id <- which(data$population.id=="FIN"&data$PC8 > 0.11)
data[pca18.id, "outlier"] <- 3
data[pca18.id,1]

pca18b.id <- which(data$population.id=="FIN"&data$PC8 < -0.1)
data[pca18b.id, "outlier"] <- 3
data[pca18b.id,1]

x11(10,7, pointsize = 18)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
layout.show(2)
plot(data[,ncol(phenotype)+1], data[,ncol(phenotype)+8], xlab="First principal component",
     ylab="Eighth principal component", pch=data$outlier, bg=col[unclass(data$population.id)],main="PCA-BCH-808samples-HapMap-206samples-492914snps")
text(x=data[pca18.id[c(1,2,4,5)],ncol(phenotype)+1]-0.003, y=data[pca18.id[c(1,2,4,5)],ncol(phenotype)+8], labels=data[pca18.id[c(1,2,4,5)],1], cex=0.6, col=col[unclass(data[pca18.id[c(1,2,4,5)],"population.id"])])
text(x=data[pca18.id[c(3)],ncol(phenotype)+1]-0.003, y=data[pca18.id[c(3)],ncol(phenotype)+8]-0.007, labels=data[pca18.id[c(3)],1], cex=0.6, col=col[unclass(data[pca18.id[c(3)],"population.id"])])
text(x=data[pca18b.id[c(2,3,4,6)],ncol(phenotype)+1]-0.003, y=data[pca18b.id[c(2,3,4,6)],ncol(phenotype)+8], labels=data[pca18b.id[c(2,3,4,6)],1], cex=0.6, col=col[unclass(data[pca18b.id[c(2,3,4,6)],"population.id"])])
text(x=data[pca18b.id[c(1,5)],ncol(phenotype)+1]+0.003, y=data[pca18b.id[c(1,5)],ncol(phenotype)+8], labels=data[pca18b.id[c(1,5)],1], cex=0.6, col=col[unclass(data[pca18b.id[c(1,5)],"population.id"])])
text(x=data[pca13.id[1],ncol(phenotype)+1]-0.003, y=data[pca13.id[1],ncol(phenotype)+8], labels=data[pca13.id[1],1], cex=0.6, col=col[unclass(data[pca13.id[1],"population.id"])])
text(x=data[pca13.id[2],ncol(phenotype)+1]-0.003, y=data[pca13.id[2],ncol(phenotype)+8], labels=data[pca13.id[2],1], cex=0.6, col=col[unclass(data[pca13.id[2],"population.id"])])
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(data$population.id)), pt.bg = col,text.col = "black", pch = 21,bg = 'white')


#PC1 vs PC9

pca19.id <- which(data$population.id=="FIN"&data$PC9 > 0.27)
data[pca19.id, "outlier"] <- 3
data[pca19.id,1]

x11(10,7, pointsize = 18)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
layout.show(2)
plot(data[,ncol(phenotype)+1], data[,ncol(phenotype)+9], xlab="First principal component",
     ylab="Nineth principal component", pch=data$outlier, bg=col[unclass(data$population.id)],main="PCA-BCH-808samples-HapMap-206samples-492914snps")
text(x=data[pca19.id,ncol(phenotype)+1]-0.003, y=data[pca19.id,ncol(phenotype)+9], labels=data[pca19.id,1], cex=0.6, col=col[unclass(data[pca19.id,"population.id"])])
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(data$population.id)), pt.bg = col,text.col = "black", pch = 21,bg = 'white')


#PC1 vs PC10

pca110.id <- which(data$population.id=="FIN"&data$PC10 < -0.3)
data[pca110.id, "outlier"] <- 3
data[pca110.id,1]

x11(10,7, pointsize = 18)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
layout.show(2)
plot(data[,ncol(phenotype)+1], data[,ncol(phenotype)+10], xlab="First principal component",
     ylab="Tenth principal component", pch=data$outlier, bg=col[unclass(data$population.id)],main="PCA-BCH-808samples-HapMap-206samples-492914snps")
text(x=data[pca110.id,ncol(phenotype)+1]-0.003, y=data[pca110.id,ncol(phenotype)+10], labels=data[pca110.id,1], cex=0.6, col=col[unclass(data[pca110.id,"population.id"])])
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(data$population.id)), pt.bg = col,text.col = "black", pch = 21,bg = 'white')


#Scatter plot of first 10 Principal components

windows()
GDD(file=paste("PCA-BCH-808samples-HapMap-206samples-492914snps-PC1-PC5.png"),height=800,width=800,ps=18)
pairs(data[,(ncol(phenotype)+1):(ncol(phenotype)+5)], pch=data$outlier, col=col[unclass(data$population.id)], main="PCA-BCH-808samples-HapMap-206samples-492914snps")
dev.off()

windows()
GDD(file=paste("PCA-BCH-808samples-HapMap-206samples-492914snps-PC6-PC10.png"),height=800,width=800,ps=18)
pairs(data[,(ncol(phenotype)+6):(ncol(phenotype)+10)], pch=data$outlier, col=col[unclass(data$population.id)], main="PCA-BCH-808samples-HapMap-206samples-492914snps")
dev.off()

windows()
GDD(file=paste("PCA-BCH-808samples-HapMap-206samples-492914snps-PC1-PC10.png"),height=800,width=800,ps=18)
pairs(data[,(ncol(phenotype)+1):(ncol(phenotype)+10)], pch=data$outlier, col=col[unclass(data$population.id)], main="PCA-BCH-808samples-HapMap-206samples-492914snps")
dev.off()
