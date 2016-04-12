setwd("T:/HHPS/326samples-19580snps-hhps-han/20080424-pca")

#annotation loading
sample_annotation <- read.table("326samples-hhps-han.sa", header=T, sep = "\t")

##########
#HHPS-HAN#
##########
data <- read.table("326samples-19580snps-hhps-han.pca", header=T, sep = "\t")
data <- merge(sample_annotation, data, by = "sample.id")


x11(10,7, pointsize = 18)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
layout.show(2)
#col = rainbow(29)
col = c("red", "orange", "pink", "pink3", "lightblue1", "purple", "lightblue2", "lightblue3")
col = "lightblue3";
plot(data$PC1, data$PC2, xlab="PC1", ylab="PC2", type = "n", main="Han Chinese (326 Individuals, 19580 SNPs)", pch = 21, bg = col[unclass(data$population.id)])
subdata <- subset(data, population.id=="TW-HB")
points(subdata$PC1, subdata$PC2, pch = 21, bg = col[unclass(factor(subdata$population.id))])
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(subdata$population.id)), pt.bg = col,
        text.col = "black", pch = 21,
        bg = 'white')



##################
#SNP Correlations#
##################
snp_correlation <- read.table("1043samples-660918snps-hgdp.cor", header=T, sep = "\t")
snp_annotation <- read.table("660918snps-hap650y-v1-db126.mk", header=T, sep = "\t")
snp_annotation <- snp_annotation[,1:3]
snp_annotation <- merge(snp_annotation, snp_correlation, by = "snp.id")

x11(10,7, pointsize = 12)
layout(matrix(c(1:24), 12, 2, byrow = F))
layout.show(24)
par(mar=c(0,4,0,0))

for (i in c(1:22, "X", "Y"))
{
snps <- subset(snp_annotation, chromosome == i, PC1!="nan")
if(dim(snps)[1]==0)
{
plot(0, 0, ylab = paste("", i), axes = F, xlim = c(0,250000000), ylim = c(0,1), col = "red", pch = ".")
#axis(1, labels = F, tick = T, tck = 0)
#axis(2, labels = F, tick = F)
}
else
{
plot(snps$position, abs(as.numeric(snps$PC1)), ylab = paste("", i), axes = F, xlim = c(0,250000000), ylim = c(0,1), col = "red", pch = ".")
#axis(1, labels = F, tick = T, tck = 0)
#axis(2, labels = F, tick = F)
}
}
