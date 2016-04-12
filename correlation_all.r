## takes a matrix or columns of values and find correlation between the columns

#setwd("C:/Users/JM/thesis/mark_work/jasmine_ncvar/plots")
#setwd("C:/Users/JM/Desktop")
setwd("C:/Users/JM/thesis/mark_work/ekta_netsnps/correlations_1KGp1/snpdensity")

# data <- read.table("combined-523genes-1KG-p1-hg19-101123-ensemblID-grch37-signalink2009-fullist-hsaOnly.txt", header=TRUE, sep="\t")
# data <- read.table("dnds_af",header=TRUE)
data <- read.table("combined-no0deg-523genes-1KG-p1-hg19-101123-ensemblID-grch37-signalink2009-fullist-hsaOnly-S.txt",header=TRUE)
col2name = "total"
col1name = "snpdensity"
col2 <- data$total
col1 <- data$snpdensity
#col2= data[,2]*(1-data[,2])*2

## plots
# x11()
# plot(col1,col2,xlab=paste(col1name),ylab=paste(col2name))
# # par(new=T)
# x11()
# boxplot(col2 ~ col1,xlab=paste(col1name),ylab=paste(col2name),
#         main=paste("Boxplot of",col2name,"VS",col1name))

# 
# x11()
# hist(col1,freq=FALSE,xlab=paste(col1name), main=paste("Histogram of",col1name))

# x11()
# hist(col2,freq=FALSE,xlab=paste(col2name), main=paste("Histogram of",col2name))

#PEARSON
#cor(col1,col2,method="pearson")
#cor.test(col1,col2,alternative="two.sided",method="pearson")
#cor.test(col1,col2,alternative="less",method="pearson")
#cor.test(col1,col2,alternative="greater",method="pearson")

#SPEARMAN
#cor(col1,col2,method="spearman")
cor.test(col1,col2,alternative="two.sided",method="spearman",exact=TRUE)
# cor.test(col1,col2,alternative="less",method="spearman",exact=TRUE)
# cor.test(col1,col2,alternative="greater",method="spearman",exact=TRUE)

#KENDALL
# cor(col1,col2,method="kendall")
# cor.test(col1,col2,alternative="two.sided",method="kendall",exact=NULL)
# cor.test(col1,col2,alternative="less",method="kendall",exact=NULL)
# cor.test(col1,col2,alternative="greater",method="kendall",exact=NULL)

################################################################################################################
######### builds a correlation matrix and calculate the individual correlation 
######### based on Spearman's or Kendall's
# data <- read.table("son-parents.txt", header=TRUE, sep="\t")

## Spearman's rho  and  kendall's tau
# cS <- cor(data, method = "spearman")
#cK <- cor(data, method = "kendall")

# fileS = "correlation.S"
#fileK = "correlation.K"
#write.table(cK, file=fileK, quote=F, row.names=T, sep="\t")

## symnum puts up a symbolic correlation matrix
# jm <- symnum(cS, corr=TRUE, cutpoints=c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95), symbols=c("", ">0.5", ">0.6", ">0.7", ">0.8", ">0.9", ">0.95"))
#symnum(cP)

## write to files
# write.table(cS, file=fileS, quote=F, row.names=T, sep="\t")
# write.table(jm, file="correlationSymbol.S", quote=F, row.names=T, sep="\t")

## How much do they differ?
#i <- lower.tri(cS)
#j <- lower.tri(cP)
#cor(cbind(S = cS[i], K = cK[i]))

## rcorr command from Hmisc package
#x <- c(-2, -1, 0, 1, 2)
#y <- c(4,   1, 0, 1, 4)
#z <- c(1,   2, 3, 4, NA)
#v <- c(1,   2, 3, 4, 5)

#rcorr(cbind(longley$GNP.deflator,longley$GNP),type=c("spearman"))
#end
