## takes a matrix or columns of values and find correlation between the columns

setwd("C:/Documents and Settings/chenjm/Desktop/")

data <- read.table("correlation", header=TRUE, sep="\t")

## Spearman's rho  and  kendall's tau
cS <- cor(data, method = "spearman")
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
