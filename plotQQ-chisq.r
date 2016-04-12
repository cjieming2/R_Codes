setwd("C:/Documents and Settings/chenjm/Desktop/work_documentation/JJ/chinese population analysis/manuscript/ajhg submission/review1/simulation/shanghai/502626snps-someoutliersrm")
#setwd("C:/Documents and Settings/chenjm/Desktop/work_documentation/JJ/chinese population analysis/simulation analysis/simu2-chisq-removeSigni-removeConserved")


filename = "eigen-output.chisq.GC-noNA.txt"
data <- read.table(filename, header=TRUE, sep="\t")

## Doing the qq plot manually
## please remove missing data
## chisq values
size <- dim(data)[1]
p <- ((1:size)-0.5)/size

## expected quantiles
chisq.quantile <- qchisq(p,df=1)

## observed quantiles
obschisq.quantile <- sort(data$Chisq)
obschisqgc.quantile <- sort(data$Chisq.gc)
obseigen.quantile <- sort(data$EIGENSTRAT)
obseigengc.quantile <- sort(data$EIGENSTRAT.gc)

## for labelling
sortedchi <- data[order(data$Chisq, decreasing=FALSE),] ## order works on the whole matrix
sortedchigc <- data[order(data$Chisq.gc, decreasing=FALSE),] ## order works on the whole matrix
sortedeigen <- data[order(data$EIGENSTRAT, decreasing=FALSE),] ## order works on the whole matrix
sortedeigengc <- data[order(data$EIGENSTRAT.gc, decreasing=FALSE),] ## order works on the whole matrix

## plotting
x11()
plot(chisq.quantile, obschisq.quantile,
#plot(obschisq.quantile, chisq.quantile,
	xlab="Expected quantiles", ylab="Observed quantiles", main=paste("chisq qqplot: ",filename,sep=""))
abline(a=0, b=1) # intercept,a = 0,0 b is slope x=y b=1
## identifying
#id <- which (obschisq.quantile > 20)
#points(chisq.quantile[id], obschisq.quantile[id], col="red")
#text(chisq.quantile[id], obschisq.quantile[id], labels=sortedchi$SNP_ID[id], pos=2, offset=0.5, col="red")

x11()
plot(chisq.quantile, obschisqgc.quantile,
	xlab="Expected quantiles", ylab="Observed quantiles", main=paste("chisq qqplot+gc: ",filename,sep=""))
abline(a=0, b=1) # intercept,a = 0,0 b is slope x=y b=1
## identifying
#id <- which (obschisqgc.quantile > 20)
#points(chisq.quantile[id], obschisqgc.quantile[id], col="red")
#text(chisq.quantile[id], obschisqgc.quantile[id], labels=sortedchigc$SNP_ID[id], pos=2, offset=0.5, col="red")

x11()
plot(chisq.quantile, obseigen.quantile,
	xlab="Expected quantiles", ylab="Observed quantiles", main=paste("chisq qqplot-eigen: ",filename,sep=""))
abline(a=0, b=1) # intercept,a = 0,0 b is slope x=y b=1
## identifying
#id <- which (obseigen.quantile > 20)
#points(chisq.quantile[id], obseigen.quantile[id], col="red")
#text(chisq.quantile[id], obseigen.quantile[id], labels=sortedeigen$SNP_ID[id], pos=2, offset=0.5, col="red")

x11()
plot(chisq.quantile, obseigengc.quantile,
	xlab="Expected quantiles", ylab="Observed quantiles", main=paste("chisq qqplot-eigen+gc: ",filename,sep=""))
abline(a=0, b=1) # intercept,a = 0,0 b is slope x=y b=1
## identifying
#id <- which (sortedeigengc$SNP_ID == "rs12141395" | sortedeigengc$SNP_ID == "rs7641148" | sortedeigengc$SNP_ID == "rs4980598")
#id <- which (obseigengc.quantile > 20)
#points(chisq.quantile[id], obseigengc.quantile[id], col="red")
#points(chisq.quantile[id2], obseigengc.quantile[id2], col="blue")
#text(chisq.quantile[id], obseigengc.quantile[id], labels=sortedeigengc$SNP_ID[id], pos=2, offset=0.5, col="red")
#text(chisq.quantile[id2], obseigengc.quantile[id2], labels=sortedeigengc$SNP_ID[id2], pos=2, offset=0.5, col="blue")