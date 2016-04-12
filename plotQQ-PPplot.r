#setwd("C:/Documents and Settings/chenjm/Desktop/work_documentation/JJ/chinese population analysis/simulation analysis/simu1-removeSigni-removeConserved")
setwd("C:/Documents and Settings/chenjm/Desktop/work_documentation/JJ/chinese population analysis/manuscript/ajhg submission/review1/simulation")


filename = "trimmed-100N-100N-new..assoc.tsv-noNA"
data <- read.table(filename, header=TRUE, sep="\t")

## Doing the qq plot manually
## chisq values
size <- dim(data)[1]
p <- ((1:size)-0.5)/size ## 0.5 is not fixed
chisq.quantile <- qchisq(p,df=1)
obschisq.quantile <- sort(data$CHISQ)
sortedchi <- data[order(data$CHISQ, decreasing=FALSE),] ## order works on the whole matrix
x11()
plot(chisq.quantile, obschisq.quantile,
	xlab="Expected quantiles", ylab="Observed quantiles", main=paste("chisq qqplot: ",filename,sep=""))
abline(a=0, b=1) # intercept,a = 0,0 b is slope x=y b=1
## identifying
#id <- which (obschisq.quantile > 20)
#points(chisq.quantile[id], obschisq.quantile[id], col="red")
#text(chisq.quantile[id], obschisq.quantile[id], labels=sortedchi$SNP[id], pos=2, offset=0.5, col="red")

## convert chisq to p
##1-pchisq(chisqstat, df)

## pp-plot - NEED TO REMOVE MONOMORPHIC SNPs
## pvalues
size <- dim(data)[1]
expected.quantile <- (-log10(qunif(seq(from=1,to=size,by=1)/(size+1))))
observed.quantile <- (-log10(sort(data$P))) ## sort works on 1 col
sortedP <- data[order(data$P, decreasing=FALSE),] ## order works on the whole matrix
##hist(observed.quantile)
x11()
plot(expected.quantile, observed.quantile, xlab="-log(expected P values)", 
	ylab="-log(observed P values)", main=paste("QQ plot p-values: ",filename, sep=""))
abline(a=0, b=1) # intercept,a = 0,0 b is slope x=y b=1

## identifying
#id <- which (observed.quantile > 5)
#points(expected.quantile[id], observed.quantile[id], col="red")
#text(expected.quantile[id], observed.quantile[id], labels=sortedP$SNP[id], pos=2, offset=0.5, col="red")