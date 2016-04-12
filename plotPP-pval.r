setwd("C:/Documents and Settings/chenjm/Desktop/work_documentation/JJ/chinese population analysis/manuscript/ajhg submission/review1/simulation/shanghai/502626-new-on-PC")
#setwd("C:/Documents and Settings/chenjm/Desktop/work_documentation/JJ/chinese population analysis/simulation analysis/simu2-chisq-removeSigni-removeConserved")


filename = "eigenstrat.poly.txt.chi2p"
#filename = "results-sim2-500N-cases-500S-controls.chisq.GC.txt.chi2p"
data <- read.table(filename, header=TRUE, sep="\t")

## pp-plot - NEED TO REMOVE MONOMORPHIC SNPs
## pvalues
size <- dim(data)[1]

expected.quantile <- (-log10(qunif(seq(from=1,to=size,by=1)/(size+1))))

## observed quantiles
obschisqP.quantile <- (-log10(sort(data$Chisq.p2))) ## sort works on 1 col
obschisqPgc.quantile <- (-log10(sort(data$Chisq.gc.p4))) ## sort works on 1 col
obschisqPeigen.quantile <- (-log10(sort(data$EIGENSTRAT.p3))) ## sort works on 1 col
obschisqPeigengc.quantile <- (-log10(sort(data$EIGENSTRAT.gc.p5))) ## sort works on 1 col

## for labelling
sortedP <- data[order(data$Chisq.p2, decreasing=FALSE),] ## order works on the whole matrix
sortedPgc <- data[order(data$Chisq.gc.p4, decreasing=FALSE),] ## order works on the whole matrix
sortedPeigen <- data[order(data$EIGENSTRAT.p3, decreasing=FALSE),] ## order works on the whole matrix
sortedPeigengc <- data[order(data$EIGENSTRAT.gc.p5, decreasing=FALSE),] ## order works on the whole matrix

##hist(observed.quantile)

x11()
plot(expected.quantile, obschisqP.quantile, xlab="-log(expected P values)", 
	ylab="-log(observed P values)", main=paste("QQ plot p-values: ",filename, sep=""))
#	ylab="-log(observed P values)", main="")

abline(a=0, b=1) # intercept,a = 0,0 b is slope x=y b=1
## identifying
#id <- which (obschisqP.quantile > 5)
#points(expected.quantile[id], obschisqP.quantile[id], col="red")
#text(expected.quantile[id], obschisqP.quantile[id], labels=sortedP$SNP_ID[id], pos=2, offset=0.5, col="red")

x11()
plot(expected.quantile, obschisqPgc.quantile, xlab="-log(expected P values)", 
	ylab="-log(observed P values)", main=paste("QQ plot p-values+gc: ",filename, sep=""))
#	ylab="-log(observed P values)", main="")

abline(a=0, b=1) # intercept,a = 0,0 b is slope x=y b=1
## identifying
#id <- which (obschisqPgc.quantile > 5)
#points(expected.quantile[id], obschisqPgc.quantile[id], col="red")
#text(expected.quantile[id], obschisqPgc.quantile[id], labels=sortedPgc$SNP_ID[id], pos=2, offset=0.5, col="red")

x11()
plot(expected.quantile, obschisqPeigen.quantile, xlab="-log(expected P values)", 
	ylab="-log(observed P values)", main=paste("QQ plot p-values+eigen: ",filename, sep=""))
#	ylab="-log(observed P values)", main="")

abline(a=0, b=1) # intercept,a = 0,0 b is slope x=y b=1
## identifying
#id <- which (obschisqPeigen.quantile > 5)
#points(expected.quantile[id], obschisqPeigen.quantile[id], col="red")
#text(expected.quantile[id], obschisqPeigen.quantile[id], labels=sortedPeigen$SNP_ID[id], pos=2, offset=0.5, col="red")

x11()
plot(expected.quantile, obschisqPeigengc.quantile, xlab="-log(expected P values)", 
	ylab="-log(observed P values)", main=paste("QQ plot p-values+eigen+gc: ",filename, sep=""))
#	ylab="-log(observed P values)", main="")

abline(a=0, b=1) # intercept,a = 0,0 b is slope x=y b=1
## identifying
#id <- which (obschisqPeigengc.quantile > 5)
#id <- which (sortedPeigengc$SNP_ID == "rs4980598" | sortedPeigengc$SNP_ID == "rs7641148" | sortedPeigengc$SNP_ID == "rs12141395")
#points(expected.quantile[id], obschisqPeigengc.quantile[id], col="red")
#text(expected.quantile[id], obschisqPeigengc.quantile[id], labels=sortedPeigengc$SNP_ID[id], pos=2, offset=0.5, col="red")

