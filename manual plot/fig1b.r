#setwd("C:/Documents and Settings/chenjm/Desktop/work_documentation/Paul/hhp-qc-new/1967samples-noSEA-panelA/2-east-asia")
setwd("E:/GIS@Work/work_documentation/Paul/w-hhp-qc-new/1967samples-noSEA-panelA/2-east-asia")

xComp="PC1"								#PC Component to plot on X Axis
yComp="PC2"								#PC Component to plot on Y Axis
baseFileName = "839samples-11125snps-es-hhp-combined-noSEA"		#Base file name of the .PCA and .EVAL file from fpca
annotFilename="839samples-hhp-pc1-lt0.015.sa" ##Annotation File

symbol = 3

pop6 = hesc
pop6n = "hESC"
pop6c = "red" ### note that pch = 21 in this one

pop1 = ami
pop1n = "Ami"
pop1c = "pink"

pop2 = atayal
pop2n = "Atayal"
pop2c = "red"

pop3 = cambodian
pop3n = "Cambodian"
pop3c = "lightblue"

pop4 = chb
pop4n = "CHB"
pop4c = "orange"

pop5 = hanchinese
pop5n = "Han Chinese"
pop5c = "tomato1"

pop7 = indian
pop7n = "Indian"
pop7c = "darkgreen"

pop8 = japanese
pop8n = "Japanese"
pop8c = "darkolivegreen2"

pop9 = jpt
pop9n = "JPT"
pop9c = "green"

pop10 = korean
pop10n = "Korean"
pop10c = "cyan"

pop11 = nonhanchinese
pop11n = "non-Han Chinese"
pop11c = "tan"

pop12 = ryukyuan
pop12n = "Ryukyuan"
pop12c = "brown"

pop13 = yakut
pop13n = "Yakut"
pop13c = "royalblue"

###############################################################################
###############################################################################
title = baseFileName	#Title of plot
pcaFilename = paste(baseFileName, ".pca", sep="")	#Change this if .pca file name diff from base file name
evalFilename= paste(baseFileName, ".eval", sep="")	#Change this if .eval file name diff from base file name
amtVariance <- read.table(evalFilename, sep="\t", header=T)
pca <- read.table(pcaFilename, header=T, sep="\t")
annot <- read.table(annotFilename, sep="\t", header=T, stringsAsFactors=F)

pca <- merge(annot, pca, by="sample.id")

pop1 <- which(pca$population.abbreviation == pop1n)
pop2 <- which(pca$population.abbreviation == pop2n)
pop3 <- which(pca$population.abbreviation == pop3n)
pop4 <- which(pca$population.abbreviation == pop4n)
pop5 <- which(pca$population.abbreviation == pop5n)
pop6 <- which(pca$population.abbreviation == pop6n)
pop7 <- which(pca$population.abbreviation == pop7n)
pop8 <- which(pca$population.abbreviation == pop8n)
pop9 <- which(pca$population.abbreviation == pop9n)
pop10 <- which(pca$population.abbreviation == pop10n)
pop11 <- which(pca$population.abbreviation == pop11n)
pop12 <- which(pca$population.abbreviation == pop12n)
pop13 <- which(pca$population.abbreviation == pop13n)

## max plotting region
par(pty = "m") 

xlabel = paste(xComp, " (", (amtVariance$percentage.of.variance[which(amtVariance$PC == xComp)]) * 100, " %)")
ylabel = paste(yComp, " (", (amtVariance$percentage.of.variance[which(amtVariance$PC == yComp)]) * 100, " %)")

#x11()
plot(-pca[[xComp]], pca[[yComp]], xlab = xlabel, ylab = ylabel, type = "n")

points(-pca[pop1,xComp], pca[pop1,yComp], pch = symbol, col = pop1c)
points(-pca[pop2,xComp], pca[pop2,yComp], pch = symbol, col = pop2c)
points(-pca[pop3,xComp], pca[pop3,yComp], pch = symbol, col = pop3c)
points(-pca[pop4,xComp], pca[pop4,yComp], pch = symbol, col = pop4c)
points(-pca[pop5,xComp], pca[pop5,yComp], pch = symbol, col = pop5c)
points(-pca[pop7,xComp], pca[pop7,yComp], pch = symbol, col = pop7c)
points(-pca[pop8,xComp], pca[pop8,yComp], pch = symbol, col = pop8c)
points(-pca[pop9,xComp], pca[pop9,yComp], pch = symbol, col = pop9c)
points(-pca[pop10,xComp], pca[pop10,yComp], pch = symbol, col = pop10c)
points(-pca[pop11,xComp], pca[pop11,yComp], pch = symbol, col = pop11c)
points(-pca[pop12,xComp], pca[pop12 ,yComp], pch = symbol, col = pop12c)
points(-pca[pop13,xComp], pca[pop13,yComp], pch = symbol, col = pop13c)
points(-pca[pop6,xComp], pca[pop6,yComp], pch = 21, col= pop6c, bg="black", cex=1.5) #################################################

#################################################

legend(locator(1), 
c("hESC","Ami","Atayal","Cambodian","CHB","Han Chinese",
"Indian","Japanese","JPT","Korean","non-Han Chinese","Ryukyuan","Yakut"), 
col = c("red","pink","red","lightblue","orange","tomato1","darkgreen",
"darkolivegreen2","green","cyan","tan","brown","royalblue"), 
text.col="black", 
pch = c(21,rep(3,times=12)),
pt.bg="black", bty = "o", cex=1.5, pt.cex = 2, pt.lwd=3, ncol=3)