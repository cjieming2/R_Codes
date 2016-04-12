#setwd("C:/Documents and Settings/chenjm/Desktop/work_documentation/Paul/hhp-qc-new/1967samples-noSEA-panelA/5-the rest")
setwd("E:/GIS@Work/work_documentation/Paul/w-hhp-qc-new/1967samples-noSEA-panelA/5-the rest")

xComp="PC2"								#PC Component to plot on X Axis
yComp="PC1"								#PC Component to plot on Y Axis
baseFileName = "715samples-11273snps-es-hhp-combined-noSEA"		#Base file name of the .PCA and .EVAL file from fpca
annotFilename="715samples-the-rest-es-hhp.sa" ##Annotation File

flipX = 1
flipY = 1

symbol = 3

pop6 = "hesc"
pop6n = "hESC"
pop6c = "red" ### note that pch = 21 in this one

pop1 = "ceu"
pop1n = "CEU"
pop1c = "cyan"

pop2 = "hgdpcentralsouthasia"
pop2n = "HGDP-CentralSouth-Asia"
pop2c = "green"

pop3 = "hgdpeurope"
pop3n = "HGDP-Europe"
pop3c = "blue"

pop4 = "hgdpmiddleeast"
pop4n = "HGDP-Middle-East"
pop4c = "brown"

pop5 = "pasnpcentralsouthasia"
pop5n = "PASNP-CentralSouth-Asia"
pop5c = "darkgreen"

pop7 = "pasnpeastasia"
pop7n = "PASNP-East-Asia"
pop7c = "tan"

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

## max plotting region
par(pty = "m") 

xlabel = paste(xComp, " (", (amtVariance$percentage.of.variance[which(amtVariance$PC == xComp)]) * 100, " %)")
ylabel = paste(yComp, " (", (amtVariance$percentage.of.variance[which(amtVariance$PC == yComp)]) * 100, " %)")

#x11()

if(flipX == 1 & flipY == 0)
{
	plot(-pca[[xComp]], pca[[yComp]], xlab = xlabel, ylab = ylabel, type = "n")

	points(-pca[pop1,xComp], pca[pop1,yComp], pch = symbol, col = pop1c)
	points(-pca[pop2,xComp], pca[pop2,yComp], pch = symbol, col = pop2c)
	points(-pca[pop3,xComp], pca[pop3,yComp], pch = symbol, col = pop3c)
	points(-pca[pop4,xComp], pca[pop4,yComp], pch = symbol, col = pop4c)
	points(-pca[pop5,xComp], pca[pop5,yComp], pch = symbol, col = pop5c)
	points(-pca[pop7,xComp], pca[pop7,yComp], pch = symbol, col = pop7c)
	points(-pca[pop6,xComp], pca[pop6,yComp], pch = 21, col= pop6c, bg="black",cex=1.5) #################################################
} else
if(flipY == 1 & flipX == 0)
{
	plot(pca[[xComp]], -pca[[yComp]], xlab = xlabel, ylab = ylabel, type = "n")

	points(pca[pop1,xComp], -pca[pop1,yComp], pch = symbol, col = pop1c)
	points(pca[pop2,xComp], -pca[pop2,yComp], pch = symbol, col = pop2c)
	points(pca[pop3,xComp], -pca[pop3,yComp], pch = symbol, col = pop3c)
	points(pca[pop4,xComp], -pca[pop4,yComp], pch = symbol, col = pop4c)
	points(pca[pop5,xComp], -pca[pop5,yComp], pch = symbol, col = pop5c)
	points(pca[pop7,xComp], -pca[pop7,yComp], pch = symbol, col = pop7c)
	points(pca[pop6,xComp], -pca[pop6,yComp], pch = 21, col= pop6c, bg="black",cex=1.5) #################################################
} else
if((flipX == 1) & (flipY == 1))
{
	plot(-pca[[xComp]], -pca[[yComp]], xlab = xlabel, ylab = ylabel, type = "n")

	points(-pca[pop1,xComp], -pca[pop1,yComp], pch = symbol, col = pop1c)
	points(-pca[pop2,xComp], -pca[pop2,yComp], pch = symbol, col = pop2c)
	points(-pca[pop3,xComp], -pca[pop3,yComp], pch = symbol, col = pop3c)
	points(-pca[pop4,xComp], -pca[pop4,yComp], pch = symbol, col = pop4c)
	points(-pca[pop5,xComp], -pca[pop5,yComp], pch = symbol, col = pop5c)
	points(-pca[pop7,xComp], -pca[pop7,yComp], pch = symbol, col = pop7c)
	points(-pca[pop6,xComp], -pca[pop6,yComp], pch = 21, col= pop6c, bg="black",cex=1.5) #################################################

} else
if((flipX == 0) & (flipY == 0))
{
	plot(pca[[xComp]], pca[[yComp]], xlab = xlabel, ylab = ylabel, type = "n")

	points(pca[pop1,xComp], pca[pop1,yComp], pch = symbol, col = pop1c)
	points(pca[pop2,xComp], pca[pop2,yComp], pch = symbol, col = pop2c)
	points(pca[pop3,xComp], pca[pop3,yComp], pch = symbol, col = pop3c)
	points(pca[pop4,xComp], pca[pop4,yComp], pch = symbol, col = pop4c)
	points(pca[pop5,xComp], pca[pop5,yComp], pch = symbol, col = pop5c)
	points(pca[pop7,xComp], pca[pop7,yComp], pch = symbol, col = pop7c)
	points(pca[pop6,xComp], pca[pop6,yComp], pch = 21, col= pop6c, bg="black",cex=1.5) #################################################
}

#################################################

legend(locator(1), 
c("hESC","CEU","HGDP-CentralSouth-Asia","HGDP-Europe","HGDP-Middle-East",
"PASNP-CentralSouth-Asia","PASNP-East-Asia"), 
col = c("red","cyan","green","blue","brown","darkgreen","tan"), 
text.col="black", 
pch = c(21,rep(3,times=6)),
pt.cex = 2, pt.bg="black", pt.lwd = 3, bty = "o", cex = 1.5)