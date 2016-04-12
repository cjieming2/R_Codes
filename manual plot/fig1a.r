#setwd("C:/Documents and Settings/chenjm/Desktop/work_documentation/Paul/hhp-qc-new/1967samples-noSEA-panelA")
setwd("E:/GIS@Work/work_documentation/Paul/w-hhp-qc-new/1967samples-noSEA-panelA")

xComp="PC2"								#PC Component to plot on X Axis
yComp="PC1"								#PC Component to plot on Y Axis
baseFileName = "1967samples-11279snps-es-hhp-combined-noSEA"		#Base file name of the .PCA and .EVAL file from fpca
title = baseFileName	#Title of plot
annotFilename="1967samples-hhp-es-noSEA.sa" ##Annotation File

pcaFilename = paste(baseFileName, ".pca", sep="")	#Change this if .pca file name diff from base file name
evalFilename= paste(baseFileName, ".eval", sep="")	#Change this if .eval file name diff from base file name
amtVariance <- read.table(evalFilename, sep="\t", header=T)
pca <- read.table(pcaFilename, header=T, sep="\t")
annot <- read.table(annotFilename, sep="\t", header=T, stringsAsFactors=F)

pca <- merge(annot, pca, by="sample.id")

hesc <- which(pca$population.abbreviation == "hESC")
ceu <- which(pca$population.abbreviation == "CEU")
chb <- which(pca$population.abbreviation == "CHB")
jpt <- which(pca$population.abbreviation == "JPT")
yri <- which(pca$population.abbreviation == "YRI")
hgdpafrica <- which(pca$population.abbreviation == "HGDP-Africa")
hgdpamerica <- which(pca$population.abbreviation == "HGDP-America")
hgdpeurope <- which(pca$population.abbreviation == "HGDP-Europe")
hgdpmiddleeast <- which(pca$population.abbreviation == "HGDP-Middle-East")
hgdpoceania <- which(pca$population.abbreviation == "HGDP-Oceania")
hgdpcentralsouthasia <- which(pca$population.abbreviation == "HGDP-CentralSouth-Asia")
hgdpeastasia <- which(pca$population.abbreviation == "HGDP-East-Asia")
pasnpcentralsouthasia <- which(pca$population.abbreviation == "PASNP-CentralSouth-Asia")
pasnpeastasia <- which(pca$population.abbreviation == "PASNP-East-Asia")

## max plotting region
par(pty = "m") 

xlabel = paste(xComp, " (", (amtVariance$percentage.of.variance[which(amtVariance$PC == xComp)]) * 100, " %)")
ylabel = paste(yComp, " (", (amtVariance$percentage.of.variance[which(amtVariance$PC == yComp)]) * 100, " %)")

#x11()
#plot(pca[[xComp]], pca[[yComp]], xlab = xlabel, ylab = ylabel, type = "n")
plot(pca[[xComp]], pca[[yComp]], xaxt='n', yaxt='n', axes=F, ann=F, type = "n")

symbol = 3
points(pca[ceu,xComp], pca[ceu,yComp], pch = symbol, col = "cyan")
points(pca[chb,xComp], pca[chb,yComp], pch = symbol, col = "orange")
points(pca[hgdpafrica,xComp], pca[hgdpafrica,yComp], pch = symbol, col = "pink")
points(pca[yri,xComp], pca[yri,yComp], pch = symbol, col = "red")
points(pca[hgdpamerica,xComp], pca[hgdpamerica,yComp], pch = symbol, col = "purple")
points(pca[hgdpeurope,xComp], pca[hgdpeurope,yComp], pch = symbol, col = "blue")
points(pca[hgdpmiddleeast,xComp], pca[hgdpmiddleeast,yComp], pch = symbol, col = "brown")
points(pca[hgdpoceania,xComp], pca[hgdpoceania,yComp], pch = symbol, col = "magenta")
points(pca[hgdpcentralsouthasia,xComp], pca[hgdpcentralsouthasia,yComp], pch = symbol, col = "green")
points(pca[hgdpeastasia,xComp], pca[hgdpeastasia,yComp], pch = symbol, col = "darkolivegreen2")
points(pca[pasnpcentralsouthasia,xComp], pca[pasnpcentralsouthasia ,yComp], pch = symbol, col = "darkgreen")
points(pca[pasnpeastasia,xComp], pca[pasnpeastasia,yComp], pch = symbol, col = "tan")
points(pca[jpt,xComp], pca[jpt,yComp], pch = symbol, col = "yellow")
points(pca[hesc,xComp], pca[hesc,yComp], pch = 21, col= "red", bg="black", cex=1.5)

legend(locator(1), c("hESC", "CEU", "CHB", "JPT", "YRI", "HGDP-Africa", "HGDP-America", 
"HGDP-Europe", "HGDP-Middle-East", "HGDP-Oceania", "HGDP-CentralSouth-Asia", 
"HGDP-East-Asia", "PASNP-CentralSouth-Asia", "PASNP-East-Asia"), 
col = c("red","cyan","orange","yellow","red","pink","purple","blue","brown","magenta",
"green","darkolivegreen2","darkgreen","tan"), text.col="black", pch = c(21,rep(3,times=13)),
pt.bg="black", bty = "o", cex=1.5, pt.cex = 2, pt.lwd=3)