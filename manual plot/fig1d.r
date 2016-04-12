#setwd("C:/Documents and Settings/chenjm/Desktop/work_documentation/Paul/w-hhp-jnov-QC-new")
setwd("E:/GIS@Work/work_documentation/Paul/w-jnov-QC-new")

xComp="PC2"								#PC Component to plot on X Axis
yComp="PC1"								#PC Component to plot on Y Axis
baseFileName = "1448samples-55972snps-merged-es-jnov"		#Base file name of the .PCA and .EVAL file from fpca
annotFilename="1448samples-merged-es-jnov.sa" ##Annotation File

flipX = 1
flipY = 0

symbol = 1
pop6 = "hesc"
pop6n = "hESC"
pop6c = "red" ### note that pch = 21 in this one

pop1 = "albania"
pop1n = "Albania"
pop1c = "lightblue"

pop2 = "austria"
pop2n = "Austria"
pop2c = "brown"

pop3 = "belgium"
pop3n = "Belgium"
pop3c = "red"

pop4 = "bosniaandherzegovina"
pop4n = "Bosnia-and-Herzegovina"
pop4c = "blue"

pop5 = "bulgaria"
pop5n = "Bulgaria"
pop5c = "blue"

pop7 = "croatia"
pop7n = "Croatia"
pop7c = "blue"

pop8 = "cyprus"
pop8n = "Cyprus"
pop8c = "red"

pop9 = "czech"
pop9n = "Czech-Republic"
pop9c = "green"

pop10 = "denmark"
pop10n = "Denmark"
pop10c = "lightblue"

pop11 = "finland"
pop11n = "Finland"
pop11c = "darkgreen"

pop12 = "france"
pop12n = "France"
pop12c = "orange"

pop13 = "germany"
pop13n = "Germany"
pop13c = "orange"

pop14 = "greece"
pop14n = "Greece"
pop14c = "blue"

pop15 = "hungary"
pop15n = "Hungary"
pop15c = "orange"

pop16 = "ireland"
pop16n = "Ireland"
pop16c = "red"

pop17 = "italy"
pop17n = "Italy"
pop17c = "brown"

pop18 = "kosovo"
pop18n = "Kosovo"
pop18c = "blue"

pop19 = "latvia"
pop19n = "Latvia"
pop19c = "green"

pop20 = "macedonia"
pop20n = "Macedonia"
pop20c = "lightblue"

pop21 = "netherlands"
pop21n = "Netherlands"
pop21c = "orange"

pop22 = "norway"
pop22n = "Norway"
pop22c = "lightblue"

pop23 = "poland"
pop23n = "Poland"
pop23c = "lightgreen"

pop24 = "portugal"
pop24n = "Portugal"
pop24c = "magenta"

pop25 = "romania"
pop25n = "Romania"
pop25c = "darkgreen"

pop26 = "russia"
pop26n = "Russian-Federation"
pop26c = "green"

pop27 = "scotland"
pop27n = "Scotland"
pop27c = "pink"

pop28 = "serbiaandmontenegro"
pop28n = "Serbia-and-Montenegro"
pop28c = "lightblue"

pop29 = "slovakia"
pop29n = "Slovakia"
pop29c = "darkgoldenrod1"

pop30 = "slovenia"
pop30n = "Slovenia"
pop30c = "blue"

pop31 = "spain"
pop31n = "Spain"
pop31c = "purple"

pop32 = "sweden"
pop32n = "Sweden"
pop32c = "lightseagreen"

pop33 = "switzerland"
pop33n = "Switzerland"
pop33c = "red"

pop34 = "turkey"
pop34n = "Turkey"
pop34c = "green"

pop35 = "ukraine"
pop35n = "Ukraine"
pop35c = "green"

pop36 = "uk"
pop36n = "United-Kingdom"
pop36c = "pink"

pop37 = "yugoslavia"
pop37n = "Yugoslavia"
pop37c = "blue"

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
pop14 <- which(pca$population.abbreviation == pop14n)
pop15 <- which(pca$population.abbreviation == pop15n)
pop16 <- which(pca$population.abbreviation == pop16n)
pop17 <- which(pca$population.abbreviation == pop17n)
pop18 <- which(pca$population.abbreviation == pop18n)
pop19 <- which(pca$population.abbreviation == pop19n)
pop20 <- which(pca$population.abbreviation == pop20n)
pop21 <- which(pca$population.abbreviation == pop21n)
pop22 <- which(pca$population.abbreviation == pop22n)
pop23 <- which(pca$population.abbreviation == pop23n)
pop24 <- which(pca$population.abbreviation == pop24n)
pop25 <- which(pca$population.abbreviation == pop25n)
pop26 <- which(pca$population.abbreviation == pop26n)
pop27 <- which(pca$population.abbreviation == pop27n)
pop28 <- which(pca$population.abbreviation == pop28n)
pop29 <- which(pca$population.abbreviation == pop29n)
pop30 <- which(pca$population.abbreviation == pop30n)
pop31 <- which(pca$population.abbreviation == pop31n)
pop32 <- which(pca$population.abbreviation == pop32n)
pop33 <- which(pca$population.abbreviation == pop33n)
pop34 <- which(pca$population.abbreviation == pop34n)
pop35 <- which(pca$population.abbreviation == pop35n)
pop36 <- which(pca$population.abbreviation == pop36n)
pop37 <- which(pca$population.abbreviation == pop37n)

## max plotting region
par(pty = "m") 

xlabel = paste(xComp, " (", (amtVariance$percentage.of.variance[which(amtVariance$PC == xComp)]) * 100, " %)")
ylabel = paste(yComp, " (", (amtVariance$percentage.of.variance[which(amtVariance$PC == yComp)]) * 100, " %)")

#laying the rules for graphwork: eg x11 is font
x11(16,9, pointsize = 18)
layout(matrix(c(rep(1,40), rep(2,10)), 1, 50, byrow = TRUE))
layout.show(2)

if(flipX == 1 & flipY == 0)
{
	plot(-pca[[xComp]], pca[[yComp]], xlab = xlabel, ylab = ylabel, type = "n")

	points(-pca[pop1,xComp], pca[pop1,yComp], pch = symbol, col = pop1c)
	points(-pca[pop2,xComp], pca[pop2,yComp], pch = symbol, col = pop2c)
	points(-pca[pop3,xComp], pca[pop3,yComp], pch = symbol, col = pop3c)
	points(-pca[pop4,xComp], pca[pop4,yComp], pch = symbol, col = pop4c)
	points(-pca[pop5,xComp], pca[pop5,yComp], pch = symbol, col = pop5c)
	points(-pca[pop7,xComp], pca[pop7,yComp], pch = symbol, col = pop7c)
	points(-pca[pop8,xComp],pca[pop8,yComp],pch=symbol,col=pop8c)
	points(-pca[pop9,xComp],pca[pop9,yComp],pch=symbol,col=pop9c)
	points(-pca[pop10,xComp],pca[pop10,yComp],pch=symbol,col=pop10c)
	points(-pca[pop11,xComp],pca[pop11,yComp],pch=symbol,col=pop11c)
	points(-pca[pop12,xComp],pca[pop12,yComp],pch=symbol,col=pop12c)
	points(-pca[pop13,xComp],pca[pop13,yComp],pch=symbol,col=pop13c)
	points(-pca[pop14,xComp],pca[pop14,yComp],pch=symbol,col=pop14c)
	points(-pca[pop15,xComp],pca[pop15,yComp],pch=symbol,col=pop15c)
	points(-pca[pop16,xComp],pca[pop16,yComp],pch=symbol,col=pop16c)
	points(-pca[pop17,xComp],pca[pop17,yComp],pch=symbol,col=pop17c)
	points(-pca[pop18,xComp],pca[pop18,yComp],pch=symbol,col=pop18c)
	points(-pca[pop19,xComp],pca[pop19,yComp],pch=symbol,col=pop19c)
	points(-pca[pop20,xComp],pca[pop20,yComp],pch=symbol,col=pop20c)
	points(-pca[pop21,xComp],pca[pop21,yComp],pch=symbol,col=pop21c)
	points(-pca[pop22,xComp],pca[pop22,yComp],pch=symbol,col=pop22c)
	points(-pca[pop23,xComp],pca[pop23,yComp],pch=symbol,col=pop23c)
	points(-pca[pop24,xComp],pca[pop24,yComp],pch=symbol,col=pop24c)
	points(-pca[pop25,xComp],pca[pop25,yComp],pch=symbol,col=pop25c)
	points(-pca[pop26,xComp],pca[pop26,yComp],pch=symbol,col=pop26c)
	points(-pca[pop27,xComp],pca[pop27,yComp],pch=symbol,col=pop27c)
	points(-pca[pop28,xComp],pca[pop28,yComp],pch=symbol,col=pop28c)
	points(-pca[pop29,xComp],pca[pop29,yComp],pch=symbol,col=pop29c)
	points(-pca[pop30,xComp],pca[pop30,yComp],pch=symbol,col=pop30c)
	points(-pca[pop31,xComp],pca[pop31,yComp],pch=symbol,col=pop31c)
	points(-pca[pop32,xComp],pca[pop32,yComp],pch=symbol,col=pop32c)
	points(-pca[pop33,xComp],pca[pop33,yComp],pch=symbol,col=pop33c)
	points(-pca[pop34,xComp],pca[pop34,yComp],pch=symbol,col=pop34c)
	points(-pca[pop35,xComp],pca[pop35,yComp],pch=symbol,col=pop35c)
	points(-pca[pop36,xComp],pca[pop36,yComp],pch=symbol,col=pop36c)
	points(-pca[pop37,xComp],pca[pop37,yComp],pch=symbol,col=pop37c)
	points(-pca[pop6,xComp], pca[pop6,yComp], pch = 21, col= pop6c, bg="black", cex=1.5) #################################################
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
	points(pca[pop8,xComp],pca[pop8,yComp],pch=symbol,col=pop8c)
	points(pca[pop9,xComp],pca[pop9,yComp],pch=symbol,col=pop9c)
	points(pca[pop10,xComp],pca[pop10,yComp],pch=symbol,col=pop10c)
	points(pca[pop11,xComp],pca[pop11,yComp],pch=symbol,col=pop11c)
	points(pca[pop12,xComp],pca[pop12,yComp],pch=symbol,col=pop12c)
	points(pca[pop13,xComp],pca[pop13,yComp],pch=symbol,col=pop13c)
	points(pca[pop14,xComp],pca[pop14,yComp],pch=symbol,col=pop14c)
	points(pca[pop15,xComp],pca[pop15,yComp],pch=symbol,col=pop15c)
	points(pca[pop16,xComp],pca[pop16,yComp],pch=symbol,col=pop16c)
	points(pca[pop17,xComp],pca[pop17,yComp],pch=symbol,col=pop17c)
	points(pca[pop18,xComp],pca[pop18,yComp],pch=symbol,col=pop18c)
	points(pca[pop19,xComp],pca[pop19,yComp],pch=symbol,col=pop19c)
	points(pca[pop20,xComp],pca[pop20,yComp],pch=symbol,col=pop20c)
	points(pca[pop21,xComp],pca[pop21,yComp],pch=symbol,col=pop21c)
	points(pca[pop22,xComp],pca[pop22,yComp],pch=symbol,col=pop22c)
	points(pca[pop23,xComp],pca[pop23,yComp],pch=symbol,col=pop23c)
	points(pca[pop24,xComp],pca[pop24,yComp],pch=symbol,col=pop24c)
	points(pca[pop25,xComp],pca[pop25,yComp],pch=symbol,col=pop25c)
	points(pca[pop26,xComp],pca[pop26,yComp],pch=symbol,col=pop26c)
	points(pca[pop27,xComp],pca[pop27,yComp],pch=symbol,col=pop27c)
	points(pca[pop28,xComp],pca[pop28,yComp],pch=symbol,col=pop28c)
	points(pca[pop29,xComp],pca[pop29,yComp],pch=symbol,col=pop29c)
	points(pca[pop30,xComp],pca[pop30,yComp],pch=symbol,col=pop30c)
	points(pca[pop31,xComp],pca[pop31,yComp],pch=symbol,col=pop31c)
	points(pca[pop32,xComp],pca[pop32,yComp],pch=symbol,col=pop32c)
	points(pca[pop33,xComp],pca[pop33,yComp],pch=symbol,col=pop33c)
	points(pca[pop34,xComp],pca[pop34,yComp],pch=symbol,col=pop34c)
	points(pca[pop35,xComp],pca[pop35,yComp],pch=symbol,col=pop35c)
	points(pca[pop36,xComp],pca[pop36,yComp],pch=symbol,col=pop36c)
	points(pca[pop37,xComp],pca[pop37,yComp],pch=symbol,col=pop37c)
	points(pca[pop6,xComp], pca[pop6,yComp], pch = 21, col= pop6c, bg="black") #################################################
}

#################################################
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(locator(1), 
c("hESC",
"South-eastern Europe",
"Albania",
"Bosnia-and-Herzegovina",
"Bulgaria",
"Croatia",
"Greece",
"Kosovo",
"Macedonia",
"Romania",
"Serbia-and-Montenegro",
"Slovenia",
"Turkey",
"Yugoslavia",
"South",
"Cyprus",
"Italy",
"South-western Europe",
"Portugal",
"Spain",
"North-eastern Europe",
"Czech-Republic",
"Finland",
"Hungary",
"Latvia",
"Poland",
"Russian-Federation",
"Slovakia",
"Ukraine",
"North",
"Norway",
"Denmark",
"Sweden",
"North-western Europe",
"Ireland",
"Netherlands",
"Scotland",
"United-Kingdom",
"Central",
"Austria",
"Belgium",
"France",
"Germany",
"Switzerland"
), 
col = c("red",
"white",
"lightblue",
"blue",
"blue",
"blue",
"blue",
"blue",
"lightblue",
"darkgreen",
"lightblue",
"blue",
"green",
"blue",
"white",
"red",
"brown",
"white",
"magenta",
"purple",
"white",
"green",
"darkgreen",
"orange",
"green",
"lightgreen",
"green",
"darkgoldenrod1",
"green",
"white",
"lightblue",
"lightblue",
"lightseagreen",
"white",
"red",
"orange",
"pink",
"pink",
"white",
"brown",
"red",
"orange",
"orange",
"red"
), 
text.col="black", 
pch = c(21,rep(1,times=43)),
pt.cex = 2, pt.bg="black", pt.lwd = 3, bty = "n", cex = 1.1)