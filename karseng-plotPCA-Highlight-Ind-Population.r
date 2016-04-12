#takes in sa and pca files only
#requires some changes with the colors for different sets of data

setwd("D:/DataAnalysisWorkingPath/")
#setwd("D:/DataAnalysisWorkingPath/Chinese Population PCA/")

#Edit the below parameters
#xComp and yComp with the respective component that needs to be plotted
#baseFileName with the fileName of the .pca and .eval file. If filename is not the same
#as the plot or of non-standard naming, change the respective pcaFilename and/or 
#evalFilename variable below.

##----- Set up plots here ----------------
xComp="PC1"								#PC Component to plot on X Axis
yComp="PC2"								#PC Component to plot on Y Axis
generateHighlight = 0						#1 too generate highlight of each population, 0 to skip
#baseFileName = "5479samples-519883snps-58c-nbs-ukmgc-combined-top"		#Base file name of the .PCA and .EVAL file from fpca
#baseFileName = "4704samples-52033snps-wtccc2-58c-nbs-top-snpset2"
baseFileName = "4155samples-51212snps-58c-nbs-affy6-snpset1"


#title = "1338 samples 4190235 snps (conserved region removed)"	#Title of plot
title = baseFileName
#title = "Test"
#annotFilename="1382samples-with-dialect.sa"
annotFilename="4155samples-58c-nbs-affy6.sa"

#colors <- c("orange", "green", "red", "blue", "black", "pink", "tan","brown","purple","slategray1","yellow","cyan","violetred","deeppink")	#Set of colours to use for all plots
#colors <- c("orange", "yellow", "limegreen", "red", "steelblue", "magenta", "brown", "pink", "tan","black","purple","lightskyblue","lightgreen","cyan","violetred","deeppink")
#colors=c("red2","green","orange1","royalblue","yellow3","darkslategrey","purple3","turquoise","hotpink","lightgreen","salmon","skyblue1","goldenrod","slategrey","purple1","maroon","darkgreen","orange4","darkblue","brown","gray","mediumpurple")
colors=c("orange1","green","red2","royalblue","yellow3","darkslategrey","purple3","turquoise","hotpink","lightgreen","salmon","skyblue1","goldenrod","slategrey","purple1","maroon","darkgreen","orange4","darkblue","brown","gray","mediumpurple")


#colors <- c("orange", "green", "red", "blue", "black", "pink", "tan","brown","purple","slategray1","yellow", "darkseagreen","steelblue3", "cyan")	#Set of colours used for chinese population studies (including sha/sze)
#colors <- c("blue", "darkgoldenrod1", "forestgreen", "deeppink3", "magenta3", "deepskyblue3", "lawngreen","dodgerblue4","indianred2","darkorchid","cornflowerblue", "seashell4","orangered1")	#Set of colours used for chinese population studies (including sha/sze) - N-S colouring

#colors <- c("orange", "steelblue", "magenta", "red")
#colors <- c("magenta", "cyan", "yellow", "red", "steelblue", "magenta", "brown", "pink", "tan","black","purple","lightskyblue","lightgreen","cyan","violetred","deeppink")


pcaFilename = paste(baseFileName, ".pca", sep="")	#Change this if .pca file name diff from base file name
evalFilename= paste(baseFileName, ".eval", sep="")	#Change this if .eval file name diff from base file name
##----- End set up ----------------

#DO NOT ADD ETHNICITY in SA FILE!!! cos it's got spaces, will screw up the read.table function
#annotation loading sa file
annot <- read.table(annotFilename, header=TRUE, sep="\t")

#pca data loading pca file
data <- read.table(pcaFilename, header=TRUE, sep="\t")

#merge annot and pca data together
data <- merge(annot, data, by='sample.id')

#read eval file
#evalFilename="db126-1293samples-45CHBsamples-7454snps-conserved-region.eval"
amtVariance <- read.table(evalFilename, header=TRUE, sep="\t")

#Get a sorted unique set of populations and number of populations
population <- sort(unique(data$population.abbreviation))
#population <- sort(unique(data$region))
#population <- sort(unique(data$sex))


numOfPop = length(population)

#Generate X and Y Axis label with amt of variance
xLabel = paste(xComp, " (", (amtVariance$percentage.of.variance[which(amtVariance$PC == xComp)] * 100), " %)")
yLabel = paste(yComp, " (", (amtVariance$percentage.of.variance[which(amtVariance$PC == yComp)] * 100), " %)")

if (generateHighlight == 1)
{
	#Go through each population and plot PCA with each population highlighted and remaining populations greyed
	for (i in numOfPop:1)
	{
		#Init colour set for highlighting
		fgColors="grey85"
		if (numOfPop > 1)
			for (j in 1:numOfPop)
				fgColors[j] = "grey85"
		fgColors[i] = colors[i]
		pca.plot(data[[xComp]], data[[yComp]], xlab=xLabel, ylab=yLabel , main=paste(title, " - ", population[i]), fgColors, data$population.abbreviation, i)
#		pca.plot(data[[xComp]], data[[yComp]], xlab=xLabel, ylab=yLabel , main=paste(title, " - ", population[i]), fgColors, data$region, i)
#		pca.plot(data[[xComp]], data[[yComp]], xlab=xLabel, ylab=yLabel , main=paste(title, " - ", population[i]), fgColors, data$sex, i)
	}
}

pca.plot(data[[xComp]], data[[yComp]], xlab=xLabel, ylab=yLabel , main=title, colors, data$population.abbreviation)
#pca.plot(data[[xComp]], data[[yComp]], xlab=xLabel, ylab=yLabel , main=title, colors, data$region)
#pca.plot(data[[xComp]], data[[yComp]], xlab=xLabel, ylab=yLabel , main=title, colors, data$sex)




#pca.plot function
pca.plot <- function (x, y, xlab, ylab, main, col, categories, highlightIdx=0, lweight=40, rweight=10, ...)
{
	#laying the rules for graphwork: eg x11 is font
	x11(16,9, pointsize = 18)
	layout(matrix(c(rep(1,lweight), rep(2,rweight)), 1, 50, byrow = TRUE))
	layout.show(2)

	if(highlightIdx == 0)
	{
		plot(x, y, xlab=xlab,ylab=ylab, main=main, pch=21, col="black", bg=col[unclass(factor(categories))])
#		abline(v=0.0025, col="red")

#		abline(v = -0.0004, col="red")
#		text(-0.0004,-0.09,labels="PC1 = -0.0004",pos=4,offset=0.5, col="red")

#		abline(h = -0.025, col="red")
#		text(0.22,-0.025,labels="PC2 = -0.025",pos=1,offset=0.5, col="red")


		#which(data$population.abbreviation == "58c" & data$PC1 > -0.0011)
		
#		id = which(data$population.abbreviation == "Shanghai" & (data$PC2 < -0.05 | data$PC3 > 0.05))
#		id = which(data$PC5 > 0.2)
#		text(data$PC4[id],data$PC5[id],labels=data$sample.id[id],pos=2,offset=0.5, col="red")

#		id = which(data$PC4 < -0.3)
#		text(data$PC4[id],data$PC5[id],labels=data$sample.id[id],pos=4,offset=0.5, col="red")

	}
	else
	{
		plot(x, y, xlab=xlab,ylab=ylab, main=main, pch=19, col=col[unclass(factor(categories))])
		id = which(data$population.abbreviation == population[highlightIdx])
#		id = which(data$region == population[highlightIdx])
#		id = which(data$sex == population[highlightIdx])

		points(x[id], y[id], pch=21, col="black", bg=colors[highlightIdx])

#		abline(v= -0.00, col="red")
#		text(-0.02,-0.03,labels="PC2 = 0.05",pos=4,offset=0.5, col="red")
#		abline(v = -0.0031, col="red")
#		text(-0.0031,-0.09,labels="PC1 = -0.0031",pos=4,offset=0.5, col="red")
	}

	par(mar = c(0,0,0,0))
	plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
	legend(-1, 0.8, sort(unique(categories)), pt.bg = col,
      	 text.col = "black", pch = 21,
        	 bg = 'white')
}
