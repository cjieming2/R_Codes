library(GDD)


#takes in sa and pca files only
#requires some changes with the colors for different sets of data

setwd("C:/Documents and Settings/chenjm/Desktop/work_documentation/JJ/surnames/2-pca/5-pca on 3971sa-373278snps no PC1-cor 0.17 snps")

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
		#plot(x, y, xlab=xlab,ylab=ylab, main=main, pch=20, bg="lightgrey", col=col[unclass(factor(categories))])
		
		## brings forward the points desired
		#id = which(data$population.abbreviation == "Sichuan")
		#id3 = which(data$population.abbreviation == "Shanghai")
		#id2 = which(data$population.abbreviation == "Guangdong-Cantonese")
		#id4 = which(data$population.abbreviation == "Guangdong-Hakka")
		#id5 = which(data$population.abbreviation == "Guangdong-Teochew")
		#points(x[id], y[id], pch=21, col="black", bg="peru")
		#points(x[id3], y[id3], pch=21, col="black", bg="orange1")
		#points(x[id2], y[id2], pch=21, col="black", bg="orange1")
		#points(x[id4], y[id4], pch=21, col="black", bg="yellow")
		#points(x[id5], y[id5], pch=21, col="black", bg="purple")
		id = which(data$surname == population[highlightIdx])
		points(x[id], y[id], pch=21, col="black", bg=colors[highlightIdx])


		
	}
	else
	{
		#plot(x, y, xlab=xlab,ylab=ylab, main=main, pch=21, col="black", bg=col[unclass(factor(categories))])
		plot(x, y, xlab=xlab,ylab=ylab, main=main, pch=19, col=col[unclass(factor(categories))])
		id = which(data$surname == population[highlightIdx])
		points(x[id], y[id], pch=21, col="black", bg=colors[highlightIdx])
		text(x[id], y[id], labels=data$province[id], pos=2, offset=0.5, col="red")
	}
	
	## add line
	#abline(v=0.1, col = "red")
	#abline(h=0.15, col = "red")
	#text(x=0.15,y=-0.1, labels="PC3=0.1", col="red")
	#text(x=-0.05,y=0.16, labels="PC1=0.15", col="red")
	
	# highlight points in red
	#id <- which (data$sample.id == "AD788" | data$sample.id == "CH-SC6" | data$sample.id == "AD927" | data$sample.id == "AD1140" | data$sample.id == "AD898" | data$sample.id == "AD829")
	#id2 <- which (data$sample.id == "CH-SC6")
	#id3 <- which (x < -0.4)
	#id4 <- which (y > 0.28 | y < -0.2)
	#id5 <- which (data$sample.id == "AD374" | data$sample.id == "AD1157")
	#points(x[id], y[id], col="red")
	#points(x[id2], y[id2], col="purple")
	#points(x[id3], y[id3], col="red")
	#points(x[id4], y[id4], col="red")
	#points(x[id5], y[id5], col="red")

	# add label sample-id on the points
	# pos 1-bottom 2-left 3-top 4-right
	# offset distance away from point
	#text(x[id], y[id], labels=data$sample.id[id], pos=2, offset=0.5, col="red")	
	#text(x[id2], y[id2], labels=data$sample.id[id2], adj=c(1,1.1), offset=0.5, col="purple")	
	#text(x[id3], y[id3], labels=data$sample.id[id3], pos=4, adj=c(1,2), offset=0.5, col="black")	
	#text(x[id4], y[id4], labels=data$sample.id[id4], pos=2, adj=c(0,1.5), offset=0.5, col="black")	
	#text(x[id5], y[id5], labels=data$sample.id[id5], pos=2, adj=c(0,1.5), offset=0.5, col="black")

	par(mar = c(0,0,0,0))
	plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
	legend(-1, 0.8, sort(unique(categories)), pt.bg = col,
	#legend(-1, 0.8, poplegend, pt.bg = colorslegend,
      	 text.col = "black", pch = 21,
        	 bg = 'white')
}

#Edit the below parameters
#xComp and yComp with the respective component that needs to be plotted
#baseFileName with the fileName of the .pca and .eval file. If filename is not the same
#as the plot or of non-standard naming, change the respective pcaFilename and/or 
#evalFilename variable below.

##----- Variables to change ----------------
xComp="PC2"								#PC Component to plot on X Axis
yComp="PC1"								#PC Component to plot on Y Axis
generateHighlight = 1						#1 to generate highlight of each population, 0 to skip
baseFileName = "db126-3970samples-371946snps-surnames-removedforpca-lt-0.17wPC1"		#Base file name of the .PCA and .EVAL file from fpca
title = baseFileName	#Title of plot

annotFilename="split-11-95surnames.list.sa" ##Annotation File

colors <- c("green", "red", "blue", "orange", "yellow", "pink", "purple","brown","cornflowerblue","slategray1","violetred","tan","cyan","deeppink","black","lightblue3","goldenrod","magenta","khaki","indianred3","lightgreen") ##Set of colours to use for all plots

#colorslegend <- c("lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey","lightgrey", "yellow", "black") ## multicoloured but name in order of N-S + Sichuan shanghai + CHB Spore GREYed

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
population <- sort(unique(data$surname))
#poplegend <- c("Liaoning", "Henan", "Hebei", "Shandong", "Anhui", "Hubei", "Jiangsu", "Hunan", "Guangdong") ## no shanghai, sichuan
#poplegend <- c("Liaoning", "Henan", "Hebei", "Shandong", "Anhui", "Hubei", "Jiangsu", "Shanghai", "Sichuan", "Hunan", "Guangdong-Cantonese", "Guangdong-Hakka", "Guangdong-Teochew","CHB", "Singapore")

numOfPop = length(population)

#Generate X and Y Axis label with amt of variance
xLabel = paste(xComp, " (", (amtVariance$percentage.of.variance[which(amtVariance$PC == xComp)]) * 100, " %)")
yLabel = paste(yComp, " (", (amtVariance$percentage.of.variance[which(amtVariance$PC == yComp)]) * 100, " %)")


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

		pca.plot(data[[xComp]], data[[yComp]], xlab=xLabel, ylab=yLabel , main=paste(title, " - ", population[i]), fgColors, data$surname, i)
		#pca.plot(data[[xComp]], data[[yComp]], xlab=xLabel, ylab=yLabel , main="", fgColors, data$province, i)
	}
}


pca.plot(data[[xComp]], data[[yComp]], xlab=xLabel, ylab=yLabel , main=title, colors, data$surname)
#pca.plot(data[[xComp]], data[[yComp]], xlab=xLabel, ylab=yLabel, main='' , colors, data$province)

#x11()
#pairs(data[,4:9], main="", pch = 20, col = colors[unclass(factor(data$surname))])