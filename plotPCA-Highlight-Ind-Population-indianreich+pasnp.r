#takes in sa and pca files only
#requires some changes with the colors for different sets of data

setwd("C:/Documents and Settings/chenjm/Desktop/work_documentation/Mark/ihhp/pca")
#setwd("C:/Documents and Settings/chenjm/Desktop/work_documentation/JJ/chinese population analysis/manuscript/ajhg submission/review1/pca/east-west")

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
		id = which(data$population.abbreviation == "CEU")
		id3 = which(data$population.abbreviation == "CHB")
		id2 = which(data$population.abbreviation == "JPT")
		id4 = which(data$population.abbreviation == "YRI")
		#id5 = which(data$population.abbreviation == "Mala")
		#id6 = which(data$population.abbreviation == "Naidu")
		#id7 = which(data$population.abbreviation == "Velama")
		#id8 = which(data$population.abbreviation == "Vysya")


		points(x[id], y[id], pch=21, col="black", bg="cyan")
		points(x[id2], y[id2], pch=21, col="black", bg="pink")
		points(x[id3], y[id3], pch=21, col="black", bg="red")
		points(x[id4], y[id4], pch=21, col="black", bg="brown")
		#points(x[id5], y[id5], pch=21, col="black", bg="magenta")
		#points(x[id6], y[id6], pch=21, col="black", bg="cyan")
		#points(x[id7], y[id7], pch=21, col="black", bg="purple")
		#points(x[id8], y[id8], pch=21, col="black", bg="pink")




		
	}
	else
	{
		plot(x, y, xlab=xlab,ylab=ylab, main=main, pch=19, col=col[unclass(factor(categories))])
		#id = which(data$population.abbreviation == population[highlightIdx])

		#id1 = which(data$dataset  == "CEU")
		#id2 = which(data$dataset  == "CHB")
		#id3 = which(data$dataset  == "JPT")
		#id4 = which(data$dataset  == "YRI")

		#points(x[id1], y[id1], pch=21, col="black", bg="magenta")
		#points(x[id2], y[id2], pch=21, col="black", bg="orange")
		#points(x[id3], y[id3], pch=21, col="black", bg="green")
		#points(x[id4], y[id4], pch=21, col="black", bg="royalblue")

		id = which(data$dataset == population[highlightIdx])
		points(x[id], y[id], pch=21, col="black", bg=colors[highlightIdx])

		#text(x=-0.045,y=-0.04, labels="CEU", col="black")
		#text(x=0.02,y=0.01, labels="CHB", col="black")
		#text(x=0.02,y=0.005, labels="JPT", col="black")
		#text(x=0.06,y=-0.1, labels="YRI", col="black")

		
	}
	
	## add line
	#abline(v=0.1, col = "red")
	#abline(h=-0.005, col = "red")
	#text(x=0.15,y=-0.1, labels="PC3=0.1", col="red")
	#text(x=-0.6,y=0, labels="PC1=-0.005", col="red")
	
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
generateHighlight = 0						#1 to generate highlight of each population, 0 to skip
baseFileName = "2932samples-5557snps-ihhp-db126-fwd"		#Base file name of the .PCA and .EVAL file from fpca
title = baseFileName	#Title of plot
annotFilename="2932samples-ihhp.sa" ##Annotation File
colors <- c("cyan", 
"red", 
"magenta", 
"orange", 
"yellow", 
"pink", 
"purple",
"brown",
"black",
"slategray1",
"violetred",
"tan",
"darksalmon", 
"deeppink",
"darkgreen", 
"orchid",
"green",
"antiquewhite3", 
"blue",
"darkblue", 
"peru",
"slateblue",
"thistle",
"tomato",
"rosybrown1",
"olivedrab",
"royalblue",
"coral3",
"bisque3") ##Set of 29 colours to use for all plots

#colors <- c("green","red","brown","blue","magenta","cyan","lightgrey","purple","pink","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey")

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
#population <- sort(unique(data$population.abbreviation))
population <- sort(unique(data$dataset))
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

		#pca.plot(data[[xComp]], data[[yComp]], xlab=xLabel, ylab=yLabel , main=paste(title, " - ", population[i]), fgColors, data$dataset, i)
		#pca.plot(data[[xComp]], data[[yComp]], xlab=xLabel, ylab=yLabel , main="", fgColors, data$population.abbreviation, i)
		pca.plot(data[[xComp]], data[[yComp]], xlab=xLabel, ylab=yLabel , main="", fgColors, data$dataset, i)
	}
}


pca.plot(data[[xComp]], data[[yComp]], xlab=xLabel, ylab=yLabel , main=title, colors, data$dataset)
#pca.plot(data[[xComp]], data[[yComp]], xlab=xLabel, ylab=yLabel, main='' , colors, data$population.abbreviation)

x11()
pairs(data[,6:9], main="", pch = 20, col = colors[unclass(factor(data$dataset))])
#pairs(data[,11:14], main="", pch = 20, col = colors[unclass(factor(data$dataset))])

#pairs(data[,4:8], main="", pch = 20, col = colors[unclass(factor(data$population.abbreviation))])