#takes in sa and pca files only
#requires some changes with the colors for different sets of data

#setwd("C:/Documents and Settings/chenjm/Desktop/work_documentation/JJ/chinese population analysis/manuscript/ajhg submission/review1/pca/1109china-metro-spore-dia/tagsnps")
#setwd("C:/Documents and Settings/chenjm/Desktop/work_documentation/JJ/chinese population analysis/manuscript/ajhg submission/review1/simulation/shanghai")
setwd("C:/Users/JM/thesis/mark_work/23andme/pca")

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
		#id = which(data$population.abbreviation == "CHB")
		#id = which(data$PC2 < -0.2)

		#id2 = which(data$population.abbreviation == "Shanghai")
		#id3 = which(data$population.abbreviation == "Singapore")

		#id4 = which(data$population.abbreviation == "Guangdong-Cantonese")
		#id5 = which(data$population.abbreviation == "Guangdong-Hakka")
		#id6 = which(data$population.abbreviation == "Guangdong-Teochew")
		#points(x[id], y[id], pch=21, col="black", bg="purple")
		#points(x[id3], y[id3], pch=21, col="black", bg="yellow")
		#points(x[id2], y[id2], pch=21, col="black", bg="orange1")
		#points(x[id4], y[id4], pch=21, col="black", bg="orange1")
		#points(x[id5], y[id5], pch=21, col="black", bg="yellow")
		#points(x[id6], y[id6], pch=21, col="black", bg="purple")
    id5 = which(data$population.abbreviation == "23andme_jc")

		#points(x[id2], y[id2], pch=21, col="black", bg="lightgrey")
    points(x[id5], y[id5], pch=21, col="red", bg="black", cex=2)

		#text(x[id], y[id], data$sample.id[id], pos=4, adj=c(0,1.5), offset=0.5, col="black")
    text(x[id5], y[id5], labels="23andme_jc", pos=2, adj=c(0,1.5), offset=0.5, col="black")
		
	}
	else
	{
		plot(x, y, xlab=xlab,ylab=ylab, main=main, pch=19, col=col[unclass(factor(categories))])
		id = which(data$population.abbreviation == population[highlightIdx])
		points(x[id], y[id], pch=21, col="black", bg=colors[highlightIdx])
		
		## brings forward the points desired
		#id = which(data$population.abbreviation == "Anhui")
		#id2 = which(data$population.abbreviation == "Guangdong")
		#id3 = which(data$population.abbreviation == "CHB")
		#id4 = which(data$population.abbreviation == "Singapore")
		#points(x[id], y[id], pch=21, col="black", bg="blue")
		#points(x[id2], y[id2], pch=21, col="black", bg="green")
		#points(x[id3], y[id3], pch=21, col="black", bg="yellow")
		#points(x[id4], y[id4], pch=21, col="black", bg="black")
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
	#legend(-1, 0.8, sort(unique(categories)), pt.bg = col,
	legend(-1, 0.8, poplegend, pt.bg = colorslegend,
      	 text.col = "black", pch = 21, bg = 'white')
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
#baseFileName = "1109samples-107565snps-china-metro-hap610quad-top"		#Base file name of the .PCA and .EVAL file from fpca
baseFileName = "sieved-sieved-1499samples-23625snps-ss+chinese+jpt+jc_23andme"
title = baseFileName	#Title of plot
#annotFilename="1109sample-china-metro.sa" ##Annotation File
annotFilename="sieved-1499samples-ss+chinese+jpt+jc_23andme.sa"
colors <- c("black","royalblue","greenyellow","magenta","purple","cyan","darkgreen","yellow1","red","lightpink","brown")
#colors <- c("black", "red", "blue", "orange", "yellow", "pink", "purple","brown","lightblue","slategray1","violetred","tan","cyan","deeppink","green") ##Set of colours to use for all plots
#colors <- c("blue", "green", "cornflowerblue", "cadetblue", "cyan", "seagreen", "lightblue", "deeppink","purple") ## multicoloured but name in order of N-S
#colors <- c("lightgrey", "purple", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "orange1", "lightgrey", "yellow") ## metro only
#colors <- c("royalblue", "greenyellow", "magenta", "purple", "cyan", "darkgreen", "yellow", "red", "lightpink", "brown") ## china-metro-spore
#colors <- c("lightgrey", "orange1", "yellow", "purple", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey","lightgrey", "lightgrey") ## dialect only
#colors <- c("royalblue", "greenyellow", "magenta", "purple", "cyan", "darkgreen", "yellow", "red", "lightpink", "lightgrey", "brown") ## china-metro-spore+1000shanghai


#colorslegend <- c("deeppink", "cadetblue", "cornflowerblue", "purple", "blue", "cyan", "lightblue", "seagreen","green") ## china-shanghai-sichuan
#colorslegend <- c("orangered1", "palevioletred3", "indianred1", "orange1", "orange1", "yellow1", "yellow1", "yellow1","yellow1", "greenyellow", "chartreuse3", "deepskyblue3", "lightblue2") ## china+metro+spore
colorslegend <- c("red", "magenta", "lightpink", "purple", "yellow", "royalblue", "cyan","brown", "darkgreen", "greenyellow") ## china-metro-spore
#colorslegend <- c("purple", "orange1", "yellow", "lightgrey") ## others-metro cities
#colorslegend <- c("orange1","yellow","purple","lightgrey") ## others-dialect
#colorslegend <- c("red", "magenta", "lightpink", "purple", "yellow", "royalblue", "cyan","brown", "darkgreen", "greenyellow", "lightgrey") ## china-metro-spore+1000shanghai


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
#poplegend <- c("Liaoning", "Henan", "Hebei", "Shandong", "Anhui", "Hubei", "Jiangsu", "Hunan", "Guangdong") ## no shanghai, sichuan
#poplegend <- c("Liaoning", "Henan", "Hebei", "Shandong", "Anhui", "Hubei", "Jiangsu", "Shanghai", "Sichuan", "Hunan", "Guangdong", "CHB")
#poplegend <- c("Liaoning", "Henan", "Hebei", "Shandong", "Anhui", "Hubei", "Jiangsu", "Shanghai", "Sichuan", "Hunan", "Guangdong", "CHB", "Singapore")
#poplegend <- c("Liaoning", "Henan", "Hebei", "Shandong", "Anhui", "Hubei", "Jiangsu", "Shanghai", "Sichuan", "Hunan", "Guangdong-Cantonese", "Guangdong-Hakka", "Guangdong-Teochew","CHB", "Singapore")
poplegend <- c("Liaoning", "Hebei", "Shandong", "Henan", "Jiangsu", "Anhui", "Hubei", "Sichuan", "Hunan", "Guangdong") #china-metro-spore
#poplegend <- c("CHB", "Shanghai", "Singapore", "Others") ##others-metro cities
#poplegend <- c("Guangdong-Cantonese","Guangdong-Hakka","Guangdong-Teochew","Others") ##others-dialect
#poplegend <- c("Liaoning", "Hebei", "Shandong", "Henan", "Jiangsu", "Anhui", "Hubei", "Sichuan", "Hunan", "Guangdong", "Shanghai") #china-metro-spore+1000shanghai


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

		pca.plot(data[[xComp]], data[[yComp]], xlab=xLabel, ylab=yLabel , main=paste(title, " - ", population[i]), fgColors, data$population.abbreviation, i)
		#pca.plot(data[[xComp]], data[[yComp]], xlab=xLabel, ylab=yLabel , main="", fgColors, data$population.abbreviation, i)
	}
}

#pca.plot(data[[xComp]], data[[yComp]], xlab=xLabel, ylab=yLabel , main=title, colors, data$population.abbreviation)
pca.plot(data[[xComp]], -data[[yComp]], xlab=xLabel, ylab=yLabel, main='' , colors, data$population.abbreviation)
#pca.plot(-(data[[xComp]]), -(data[[yComp]]), xlab=xLabel, ylab=yLabel, main='' , colors, data$population.abbreviation)

x11()
pairs(data[,3:12], main="", pch = 20, col = colors[unclass(factor(data$population.abbreviation))])