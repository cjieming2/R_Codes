#takes in sa and pca files only
#requires some changes with the colors for different sets of data

#setwd("C:/Documents and Settings/chenjm/Desktop/work_documentation/Mark/indian-reich-hapmap-pasnp-ihp/pca")
#setwd("C:/Documents and Settings/chenjm/Desktop/work_documentation/JJ/chinese population analysis/manuscript/ajhg submission/review1/pca/east-west")
setwd("C:/Documents and Settings/chenjm/Desktop/work_documentation/Paul/w-hhp-jnov-QC-new")
#setwd("C:/Documents and Settings/chenjm/Desktop/work_documentation/popres/john novembre")

#pca.plot function
pca.plot <- function (x, y, xlab, ylab, main, col, categories, highlightIdx=0, lweight=40, rweight=10, ...)
{
	#laying the rules for graphwork: eg x11 is font
	x11(16,9, pointsize = 18)
	layout(matrix(c(rep(1,lweight), rep(2,rweight)), 1, 50, byrow = TRUE))
	layout.show(2)

	if(highlightIdx == 0)
	{
		#plot(x, y, xlab=xlab,ylab=ylab, main=main, pch=21, col="black", bg=col[unclass(factor(categories))])
		plot(x, y, xlab=xlab,ylab=ylab, main=main, pch=20, bg="lightgrey", col=col[unclass(factor(categories))])
		
		## brings forward the points desired
		#id = which(data$population.abbreviation == "YRI")
		#id2 = which(data$population.abbreviation == "Japan")
		#id3 = which(data$population.abbreviation == "JPT")
		#id4 = which(data$population.abbreviation == "CHB")
		#id5 = which(data$population.abbreviation == "CEU")
		id6 = which(data$population.abbreviation == "hESC")
		#id6 = which(data$sample.id == "NN10")

		#points(x[id], y[id], pch=20, col="red")
		#points(x[id2], y[id2], pch=21, col="black", bg="cyan")
		#points(x[id3], y[id3], pch=20, col="yellow")
		#points(x[id4], y[id4], pch=20, col="orange")
		#points(x[id5], y[id5], pch=20, col="cyan")
		points(x[id6], y[id6], pch=21, col="red", bg="black")


		#id10 = which(data$population.abbreviation != "CEU" & data$population.abbreviation != "JPT" & data$population.abbreviation != "CHB" & data$population.abbreviation != "YRI")
		#points(x[id10], y[id10], pch=21, col="black", bg="pink")
		text(x[id6], y[id6], labels=data$sample.id[id6], pos=1, adj=c(0,1.5), offset=0.5, col="black")
		



		
	}
	else
	{
		plot(x, y, xlab=xlab,ylab=ylab, main=main, pch=19, col=col[unclass(factor(categories))])
		id = which(data$population.abbreviation == population[highlightIdx])
		#id = which(data$dataset == population[highlightIdx])
		points(x[id], y[id], pch=21, col="black", bg=colors[highlightIdx])
		id6 = which(data$population.abbreviation == "hESC")
		points(x[id6], y[id6], pch=21, col="red", bg="black")
		text(x[id6], y[id6], labels=data$sample.id[id6], pos=1, adj=c(0,1.5), offset=0.5, col="black")


	}
	
	## add line
	#abline(v=0.0564, col = "red")
	#abline(v=-0.018, col = "red")
	#abline(v=-0.04, col = "red")
	#abline(h=-0.005, col = "red")
	#text(x=0.045,y=-0.2, labels="PC2=0.0564", col="red")
	#text(x=-0.01,y=-0.2, labels="PC2=-0.018", col="red")
	#text(x=-0.05,y=-0.2, labels="PC2=-0.04", col="red")
	
	# highlight points in red
	#id <- which (data$sample.id == "AD788" | data$sample.id == "CH-SC6" | data$sample.id == "AD927" | data$sample.id == "AD1140" | data$sample.id == "AD898" | data$sample.id == "AD829")
	#id2 <- which (data$population.abbreviation == "")
	#id2 <- which (data$population.abbreviation == "California, USA" | data$population.abbreviation == "Chicago, USA" | data$population.abbreviation == "India" | data$population.abbreviation == "Iran" | data$population.abbreviation == "Korea" | data$population.abbreviation == "Sydney, Australia" | data$population.abbreviation == "USA" | data$population.abbreviation == "Czech Republic" | data$population.abbreviation == "Edinburgh, UK" | data$population.abbreviation == "Finland" | data$population.abbreviation == "Netherlands" | data$population.abbreviation == "Russia" | data$population.abbreviation == "Sweden")
	#id3 <- which (data$population.abbreviation == "hESC")


	#id3 <- which ((data$population.abbreviation == "Escells" & data$PC1 < -0.03))
	#id4 <- which (data$PC5 > 0.2 & data$PC5 < 0.4)
	#id5 <- which (data$sample.id == "HGDP01031" | data$sample.id == "HGDP0117")
	#points(x[id], y[id], col="red")
	#points(x[id2], y[id2], col="red")
	#points(x[id3], y[id3], col="red")
	#points(x[id4], y[id4], col="red")
	#points(x[id5], y[id5], col="cyan")

	# add label sample-id on the points
	# pos 1-bottom 2-left 3-top 4-right
	# offset distance away from point
	#text(x[id], y[id], labels=data$sample.id[id], pos=2, offset=0.5, col="red")	
	#text(x[id2], y[id2], labels=data$sample.id[id2], pos=2, adj=c(1,1.1), offset=0.5, col="red")	
	#text(x[id3], y[id3], labels=data$sample.id[id3], pos=2, adj=c(1,2), offset=0.5, col="black")	
	#text(x[id4], y[id4], labels=data$sample.id[id4], pos=2, adj=c(0,1.5), offset=0.5, col="black")	
	#text(x[id5], y[id5], labels=data$sample.id[id5], pos=1, adj=c(0,1.5), offset=0.5, col="cyan")

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
baseFileName = "1464samples-55972snps-merged-es-jnov"		#Base file name of the .PCA and .EVAL file from fpca
title = baseFileName	#Title of plot
annotFilename="1464samples-merged-es-jnov.sa" ##Annotation File

colors <- c("cyan",
"orange", 
"black", 
"pink", 
"purple", 
"green", 
"darkolivegreen2", 
"blue",
"brown",
"magenta",
"yellow",
"darkgreen",
"tan",
#"lightseagreen", 
"red") ##Set of 14 colours to use for all plots

colors <- c("red", 
"blue", 
"black",
"orange", 
"coral3", 
"green", 
"olivedrab",
"magenta",
"yellow",
"cyan",
#"darkgreen",
"tan",
"tomato", 
#"deeppink",
"royalblue", 
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

#novembre original
colors <- c("lightblue","brown","red","blue","blue","blue","red","green","lightblue","darkgreen","orange","orange","blue","orange","red","brown","blue","green","lightblue","orange","lightblue","lightgreen","magenta","darkgreen","green","pink","lightblue","brown","blue","purple","lightseagreen","red","red","red","green","green","pink","blue")
#novembre w hESC & swissfrench-german-italian
colors <- c("lightblue","brown","red","blue","blue","blue","red","green","lightblue","darkgreen","orange","orange","blue","black","orange","red","brown","blue","green","lightblue","orange","lightblue","lightgreen","magenta","darkgreen","green","pink","lightblue","brown","blue","purple","lightseagreen","red","green","green","pink","blue")
#nov-switzerland,britishisle,
#colors <- c("blue","orange","cyan","blue","blue","blue","red","green","lightseagreen","darkgreen","orange","orange","blue","black","yellow","pink","brown","blue","green","blue","red","lightseagreen","lightgreen","magenta","blue","green","pink","blue","navajowhite2","blue","purple","lightseagreen","orange","cyan","green","pink","blue")


#colorslegend <- c("black","cyan","orange","yellow","red","pink","purple","blue","brown","magenta","green","darkolivegreen2","darkgreen","tan","lightseagreen")
#colorslegend <- c("black","cyan","orange","yellow","red","pink","purple","blue","brown","magenta","green","darkolivegreen2","darkgreen","tan")


pcaFilename = paste(baseFileName, ".pca", sep="")	#Change this if .pca file name diff from base file name
evalFilename= paste(baseFileName, ".eval", sep="")	#Change this if .eval file name diff from base file name
##----- End set up ----------------

#DO NOT ADD ETHNICITY in SA FILE!!! cos it's got spaces, will screw up the read.table function
#annotation loading sa file
annot <- read.table(annotFilename, header=TRUE, sep="\t", stringsAsFactors = F)

#pca data loading pca file
data <- read.table(pcaFilename, header=TRUE, sep="\t", stringsAsFactors = F)


#merge annot and pca data together
data <- merge(annot, data, by='sample.id')

#read eval file
#evalFilename="db126-1293samples-45CHBsamples-7454snps-conserved-region.eval"
amtVariance <- read.table(evalFilename, header=TRUE, sep="\t")

#Get a sorted unique set of populations and number of populations
population <- sort(unique(data$population.abbreviation))
#population <- sort(unique(data$dataset))
#poplegend <- c("hESC","CEU","CHB","JPT","YRI","HGDP-Africa","HGDP-America","HGDP-Europe","HGDP-Middle-East","HGDP-Oceania","HGDP-CentralSouth-Asia","HGDP-East-Asia","PASNP-CentralSouth-Asia","PASNP-East-Asia","PASNP-Southeast-Asia")
#poplegend <- c("hESC","CEU","CHB","JPT","YRI","HGDP-Africa","HGDP-America","HGDP-Europe","HGDP-Middle-East","HGDP-Oceania","HGDP-CentralSouth-Asia","HGDP-East-Asia","PASNP-CentralSouth-Asia","PASNP-East-Asia")
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
		pca.plot(data[[xComp]], data[[yComp]], xlab=xLabel, ylab=yLabel , main=paste(title, " - ", population[i]), fgColors, data$population.abbreviation, i)
		#pca.plot(data[[xComp]], data[[yComp]], xlab=xLabel, ylab=yLabel , main="", fgColors, data$dataset, i)
	}
}

#pca.plot(data[[xComp]], data[[yComp]], xlab=xLabel, ylab=yLabel , main=title, colors, data$dataset)
#pca.plot(data[[xComp]], data[[yComp]], xlab=xLabel, ylab=yLabel, main=title , colors, data$population.abbreviation)
pca.plot(-(data[[xComp]]), (data[[yComp]]), xlab=xLabel, ylab=yLabel, main=title , colors, data$population.abbreviation)

#x11()
#pairs(data[,13:23], main="", pch = 20, col = colors[unclass(factor(data$population.abbreviation))])
#pairs(data[,5:8], main="", pch = 20, col = colors[unclass(factor(data$population.abbreviation))])
#pairs(data[,4:8], main="", pch = 20, col = colors[unclass(factor(data$population.abbreviation))])


