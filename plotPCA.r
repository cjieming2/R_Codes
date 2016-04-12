### takes in sa and pca files only
### requires some changes with the colors for different sets of data

setwd("C:/Documents and Settings/chenjm/Desktop/work_documentation/chinese population analysis/shanghai/postqc w hapmap2r23")
pca.plot <- function (x, y, xlab, ylab, main, col, categories, lweight=40, rweight=10, ...)
{
	#laying the rules for graphwork: eg x11 is window
	x11(10,7, pointsize = 18)
	layout(matrix(c(rep(1,lweight), rep(2,rweight)), 1, 50, byrow = TRUE))
	layout.show(2)
	plot(x, y, xlab=xlab,ylab=ylab, main=main, pch=21, bg = colors[unclass(factor(categories))])
	
	## add line
	#abline(v=0.1, col = "red")
	#abline(h=0.15, col = "red")
	#text(x=0.15,y=-0.1, labels="PC3=0.1", col="red")
	#text(x=-0.05,y=0.16, labels="PC1=0.15", col="red")

	# highlight points in red
	#id <- which (x > -0.04 & x < -0.02 & data$province == "Shanghai")
	#id2 <- which (y > 0.05 & y < 0.1)
	#id3 <- which (j > 0.35)
	points(j[id],m[id], col="red")
	points(j[id2],m[id2], col="red")
	#points(j[id3],m[id3], col="blue")


	# add label sample-id on the points
	# pos 1-bottom 2-left 3-top 4-right
	# offset distance away from point
	#text(j[id], m[id], labels=data$sample.id[id], pos=4, offset=0.5, col="black")	
	#text(j[id2], m[id2], labels=data$sample.id[id2], pos=1, offset=0.5, col="black")	
	#text(j[id3], m[id3], labels=data$sample.id[id3], pos=2, adj=c(0,1.5), offset=0.5, col="blue")	
	
	par(mar = c(0,0,0,0))
	plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")

	legend(-1, 0.8, sort(unique(categories)), pt.bg = colors,
      	 text.col = "black", pch = 21,
        	 bg = 'white')
}

### DO NOT ADD ETHNICITY in SA FILE!!! cos it's got spaces, will screw up the read.table function
### annotation loading sa file
annot <- read.table("939samples-shanghai-hapmap2r23.sa", header=TRUE, sep="\t")

baseFile = "939samples-51707snps-shanghai-hapmap2r23-forward-snpset4"

### pca data loading pca file
filename= paste(baseFile, ".pca", sep="")
data <- read.table(filename, header=TRUE, sep="\t")

### merge annot and pca data together
data <- merge(annot, data, by='sample.id')

### include in the vaRIANCE info from the .eval file
evalFile = paste(baseFile, ".eval", sep="")
varInfo <- read.table(evalFile, header=TRUE, sep="\t")

### setting the colors
colors <- c("green", "red", "blue", "orange", "black", "pink", "purple","brown","violetred","slategray1","yellow","tan","cyan","deeppink")
#colors <- c("orange", "grey", "grey", "grey", "grey", "grey", "grey","grey","grey","grey","grey","tan","cyan","deeppink")

### last field of this func change to categories according to that category
xaxis="PC1"
yaxis="PC2"
j <- data$PC1
m <- data$PC2
xlabel=paste(xaxis, "(",varInfo$percentage.of.variance[which(varInfo$PC==xaxis)]*100,"%)")
ylabel=paste(yaxis, "(",varInfo$percentage.of.variance[which(varInfo$PC==yaxis)]*100,"%)")

#EX == pca.plot(data$PC3, data$PC1, xlab="PC1", ylab="PC2", main=filename, colors, data$population.abbreviation)
pca.plot(j, m, xlab=xlabel, ylab=ylabel, main=filename, colors, data$population.abbreviation)

### pairwise plots PC1 to 5 in one page
### change the data[,3:7]
### x11() creates another window
x11()
### note that col is border, bg is backgroundS of the points pch 21
### pch20 are bullets, no border, so only col 
pairs(data[,3:7], main="", pch = 20, col = colors[unclass(factor(data$population.abbreviation))])
##pairs(data[,3:7], main="", pch = 21, bg = colors[unclass(factor(data$province))])
