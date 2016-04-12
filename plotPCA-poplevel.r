### takes in sa and pca files only
### requires some changes with the colors for different sets of data

setwd("C:/Documents and Settings/chenjm/Desktop/work_documentation/hapmap/hapmap 2 r26/pca-pop")
pca.plot <- function (x, y, xlab, ylab, main, col, categories, lweight=40, rweight=10, ...)
{
	#laying the rules for graphwork: eg x11 is window
	x11(10,7, pointsize = 30)
	layout(matrix(c(rep(1,lweight), rep(2,rweight)), 1, 50, byrow = TRUE))
	layout.show(2)
	plot(x, y, xlab=xlab,ylab=ylab, main=main, pch=21, bg = colors[unclass(factor(categories))])
	
	# highlight points in red
	#id <- which (m > -0.02 & m < -0.018 & j < -0.005 & j > -0.01)
	#id2 <- which (j < -0.05 & m > 0.12)
	#id3 <- which (j > 0.1 & m < -0.4)
	#points(data$PC1[id],data$PC2[id], col="red")

	# add label sample-id on the points
	# pos 1-bottom 2-left 3-top 4-right
	# offset distance away from point
	#text(data$PC1[id], data$PC2[id], labels=data$sample.id[id], pos=2, offset=0.5, col="black")	
	#text(j[id2], m[id2], labels=data$sample.id[id2], pos=4, offset=0.5, col="black")	
	#text(j[id3], m[id3], labels=data$sample.id[id3], pos=2, adj=c(0,1.5), offset=0.5, col="black")	
	
	par(mar = c(0,0,0,0))
	plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
	
	## cex changes font size; it's a multiplier, so decimal shrinks
	legend(-1, 0.8, sort(unique(categories)), pt.bg = colors,
      	 text.col = "black", pch = 21, cex=0.8,
        	 bg = 'white')
}

### DO NOT ADD ETHNICITY in SA FILE!!! cos it's got spaces, will screw up the read.table function
### annotation loading sa file
#annot <- read.table("462samples-hs-r26.sa", header=TRUE, sep="\t")

### pca data loading pca file
filename="renamed-db126-196samples-2896293snps-hapmap2-r26-forward-founders-only.pca"
data <- read.table(filename, header=TRUE, sep="\t")

### merge annot and pca data together
#data <- merge(annot, data, by='sample.id')

### include in the vaRIANCE info from the .eval file
varInfo <- read.table("renamed-db126-196samples-2896293snps-hapmap2-r26-forward-founders-only.eval", header=TRUE, sep="\t")

### setting the colors
colors <- c("orange", "brown", "blue", "green", "black", "red", "purple","pink","violetred","slategray1","yellow","tan","cyan","deeppink")
#colors <- c("orange", "grey", "grey", "grey", "grey", "grey", "grey","grey","grey","grey","grey","tan","cyan","deeppink")

### last field of this func change to categories according to that category
xaxis="PC1"
yaxis="PC2"
j <- data$PC1
m <- data$PC2
xlabel=paste(xaxis, "(",varInfo$percentage.of.variance[which(varInfo$PC==xaxis)]*100,"%)")
ylabel=paste(yaxis, "(",varInfo$percentage.of.variance[which(varInfo$PC==yaxis)]*100,"%)")

#EX == pca.plot(data$PC1, data$PC2, xlab="PC1", ylab="PC2", main=filename, colors, data$population.id)
pca.plot(j, m, xlab=xlabel, ylab=ylabel, main=filename, colors, data$population.id)

### pairwise plots PC1 to 5 in one page
### change the data[,3:7]
### x11() creates another window
x11()
### note that col is border, bg is backgroundS of the points pch 21
### pch20 are bullets, no border, so only col 
pairs(data[,2:6], main="", pch = 20, col = colors[unclass(factor(data$population.id))])
##pairs(data[,3:7], main="", pch = 21, bg = colors[unclass(factor(data$province))])
