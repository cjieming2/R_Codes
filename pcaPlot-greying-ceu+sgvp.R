#################################################################
### Read Sample Annotation
path <- "C:/Documents and Settings/chenjm/Desktop/work_documentation/hsa/6 pca/"
setwd(path)

sample_annotation <- read.table("327samples-CEU+SGVP.sa", header=T, sep = "\t")
baseFileName = "ceu+sgvp-sieved-793snp-sgvp-hapmap_DPYD.mk"

### highlight points in red
#id <- which (j < -0.1 & m < 0.05 & m > 0)
#id2 <- which (j < -0.05 & m > 0.12)
#id3 <- which (j > 0.1 & m < -0.4)
#points(data$PC1[id],data$PC2[id], col="red")

### add label sample-id on the points
### pos 1-bottom 2-left 3-top 4-right
### offset distance away from point
#text(j[id], m[id], labels=data$sample.id[id], pos=4, offset=0.5, col="black")	
#text(j[id2], m[id2], labels=data$sample.id[id2], pos=4, offset=0.5, col="black")	
#text(j[id3], m[id3], labels=data$sample.id[id3], pos=2, adj=c(0,1.5), offset=0.5, col="black")	

### Read data
pcaFile = paste(baseFileName, ".pca", sep="")
data <- read.table(pcaFile, header=T, sep = "\t")
data <- merge(sample_annotation, data, by = "sample.id")

### include in the vaRIANCE info from the .eval file
evalFile = paste(baseFileName, ".eval", sep="")
varInfo <- read.table(evalFile, header=TRUE, sep="\t")
xaxis="PC1"
yaxis="PC2"
xlabel=paste(xaxis, "(",varInfo$percentage.of.variance[which(varInfo$PC==xaxis)]*100,"%)")
ylabel=paste(yaxis, "(",varInfo$percentage.of.variance[which(varInfo$PC==yaxis)]*100,"%)")

### colors of each
col1 <- c("red", "green", "orange", "purple", "black", "violetred", "pink","brown","blue","slategray1","yellow","tan","cyan","deeppink")
col <- c(col1, "grey", "yellowgreen", "brown", "black", "orange", "white")

ceusgvpcol <- c(col1[1], col1[2], "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey")
#chbcol <- c("lightgrey", col1[2], "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey")

### ceu+sgvp
x11(10,7, pointsize = 20)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
layout.show(2)
plot(data$PC1, data$PC2, xlab=xlabel, ylab=ylabel, main=pcaFile, pch = 20, col = ceusgvpcol[unclass(factor(data$population.abbreviation))])

myceu <- which(data$population.abbreviation=="CEU")
mysgvpch <- which(data$population.abbreviation=="SGVP-CH")
myceucolor <- col1[1]
mysgvpchcolor <- col1[2]
points(data$PC1[myceu], data$PC2[myceu], pch=21, bg=myceucolor)
points(data$PC1[mysgvpch], data$PC2[mysgvpch], pch=21, bg=mysgvpchcolor)

par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(factor(data$population.abbreviation))), col = ceusgvpcol, text.col = "black", pch = 20, bg = 'white')

newcol <- c(col1[1], col1[2], "lightgrey", "lightgrey", "lightgrey")
x11()
pairs(data[,5:9],  main="", pch = 20, col = newcol[unclass(factor(data$population.abbreviation))])

