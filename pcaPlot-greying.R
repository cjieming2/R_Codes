#################################################################
### Read Sample Annotation
path <- "C:/Documents and Settings/chenjm/Desktop/work_documentation/hsa/6 pca/"
setwd(path)

sample_annotation <- read.table("327samples-CEU+SGVP.sa", header=T, sep = "\t")

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
infile <- "ceu+sgvp-sieved-127snp-sgvp-hapmap_F5.mk.pca";
data <- read.table(infile, header=T, sep = "\t")
data <- merge(sample_annotation, data, by = "sample.id")

### include in the vaRIANCE info from the .eval file
varInfo <- read.table("ceu+sgvp-sieved-127snp-sgvp-hapmap_F5.mk.eval", header=TRUE, sep="\t")
xaxis="PC1"
yaxis="PC2"
xlabel=paste(xaxis, "(",varInfo$percentage.of.variance[which(varInfo$PC==xaxis)]*100,"%)")
ylabel=paste(yaxis, "(",varInfo$percentage.of.variance[which(varInfo$PC==yaxis)]*100,"%)")

### colors of each
col1 <- c("orange", "green", "red", "purple", "black", "violetred", "pink","brown","blue","slategray1","yellow","tan","cyan","deeppink")
col <- c(col1, "grey", "yellowgreen", "brown", "black", "orange", "white")

ceusgvpcol <- c(col1[1], col1[2], "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey")
#chbcol <- c("lightgrey", col1[2], "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey")
#guangzhoucol <- c("lightgrey", "lightgrey", col1[3], "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey")
#hebeicol <- c("lightgrey", "lightgrey", "lightgrey", col1[4], "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey")
#henancol <- c("lightgrey", "lightgrey", "lightgrey", "lightgrey", col1[5], "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey")
#hubeicol <- c("lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", col1[6], "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey")
#hunancol <- c("lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", col1[7], "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey")
#jiangsucol <- c("lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", col1[8], "lightgrey", "lightgrey", "lightgrey", "lightgrey")
#liaoningcol <- c("lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", col1[9], "lightgrey", "lightgrey", "lightgrey")
#shandongcol <- c("lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", col1[10], "lightgrey", "lightgrey")
#singaporecol <- c("lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", col1[11], "lightgrey")

### Anhui
x11(10,7, pointsize = 20)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
layout.show(2)
myprovince <- which(data$province=="Anhui")
mycolor <- col1[1]
plot(data$PC1, data$PC2, xlab=xlabel, ylab=ylabel, main="PC1 vs PC2", pch = 20, col = anhuicol[unclass(factor(data$province))])
points(data$PC1[myprovince], data$PC2[myprovince], pch=21, bg=mycolor)
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(factor(data$province))), col = anhuicol, text.col = "black", pch = 20, bg = 'white')

### CHB
x11(10,7, pointsize = 20)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
layout.show(2)
myprovince <- which(data$province=="CHB")
mycolor <- col1[2]
plot(data$PC1, data$PC2, xlab=xlabel, ylab=ylabel, main="PC1 vs PC2", pch = 20, col = chbcol[unclass(factor(data$province))])
points(data$PC1[myprovince], data$PC2[myprovince], pch=21, bg=mycolor)
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(factor(data$province))), col = chbcol, text.col = "black", pch = 20, bg = 'white')


### Guangzhou
x11(10,7, pointsize = 20)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
layout.show(2)
myprovince <- which(data$province=="Guangzhou")
mycolor <- col1[3]
plot(data$PC1, data$PC2, xlab=xlabel, ylab=ylabel, main="PC1 vs PC2", pch = 20, col = guangzhoucol[unclass(factor(data$province))])
points(data$PC1[myprovince], data$PC2[myprovince], pch=21, bg=mycolor)
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(factor(data$province))), col = guangzhoucol, text.col = "black", pch = 20, bg = 'white')

### Hebei
x11(10,7, pointsize = 20)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
layout.show(2)
myprovince <- which(data$province=="Hebei")
mycolor <- col1[4]
plot(data$PC1, data$PC2, xlab=xlabel, ylab=ylabel, main="PC1 vs PC2", pch = 20, col = hebeicol[unclass(factor(data$province))])
points(data$PC1[myprovince], data$PC2[myprovince], pch=21, bg=mycolor)
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(factor(data$province))), col = hebeicol, text.col = "black", pch = 20, bg = 'white')

### Henan
x11(10,7, pointsize = 20)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
layout.show(2)
myprovince <- which(data$province=="Henan")
mycolor <- col1[5]
plot(data$PC1, data$PC2, xlab=xlabel, ylab=ylabel, main="PC1 vs PC2", pch = 20, col = henancol[unclass(factor(data$province))])
points(data$PC1[myprovince], data$PC2[myprovince], pch=21, bg=mycolor)
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(factor(data$province))), col = henancol, text.col = "black", pch = 20, bg = 'white')

### Hubei
x11(10,7, pointsize = 20)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
layout.show(2)
myprovince <- which(data$province=="Hubei")
mycolor <- col1[6]
plot(data$PC1, data$PC2, xlab=xlabel, ylab=ylabel, main="PC1 vs PC2", pch = 20, col = hubeicol[unclass(factor(data$province))])
points(data$PC1[myprovince], data$PC2[myprovince], pch=21, bg=mycolor)
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(factor(data$province))), col = hubeicol, text.col = "black", pch = 20, bg = 'white')

### Hunan
x11(10,7, pointsize = 20)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
layout.show(2)
myprovince <- which(data$province=="Hunan")
mycolor <- col1[7]
plot(data$PC1, data$PC2, xlab=xlabel, ylab=ylabel, main="PC1 vs PC2", pch = 20, col = hunancol[unclass(factor(data$province))])
points(data$PC1[myprovince], data$PC2[myprovince], pch=21, bg=mycolor)
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(factor(data$province))), col = hunancol, text.col = "black", pch = 20, bg = 'white')

### Jiangsu
x11(10,7, pointsize = 20)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
layout.show(2)
myprovince <- which(data$province=="Jiangsu")
mycolor <- col1[8]
plot(data$PC1, data$PC2, xlab=xlabel, ylab=ylabel, main="PC1 vs PC2", pch = 20, col = jiangsucol[unclass(factor(data$province))])
points(data$PC1[myprovince], data$PC2[myprovince], pch=21, bg=mycolor)
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(factor(data$province))), col = jiangsucol, text.col = "black", pch = 20, bg = 'white')

### Liaoning
x11(10,7, pointsize = 20)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
layout.show(2)
myprovince <- which(data$province=="Liaoning")
mycolor <- col1[9]
plot(data$PC1, data$PC2, xlab=xlabel, ylab=ylabel, main="PC1 vs PC2", pch = 20, col = liaoningcol[unclass(factor(data$province))])
points(data$PC1[myprovince], data$PC2[myprovince], pch=21, bg=mycolor)
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(factor(data$province))), col = liaoningcol, text.col = "black", pch = 20, bg = 'white')

### Shandong
x11(10,7, pointsize = 20)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
layout.show(2)
myprovince <- which(data$province=="Shandong")
mycolor <- col1[10]
plot(data$PC1, data$PC2, xlab=xlabel, ylab=ylabel, main="PC1 vs PC2", pch = 20, col = shandongcol[unclass(factor(data$province))])
points(data$PC1[myprovince], data$PC2[myprovince], pch=21, bg=mycolor)
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(factor(data$province))), col = shandongcol, text.col = "black", pch = 20, bg = 'white')

### Singapore
x11(10,7, pointsize = 20)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
layout.show(2)
myprovince <- which(data$province=="Singapore")
mycolor <- col1[11]
plot(data$PC1, data$PC2, xlab=xlabel, ylab=ylabel, main="PC1 vs PC2", pch = 20, col = singaporecol[unclass(factor(data$province))])
points(data$PC1[myprovince], data$PC2[myprovince], pch=21, bg=mycolor)
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(factor(data$province))), col = singaporecol, text.col = "black", pch = 20, bg = 'white')


#newcol <- c(col1[1], col1[2], col1[3], "black", "orange")
#x11()
#pairs(data[,3:7],  main="", pch = 20, col = newcol[unclass(factor(data$province))])

