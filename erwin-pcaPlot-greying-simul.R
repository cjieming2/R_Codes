#################################################################
### Read Sample Annotation
path <- "E:/Erwin/Report-Consolidation/Population-Analysis"
setwd(path)

sample_annotation <- read.table("7274samples.sa", header=T, sep = "\t")

### Read data
infile <- "db126-7274samples.pca";
data <- read.table(infile, header=T, sep = "\t")
data <- merge(sample_annotation, data, by = "sample.id")


col1 <- rainbow(6)
col <- c(col1, "grey", "yellowgreen", "brown", "black", "orange", "white")

anhuicol <- c(col1[1], "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey")
chbcol <- c("lightgrey", col1[2], "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey")
guangzhoucol <- c("lightgrey", "lightgrey", col1[3], "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey")
hebeicol <- c("lightgrey", "lightgrey", "lightgrey", col1[4], "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey")
henancol <- c("lightgrey", "lightgrey", "lightgrey", "lightgrey", col1[5], "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey")
hubeicol <- c("lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", col1[6], "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey")
hunancol <- c("lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "grey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey")
jiangsucol <- c("lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "yellowgreen", "lightgrey", "lightgrey", "lightgrey", "lightgrey")
liaoningcol <- c("lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "brown", "lightgrey", "lightgrey", "lightgrey")
shandongcol <- c("lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "black", "lightgrey", "lightgrey")
singaporecol <- c("lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "orange", "lightgrey")

### Guangzhou
x11(10,7, pointsize = 20)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
layout.show(2)
myprovince <- which(data$province=="Guangzhou")
mycolor <- col1[3]
plot(data$PC1, data$PC2, xlab="PC1", ylab="PC2", main="PC1 vs PC2", pch = 20, col = guangzhoucol[unclass(factor(data$province))])
points(data$PC1[myprovince], data$PC2[myprovince], pch=21, bg=mycolor)

### highlight points in red
j <- data$PC1
m <- data$PC2
id <- which (data$province=="Guangzhou" & j > -0.0025)
abline(v=-0.0025, col="red")
text(0.00,-0.02,labels="x=-0.0025", col="red")
points(data$PC1[id],data$PC2[id], col="red")

### add label sample-id on the points
### pos 1-bottom 2-left 3-top 4-right
### offset distance away from point
text(j[id], m[id], labels=data$sample.id[id], pos=3, offset=0.5, col="black")	

par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(factor(data$province))), col = guangzhoucol, text.col = "black", pch = 20, bg = 'white')

### Shandong
x11(10,7, pointsize = 20)
layout(matrix(c(1,1,1,1,2), 1, 5, byrow = TRUE))
layout.show(2)
myprovince <- which(data$province=="Shandong")
mycolor <- "black"
plot(data$PC1, data$PC2, xlab="PC1", ylab="PC2", main="PC1 vs PC2", pch = 20, col = shandongcol[unclass(factor(data$province))])
points(data$PC1[myprovince], data$PC2[myprovince], pch=21, bg=mycolor)

### highlight points in red
j <- data$PC1
m <- data$PC2
id <- which (data$province=="Shandong" & j < 0.007)
abline(v=0.007, col="red")
text(0.005,-0.02,labels="x=0.007",col="red")
points(data$PC1[id],data$PC2[id], col="red")

### add label sample-id on the points
### pos 1-bottom 2-left 3-top 4-right
### offset distance away from point
text(j[id], m[id], labels=data$sample.id[id], pos=3, offset=0.5, col="black")	

par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(factor(data$province))), col = shandongcol, text.col = "black", pch = 20, bg = 'white')

