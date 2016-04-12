setwd("G:/home/tana2/Hepatitis B/20071120-pca/db123-1940samples-728snps-hepb-stage1-stage2-hapmap")

#annotation loading
sample_annotation <- read.table("1987samples-hepb-stage1-stage2-hapmap.sa", header=T, sep = "\t")

pca.plot <- function (x, y, xlab, ylab, main, col, categories, lweight=40, rweight=10, ...)
{
	x11(10,7, pointsize = 18)
	layout(matrix(c(rep(1,lweight), rep(2,rweight)), 1, 50, byrow = TRUE))
	layout.show(2)
	plot(x, y, xlab=xlab,ylab=ylab, main=main, pch = 21, bg = col[unclass(factor(categories))])
	par(mar = c(0,0,0,0))
	plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
	legend(-1, 0.8, sort(unique(categories)), pt.bg = col,
      	 text.col = "black", pch = 21,
        	 bg = 'white')
}

########################################
#Hepatitis B Stage 1 - Stage 2 / Hapmap#
########################################
data <- read.table("db123-1940samples-728snps-hepb-stage1-stage2-hapmap.pca", header=T, sep = "\t")
data <- merge(sample_annotation, data, by = "sample.id")

col <- c("yellow", "blue", "black", "green", "purple");

PC      eigenvalue      percentage-of-variance
PC1     27.180920       0.027388
PC2     17.990905       0.018128
PC3     15.673300       0.015793
PC4     15.245204       0.015361
PC5     11.056808       0.011141
PC6     10.861598       0.010944
PC7     10.216692       0.010295
PC8     9.902375        0.009978
PC9     9.261326        0.009332

pca.plot(data$PC1, data$PC2, xlab="PC1 (2.74%)", ylab="PC2 (1.81%)", main="Hepatitis B Stage 1 - Stage 2 - Hapmap", col, data$population.id)
pca.plot(data$PC1, data$PC3, xlab="PC1 (2.74%)", ylab="PC3 (1.58%)", main="Hepatitis B Stage 1 - Stage 2 - Hapmap", col, data$population.id)
pca.plot(data$PC1, data$PC4, xlab="PC1 (2.74%)", ylab="PC4 (1.54%)", main="Hepatitis B Stage 1 - Stage 2 - Hapmap", col, data$population.id)

col <- c("yellow", "blue", "black", "grey", "green", "purple");

pca.plot(data$PC1, data$PC2, xlab="PC1 (2.74%)", ylab="PC2 (1.81%)", main="Hepatitis B Stage 1 - Stage 2 - Hapmap", col, data$stage)
pca.plot(data$PC1, data$PC3, xlab="PC1 (2.74%)", ylab="PC3 (1.58%)", main="Hepatitis B Stage 1 - Stage 2 - Hapmap", col, data$stage)
pca.plot(data$PC1, data$PC4, xlab="PC1 (2.74%)", ylab="PC4 (1.54%)", main="Hepatitis B Stage 1 - Stage 2 - Hapmap", col, data$stage)



## put histograms on the diagonal
panel.hist <- function(x, ...)
{
   usr <- par("usr"); on.exit(par(usr))
   par(usr = c(usr[1:2], 0, 1.5) )
   h <- hist(x, plot = FALSE)
   breaks <- h$breaks; nB <- length(breaks)
   y <- h$counts; y <- y/max(y)
   rect(breaks[-nB], 0, breaks[-1], y, col="cyan")
}

panel.smooth <- function (x, y, col = par("col"), bg = NA, pch = par("pch"),
   cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...)
{
   col1 = c("yellow", "blue", "black", "green", "purple")
   points(x, y, pch = ".", col = col1[unclass(data$population.id)], bg = bg, cex = cex)
   ok <- is.finite(x) & is.finite(y)
   if (any(ok))
       lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),
           col = col.smooth, ...)
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")

    text(0.5, 0.5, txt, cex = 0.5)
}

panel.text <- function(x, y, labels, cex, font, ...) 
{
    text(0.5, 0.8, labels, cex = 1)
}

x11()
pairs(data[,7:16], diag.panel=panel.hist, lower.panel=panel.smooth,  upper.panel=panel.smooth, text.panel=panel.text, main="Hepatitis B Stage 1 - Stage 2 - Hapmap (1940 individuals, 728 snps)") 
