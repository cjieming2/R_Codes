setwd("C:/Documents and Settings/chenjm/Desktop/")
ibs.plot <- function (x, y, xlab, ylab, main, col, categories, lweight=40, rweight=10, ...)
{
	#laying the rules for graphwork: eg x11 is font
	x11(10,7, pointsize = 18)
	layout(matrix(c(rep(1,lweight), rep(2,rweight)), 1, 50, byrow = TRUE))
	layout.show(2)
	plot(x, y, xlab=xlab,ylab=ylab, main=main, pch=21, col = colors[unclass(factor(categories))])
	
	### identify
	#id <- which (y < 0.59)
	id <- which (data$relationship == "FS" | data$relationship == "AV_HS")
	id2 <- which (data$relationship == "CO_GG")
	#points(x[id], y[id], col = "green")
	#text(x[id], y[id], labels=data$sample.pair.id[id],pos=2, offset=0.5)
	#text(x[id2], y[id2], labels=data$sample.pair.id[id2],pos=1, offset=0.5)


	par(mar = c(0,0,0,0))
	plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
	legend(-1, 0.8, sort(unique(categories)), pt.bg = colors,
      	 text.col = "black", pch = 21,
        	 bg = 'white')
}

#no sample id pair pls
#col1 - ibs mean 
#col2 - ibs sd
#col3 - relationship (different from rss)
data <- read.table("data.sar", header=TRUE, sep="\t")

#col is border, bg is background
colors <- c("black", "red", "blue", "purple", "green", "magenta", "lightblue1","orange")

ibs.plot(data$ibs.mean, data$ibs.stdev, xlab="IBS mean", ylab="IBS standard deviation", main="", colors, data$relationship)