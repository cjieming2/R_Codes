setwd("C:/Documents and Settings/chenjm/Desktop")

## A function to use identify to select points, and overplot the
## points with another symbol as they are selected
identifyPch <- function(x, y=NULL, n=length(x), pch=19, tag=data$sample.id, sizeofpoint = 3,...)
{
    xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
    sel <- rep(FALSE, length(x)); res <- integer(0)
    while(sum(sel) < n) {
        ans <- identify(x[!sel], y[!sel],labels = as.character(tag), col="red", ...)
        if(!length(ans)) break
        ans <- which(!sel)[ans]
        points(x[ans], y[ans], pch = pch, cex = sizeofpoint)
        sel[ans] <- TRUE
        res <- c(res, ans)
    }
    res
}

#this plots the call rates VS sample/snp heterozygotes
#low call rates high hetero means perhaps bad genotype calling/sample contmaination/excess heterozygosity
### PLEASE STILL DOUBLE CHECK THE VERACITY OF THE SAMPLES IN THE FILES!!!

filename="samples.txt"
data <- read.table(filename, header=TRUE, sep="\t")
#hetero <- (data$A_B/(data$A_A + data$A_B + data$B_B + data$N_N))*100
hetero <- (data$A_B/(data$A_A + data$A_B + data$B_B))*100
cr <- (data$sample.call.rate) * 100
#### rem to change ?_? to N_N

#col is border, bg is background
colors <- c("black", "orange", "blue", "red", "green", "magenta", "lightblue1")

numsamples <- length(hetero)
mainname = paste("zoomed samplecr VS heterozygosity plot for", numsamples," samples")
#zoomed plot
x11()
plot(hetero, cr, xlab="heterozygosity (%)", ylab="sample-call-rate (%)", main=mainname, pch=21, col="black")
identifyPch(hetero, cr, n=length(hetero), pch=4, col="red", tag=data$sample.id, sizeofpoint=2)

# add red line
#abline(v=26, col="red")

#scale of 0-100%
x11()
mainname = paste("samplecr VS heterozygosity plot for", numsamples," samples")
plot(hetero, cr, xlab="heterozygosity (%)", ylab="sample-call-rate (%)", main=mainname, pch=21, col="black", xlim=c(0,100), ylim=c(0,100))
identifyPch(hetero, cr, n=length(hetero), pch=4, col="red", tag=data$sample.id, sizeofpoint=2)

x11()
mainname = paste("samplecr VS heterozygosity plot for", numsamples," samples")
plot(hetero, cr, xlab="heterozygosity (%)", ylab="sample-call-rate (%)", main=mainname, pch=21, col="black", ylim=c(0,100))
identifyPch(hetero, cr, n=length(hetero), pch=4, col="red", tag=data$sample.id, sizeofpoint=2)

# identify desired points
# this stores the ID of the points in vector/matrix id
#id <- which (hetero < 25.5)
#id2 <- which (hetero < 29)
# highlight points in red
#points(hetero[id],cr[id], col="red")
#points(hetero[id2],cr[id2], col="red")
# add label sample-id on the points
# pos 1-bottom 2-left 3-top 4-right
# offset distance away from point
#text(hetero[id], cr[id], labels=data$sample.id[id], pos=4, offset=0.5, col="red")
#text(hetero[id2], cr[id2], labels=data$sample.id[id2], pos=4, offset=0.5, col="red")

x11()
mainname = paste("samplecr VS heterozygosity plot for", numsamples," samples")
plot(hetero, cr, xlab="heterozygosity (%)", ylab="sample-call-rate (%)", main=mainname, pch=21, col="black", xlim=c(0,100))
identifyPch(hetero, cr, n=length(hetero), pch=4, col="red", tag=data$sample.id, sizeofpoint=2)