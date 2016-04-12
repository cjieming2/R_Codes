setwd("C:/Documents and Settings/chenjm/Desktop")

filename="gafr-fst-aafd.txt"
data <- read.table(filename, header=TRUE, sep="\t")

#col is border, bg is background
colors <- c("black", "orange", "blue", "red", "green", "magenta", "lightblue1")

#zoomed plot
x11()
plot(data$Fst, data$GAFR, xlab="Fst", ylab="GAFR", main="", pch=20, col="orange")
#legend(-1, 0.8, "GAFR<2 Fst<0.1", pt.bg = "blue", text.col = "black", pch = 21, bg = 'white')

# add red line
abline(v=0.1, col="red")

# identify desired points
# this stores the ID of the points in vector/matrix id
id <- which (data$GAFR < 2 & data$Fst < 0.1)
id1 <- which (data$AAFD < 0.05 & data$Fst < 0.1)

# highlight points in red
points(data$Fst[id],data$GAFR[id], pch=20, col="blue")
points(data$Fst[id1],data$GAFR[id1], pch=20, col="green")


# add label sample-id on the points
# pos 1-bottom 2-left 3-top 4-right
# offset distance away from point
#text(hetero[id], cr[id], labels=data$sample.id[id], pos=2, offset=0.5, col="red")

#scale of 0-100%
#x11()
#plot(hetero, cr, xlab="heterozygosity (%)", ylab="sample-call-rate (%)", main="samplecr VS heterozygosity plot for 43 samples on 909622 SNPs", pch=21, col="black", xlim=c(0,100), ylim=c(0,100))