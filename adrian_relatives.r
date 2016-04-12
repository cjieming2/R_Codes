setwd("C:/Documents and Settings/tana2/Desktop/tana2/PASNP/v3/1928samples-54794snps-pasnp-hapmap/20080507-relatives")

samples <- read.table("pasnp.sar", header=T, sep = "\t")

col = c("yellow", "red", "black", "green")

populations = c("TH-MA", "CN-HM", "CEU-NA", "ID-RA", "IN-NL", "IN-SP")
populations = unique(samples$population.id)

for (population in populations)
{
population = populations[population]
#x11(10,7, pointsize = 18)
png(1000, 700, pointsize = 18, file = paste(population, ".png", sep = ""))
layout(matrix(c(rep(1,37), rep(2,13)), 1, 50, byrow = TRUE))
layout.show(2)
plot(samples$ibs.mean, samples$ibs.stdev, xlab="IBS mean", ylab="IBS standard deviation", type = "n", main=paste(population," Relative Plot"), pch = 21,  bg = col[unclass(samples$first.degree.relationship)])

sub.samples = subset(samples, first.degree.relationship != "n/a")
points(sub.samples$ibs.mean, sub.samples$ibs.stdev, pch = 4, col = col[unclass(sub.samples$first.degree.relationship)])
sub.samples = subset(samples, population.id == population)

points(sub.samples$ibs.mean, sub.samples$ibs.stdev, pch = 21, bg = col[unclass(sub.samples$first.degree.relationship)])
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, unique(samples$first.degree.relationship), pt.bg = col[unclass(unique(samples$first.degree.relationship))],
        text.col = "black", pch = 21,
        bg = 'white')

dev.off()
}

col = c("yellow", "red", "black", "green")
x11(10,7, pointsize = 18)
layout(matrix(c(rep(1,37), rep(2,13)), 1, 50, byrow = TRUE))
layout.show(2)
plot(samples$ibs.mean, samples$ibs.stdev, xlab="IBS mean", ylab="IBS standard deviation", type = "n", main="PASNP Relative Plot", pch = 21,  bg = col[unclass(samples$first.degree.relationship)])
sub.samples = subset(samples, T)
points(sub.samples$ibs.mean, sub.samples$ibs.stdev, pch = 21, bg = col[unclass(sub.samples$first.degree.relationship)])
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(samples$first.degree.relationship)), pt.bg = col[unclass(sort(unique(samples$first.degree.relationship)))],
        text.col = "black", pch = 21,
        bg = 'white')

col = c("yellow", "red", "black", "green")
x11(10,7, pointsize = 18)
layout(matrix(c(rep(1,37), rep(2,13)), 1, 50, byrow = TRUE))
layout.show(2)
plot(samples$ibs.mean, samples$ibs.stdev, xlab="IBS mean", ylab="IBS standard deviation", type = "n", main="PASNP (Without first degree relatives) Relative Plot", pch = 21,  bg = col[unclass(samples$first.degree.relationship)])
sub.samples = subset(samples, P1==1)
#sub.samples = subset(samples, T)
points(sub.samples$ibs.mean, sub.samples$ibs.stdev, pch = 21, bg = col[unclass(sub.samples$first.degree.relationship)])
par(mar = c(0,0,0,0))
plot(0, 0, axes = F, xlab= NA, ylab= NA, type = "n")
legend(-1, 0.8, sort(unique(samples$first.degree.relationship)), pt.bg = col[unclass(sort(unique(samples$first.degree.relationship)))],
        text.col = "black", pch = 21,
        bg = 'white')
