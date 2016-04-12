# This software is supplied without any warranty or guaranteed support whatsoever.
# NUS CME can not be responsible for its use, misuse, or functionality.
#
# This is example code for producing the varLD region plot of VKORC1 that is similar to 
# Figure 1 in the Bioinformatics Application Note. 
# Anybody who is moderately proficient in R should be able to modify the codes here
# to customise the plot required. 
#
# Author: rick and yy teo
# Version: 1.0
# Last Updated : 07 March 2010
######################################################################################
## global variables to modify for own use 

# specifies the folder where the output files from varLD is stored
PATH = "C:/Documents and Settings/chenjm/Desktop/psoriasis varld/"

# filename prefix before the chromosome number, i.e. CHS_INS_chr1, CHS_INS_chr2
FILENAME = "varld-3291chinese-4524europeans-genomewide-chr" 

# chromosome of region to plot
chr = 1

# percentile of the genomewide distribution to highlight 
percentile.out 	= c(0.95, 0.99, 0.999, 0.9999)

# varLD score corresponding to the stated percentiles in percentile.out ### this is taken from the output when you run the standar
varLD.threshold 	= c(1.94, 3.38, 5.40, 7.13)

# start and end base-pair coordinates of the region to plot
region.start = 150000000
region.end = 152000000


temp.in <- read.table(paste(PATH, FILENAME, chr, "_standardized.out", sep=""), sep="\t", header = T)
region.flag <- which(temp.in[,"position"] >= region.start & temp.in[,"position"] <= region.end)

x11()
plot(0, 1, xlab = "Physical position (Mb)", ylab = "Standardized score", las = 1, type = "n", xlim = c(region.start, region.end)/10^6, ylim = c(min(varLD.threshold), max(temp.in[region.flag,"standardized_score"])), bty = "n", cex.lab = 1.4, cex.main = 1.4)
points(temp.in[region.flag, "position"]/10^6, temp.in[region.flag, "standardized_score"], pch = 16, col = "red")

abline(v=150.7, lty=2, col="red")
abline(v=151.1, lty=2, col="red")

n.length.threshold <- length(varLD.threshold)
if (n.length.threshold >= 2){
   for (i in 2:n.length.threshold){
      abline(h = varLD.threshold[i], lty = 2)
      text(region.start/10^6 + 0.2, varLD.threshold[i] + 0.15, labels = paste("Top ", round((1 - percentile.out[i]) * 100, 2), "%", sep=""))
   }
}