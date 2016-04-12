# This software is supplied without any warranty or guaranteed support whatsoever.
# NUS CME can not be responsible for its use, misuse, or functionality.
#
# This R script can be used for standardizing the genome-wide scores to have 
# a mean of 0 and a variance of 1, and to output files with the standardized scores. 
# It assumes that varld scores are available for all 22 autosomal chromosomes. 
# This code also identifies the varLD scores that correspond to the stated percentiles.
#
# Author: rick and yy teo
# Version: 1.0
# Last Updated : 07 March 2010
######################################################################################
## global variables to modify for own use 
# specifies the folder containing all varld output for 22 autosomal chromosomes
PATH = "C:/Documents and Settings/chenjm/Desktop/work_documentation/JJ/psoriasis-finemapping/varld-LCE/psoriasis varld output+standardization/"
# filename prefix before the chromosome number, i.e. CHS_INS_chr1, CHS_INS_chr2 
FILENAME = "varld-3291chinese-4524europeans-genomewide-chr" 
# percentile of the genomewide distribution to highlight, default to 95%, 99%, 99.9% and 99.99%.
percentile.out = c(0.95, 0.99, 0.999, 0.9999)


varLD.out <- {}
chr.store <- {}
for (chr in 1:22){
   varLD.temp <- read.table(paste(PATH, FILENAME, chr, ".txt", sep=""), sep="\t", header = T)
   varLD.out <- rbind(varLD.out, varLD.temp)
   chr.store <- c(chr.store, rep(chr, dim(varLD.temp)[1]))
   print(paste("completed reading in unstandardized varLD output file for chromosome ", chr, sep=""))
}
varLD.mean <- mean(varLD.out[,"raw_score"])
varLD.sd <- sd(varLD.out[,"raw_score"])
standardized_score <- (varLD.out[,"raw_score"] - varLD.mean)/varLD.sd
varLD.out <- cbind(varLD.out, standardized_score)
varLD.threshold <- quantile(standardized_score, probs = percentile.out)

for (chr in 1:22){
   chr.flag <- which(chr.store == chr)
   write.table(varLD.out[chr.flag,], paste(PATH, FILENAME, chr, "_standardized.out", sep=""), sep="\t", quote=F, row.names=F)
   print(paste("completed writing out standardized varLD output file for chromosome ", chr, sep=""))
}
n.length.percentile <- length(percentile.out)
for (i in 1:n.length.percentile){
   print(paste("varLD threshold for ", percentile.out[i], " = ", varLD.threshold[i], sep=""))
} 
