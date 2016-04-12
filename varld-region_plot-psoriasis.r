### specifies the folder where the output files from varLD is stored
PATH = "C:/Documents and Settings/chenjm/Desktop/work_documentation/JJ/psoriasis-finemapping/varld-LCE/psoriasis varld output+standardization/"

### filename prefix before the chromosome number, i.e. CHS_INS_chr1, CHS_INS_chr2
FILENAME = "varld-3291chinese-4524europeans-genomewide-chr" 

### chromosome of region to plot
chr = 1

### percentile of the genomewide distribution to highlight 
percentile.out 	= c(0.95, 0.99, 0.999, 0.9999)

### varLD score corresponding to the stated percentiles in percentile.out ### this is taken from the output when you run the standar
varLD.threshold 	= c(1.94, 3.38, 5.40, 7.13)

### start and end base-pair coordinates of the region to plot
region.start = 149000000
region.end = 153000000

temp.in <- read.table(paste(PATH, FILENAME, chr, "_standardized.out", sep=""), sep="\t", header = T)
region.flag <- which(temp.in[,"position"] >= region.start & temp.in[,"position"] <= region.end)

x11()
plot(temp.in[region.flag, "position"]/10^6, temp.in[region.flag, "standardized_score"], 
	xlab = "Mean position of sliding window (Mb)", ylab = "Standardized score", 
	las = 1, type = "p", 
	xlim = c(region.start, region.end)/10^6, ylim = c(-2,6), 
	bty = "l", cex.lab = 1.4, cex.main = 1.4,
	main = "LCE gene cluster",
	col = 'black')

rect(150.7, -2.33, 151.1, 6, col="gray85", border = "transparent")
points(temp.in[region.flag, "position"]/10^6, temp.in[region.flag, "standardized_score"], pch = 16, col = "black")
#lines(temp.in[region.flag, "position"]/10^6, temp.in[region.flag, "standardized_score"], pch = 16, col = "black")

usr <- par("usr")
#abline(v=150.7, lty=2, col="black")
#abline(v=151.1, lty=2, col="black")

abline(h = 1.94, lty = 2)
text(region.start/10^6 + 0.2, 1.94 + 0.15, labels = paste("Top ", round((1 - 0.95) * 100, 2), "%", sep=""))
abline(h = 3.38, lty = 2)
text(region.start/10^6 + 0.2, 3.38 + 0.15, labels = paste("Top ", round((1 - 0.99) * 100, 2), "%", sep=""))
abline(h = 5.4, lty = 2)
text(region.start/10^6 + 0.2, 5.4 + 0.15, labels = paste("Top ", round((1 - 0.999) * 100, 2), "%", sep=""))

top5 <- which(temp.in[,"standardized_score"] >= 1.94)
points(temp.in[top5, "position"]/10^6, temp.in[top5, "standardized_score"], pch = 16, col = "green")
#clip(usr[1], usr[2], 1.94, 3.38)
#lines(temp.in[region.flag, "position"]/10^6, temp.in[region.flag, "standardized_score"], col = 'green')

top1 <- which(temp.in[,"standardized_score"] >= 3.38)
points(temp.in[top1, "position"]/10^6, temp.in[top1, "standardized_score"], pch = 16, col = "red")
#clip(usr[1], usr[2], 3.38, 5.40)
#lines(temp.in[region.flag, "position"]/10^6, temp.in[region.flag, "standardized_score"], col = 'red')


