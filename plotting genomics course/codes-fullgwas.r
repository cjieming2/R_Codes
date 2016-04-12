## ========================================================================== ##
## Genetics and Genomics of Complex Phenotypes April 2010
## Day5_Session2 (Full GWAS Practical)
## ========================================================================== ##

# 1. Reading in the data and plotting missingness versus heterozygosity #
het <- read.table("D:/GGCP_April2010/Practicals/Day5_Session2/full_gwas.het", header=T, stringsAsFactors=F)
missing <- read.table("D:/GGCP_April2010/Practicals/Day5_Session2/full_gwas.imiss", header=T, stringsAsFactors=F)
dat.merge <- merge(het, missing, by=c("FID","IID"), sort=F)
dat.merge$heterozygosity <- 1-dat.merge$O.HOM./dat.merge$N.NM.
dat.merge$missingness <- dat.merge$N_MISS/dat.merge$N_GENO
plot(dat.merge[,"heterozygosity"]*100, dat.merge[,"missingness"]*100,
     xlab="Heterozygosity (%)", ylab="Sample Missingness (%)", pch=20)
abline(h=5, lty=2)
abline(v=25, lty=2)
abline(v=34, lty=2)
drop.id <- which((dat.merge$missingness>0.05) | (dat.merge$heterozygosity<0.25 | dat.merge$heterozygosity>0.34))
points(dat.merge[drop.id,"heterozygosity"]*100, dat.merge[drop.id,"missingness"]*100,, col="red", pch=20)
length(drop.id)

# 2. Population Structure #
# With Hapmap #
hapmap <- read.table("D:/GGCP_April2010/Practicals/Day5_Session2/196hapmap-558case-control-pairs-obesecase-26PCA-104545autosomalsnps.pca",
                     header=T, sep="\t", stringsAsFactors=F)
dat.plot <- hapmap
dat.plot$pop.id <- 0
dat.plot[grep("CEU",dat.plot$sample.id),"pop.id"] <- 6
dat.plot[grep("JPT",dat.plot$sample.id),"pop.id"] <- 5
dat.plot[grep("CHB",dat.plot$sample.id),"pop.id"] <- 7
dat.plot[grep("YRI",dat.plot$sample.id),"pop.id"] <- 9
dat.plot[grep("SP",dat.plot$sample.id),"pop.id"] <- 2
dat.plot[which(dat.plot$pop.id==0),"pop.id"] <- 3
#  1st with the 2nd component  #
plot(dat.plot[,"PC1"], dat.plot[,"PC2"], xlab="first principal component",
     ylab="second principal component", pch=20, col=dat.plot$pop.id, main="")
legend(locator(1), pch=20, col=c(6,7,5,9,2,3),
      c("CEU","CHB","JPT","YRI","controls","cases"), bty="n", cex=0.8)
identify(dat.plot[,"PC1"], dat.plot[,"PC2"], labels = as.character(dat.plot[,"sample.id"]))
points(dat.plot[which((dat.plot[,"pop.id"]==2 | dat.plot[,"pop.id"]==3) & dat.plot[,"PC1"]>0.01), "PC1"],
       dat.plot[which((dat.plot[,"pop.id"]==2 | dat.plot[,"pop.id"]==3) & dat.plot[,"PC1"]>0.01), "PC2"], pch=21, cex=1.2)
hapmap.outlier <- dat.plot[which((dat.plot[,"pop.id"]==2 | dat.plot[,"pop.id"]==3) & dat.plot[,"PC1"]>0.01), "sample.id"]

# With SGVP #
sgvp <- read.table("D:/GGCP_April2010/Practicals/Day5_Session2/268sgvp-558case-control-pairs-obesecase-26PCA-104545autosomalsnps.pca",
                     header=T, sep="\t", stringsAsFactors=F)
samples <- read.table("D:/GGCP_April2010/Practicals/Day5_Session2/samples-information.txt",
                      header=T, sep="\t", stringsAsFactors=F)

info <- samples[,c("IID","population")]
names(info) <- c("sample.id","pop.id")

dat.plot <- merge(sgvp, info, by="sample.id", sort=F, all.x=T)
dat.plot[which(dat.plot[,"pop.id"]=="CHS"), "pop.id"] <- 2
dat.plot[which(dat.plot[,"pop.id"]=="MAS"), "pop.id"] <- 4
dat.plot[which(dat.plot[,"pop.id"]=="INS"), "pop.id"] <- 3
dat.plot[grep("SP",dat.plot$sample.id),"pop.id"] <- 5
dat.plot[which(is.na(dat.plot$pop.id)),"pop.id"] <- 6
#  1st with the 2nd component  #
plot(dat.plot[,"PC1"], dat.plot[,"PC2"], xlab="first principal component",
     ylab="second principal component", pch=20, col=dat.plot$pop.id, main="")
points(dat.plot[which(dat.plot[,"pop.id"]==2),"PC1"], dat.plot[which(dat.plot[,"pop.id"]==2),"PC2"],
       col=2, pch=20)
legend(locator(1), pch=20,col=c(2,4,3,5,6), c("CHS", "MAS", "INS","controls","cases"), bty="n", cex=0.8)
identify(dat.plot[,"PC1"], dat.plot[,"PC2"], labels = as.character(dat.plot[,"sample.id"]))
# outliers
points(dat.plot[which((dat.plot[,"pop.id"] > 4) & (dat.plot[,"PC1"] > 0.06)), "PC1"],
       dat.plot[which((dat.plot[,"pop.id"] > 4) & (dat.plot[,"PC1"] > 0.06)), "PC2"], pch=21, cex=1.2)
points(dat.plot[which((dat.plot[,"pop.id"] > 4) & (dat.plot[,"PC2"] < -0.040)), "PC1"],
       dat.plot[which((dat.plot[,"pop.id"] > 4) & (dat.plot[,"PC2"] < -0.040)), "PC2"], pch=21, cex=1.2)

                                                      
# 3. Association testing and plotting #
# PP-plot #
dat.plot <- read.table("D:/GGCP_April2010/Practicals/Day5_Session2/full_gwas_qc.assoc.logistic", header=T,stringsAsFactors=F)
dat.snp.plot <- dat.plot[grep("ADD", dat.plot[,"TEST"]),]
size <- nrow(dat.snp.plot)
p.exp <- (-log10(qunif(seq(from=1, to=size, by=1)/(size+1))))
p.obs <- (-log10(sort(dat.snp.plot[,"P"])))
plot(p.exp, p.obs, xlab="-log(expected P values)", pch=20,
     ylab="-log(observed Pvalues)", main=paste("Diabetes case control ",nrow(dat.snp.plot), " SNPs"))
abline(a=0, b=1)
# Manhattan-plot #
chr.no <- unique(dat.snp.plot$CHR)
dist.new <- NULL
dist.new <- dat.snp.plot[which(dat.snp.plot$CHR==chr.no[1]), "BP"]*10^(-10)
for (i in 2:length(chr.no))
    {
    dist.chromosome.tmp <- dat.snp.plot[which(dat.snp.plot$CHR==chr.no[i]), "BP"]*10^(-10) + max(dist.new)
    dist.new <- append(dist.new, dist.chromosome.tmp)
    }
chr.mean <- NULL
for (i in 1:22)
    {
    chromosome <- which(dat.snp.plot$CHR==i)
    chromosome.mean.tmp <- mean(dist.new[chromosome])
    chr.mean <- rbind(chr.mean, chromosome.mean.tmp)
    }
plot(dist.new, (-log10(dat.snp.plot[which(dat.snp.plot$CHR>0),"P"])),ylim=c(0,10), type="p", pch=20,
     col=dat.snp.plot$CHR, axes=F, xlab="", ylab="",
     main=paste("Diabetes case control ",nrow(dat.snp.plot), " SNPs"))
axis(side=1, at=chr.mean, labels=seq(from=1, to=22, by=1), tick=F, cex.lab=0.8, cex.axis=0.7)
axis(side=2, cex.lab=0.8,cex.axis=0.7)
mtext("chromosome",side=1,line=2,cex=0.8)
mtext("-log10(pvalues)",side=2,line=2,cex=0.8)
box(which="plot")
abline(h= -log10(0.00000005))


# 4. Extracting candidate region for visualisation with imputed data #
dat <- read.table("D:/GGCP_April2010/Practicals/Day5_Session2/candidate-region-plotting.txt",
                  sep="\t", header=T, stringsAsFactors=F)
snp = "rs8050136"
locusname = "FTO"
chr = "16"
locus = dat
range = 9
best.pval = locus[which(locus$SNP==snp),"PVAL"]

hit <- locus[which(locus$SNP==snp),]

# size of the region
min.pos <- min(locus$POS) - 10000
max.pos <- max(locus$POS) + 10000
size.pos <- max.pos - min.pos
center.pos <- min.pos + ( size.pos / 2 )
center.100kb.pos <- round(center.pos / 100000) * 100000
offset.100kb.pos <- round((size.pos/3) / 100000) * 100000
# range of y-axis
# this dedicates 33% of the yaxis to the genes, labels, recomb rate
offset <- ( range * 4 / 3 ) - range
big.range <- range + offset 
ystart.gene <- - offset
ystart.recomb <- - offset + (big.range / 8)
# recombination rate 
recomb <- read.table(paste("D:/GGCP_April2010/Practicals/Day5_Session2/genetic_map_chr", chr, "_b36.txt", sep=""), header=T)
keep.recomb <- subset(recomb, recomb[,1] > min.pos & recomb[,1] < max.pos)
# genes in the region
genelist <- read.table(paste("D:/GGCP_April2010/Practicals/Day5_Session2/known_genes_build36_240310_chr", chr, ".txt", sep=""), header=T)
genes.in.locus <- subset(genelist, ( genelist$START > min.pos & genelist$START < max.pos ) | ( genelist$STOP > min.pos & genelist$STOP < max.pos) )
print(genes.in.locus)

# start plot with recombination rate (in background)
par(mar=c(4,4,3,4))
plot(keep.recomb[,1], ystart.recomb + ( ( keep.recomb[,2] / 60 ) * ( 6 * big.range / 8 )), type="l", col="lightblue", lwd=1, xlim=c(min.pos, max.pos), ylim=c(-offset,range), xlab="", ylab="", main=locusname, axes=F)
# axes, titles and legends
mtext(paste("Chromosome", chr, "position (kb)", sep=" "), side=1, line=2.5)
axis(1, at=c(center.100kb.pos - offset.100kb.pos, center.100kb.pos, center.100kb.pos + offset.100kb.pos), labels=c((center.100kb.pos - offset.100kb.pos) / 1000, center.100kb.pos / 1000, (center.100kb.pos + offset.100kb.pos) / 1000), las=1) 
axis(2, at=seq(0,range,2), labels=seq(0,range,2), las=1) 
mtext("Observed (-logP)", side=2, at=(range/2), line=2)
axis(4, at=c( ystart.recomb, ystart.recomb + (big.range / 4), ystart.recomb + ( 2 * big.range / 4), ystart.recomb + ( 3 * big.range / 4 ) ), labels=c("0","20","40","60"), las=1)
mtext("Recombination rate (cM/Mb)", side=4, at=(-offset+big.range/2), line=2)
box()
lines(c(min.pos, max.pos), c(0,0), lty="dotted", lwd=1, col="black")
# this is the hit
points(hit$POS, -(log10(hit$PVAL)), pch=23, cex=2.5, bg="red")
if ( -(log10(best.pval)) < range ) {
	points(hit$POS, -(log10(best.pval)), pch=23, cex=2.5, bg="blue")
	text(hit$POS, -(log10(best.pval)), labels=c(paste("P=",best.pval,sep="")), pos=4, offset=2)
} else {
	points(hit$POS, range, pch=23, cex=2.5, bg="blue")
	text(hit$POS, range, labels=c(paste("P=",best.pval,sep="")), pos=4, offset=1)
}

# genotyped markers
markers.typed <- subset(locus, (row.names(locus) != snp & locus$TYPE == "typed"))
# imputed SNPs
markers.imputed <- subset(locus, (row.names(locus) != snp & locus$TYPE == "imputed"))
# plot the imputed SNPs
points(markers.imputed$POS, -(log10(markers.imputed$PVAL)), pch=20, cex=0.8, col="red")
# plot the genotyped markers
points(markers.typed$POS, -(log10(markers.typed$PVAL)), pch=20, cex=0.8, col="black")
# plot the genes
genes.tmp <- apply(genes.in.locus, 1, function(x) {paste(x, collapse="")})
unique.genes <- which(duplicated(genes.tmp)==T)
genes.in.locus.new <- genes.in.locus[-unique.genes, ]
if (locusname %in% genes.in.locus.new[,"GENE"]) {
  genes.in.locus.snp <- genes.in.locus.new[which(genes.in.locus.new[,"GENE"]%in%locusname),]
  for ( i in 1:nrow(genes.in.locus.snp) ) { 
  	if ( genes.in.locus.snp[i,]$STRAND == "+" ) {
  		arrows(max(genes.in.locus.snp[i,]$START, min.pos), -offset+i/nrow(genes.in.locus.snp), min(genes.in.locus.snp[i,]$STOP, max.pos), -offset+i/nrow(genes.in.locus.snp), length=0.05, lwd=2, code=2, lty="solid", col="darkgreen")
  	} else {		
  		arrows(max(genes.in.locus.snp[i,]$START, min.pos), -offset+i/nrow(genes.in.locus.snp), min(genes.in.locus.snp[i,]$STOP, max.pos), -offset+i/nrow(genes.in.locus.snp), length=0.05, lwd=2, code=1, lty="solid", col="darkgreen")
  	}
  	if ( ! is.na(genes.in.locus.snp[i,]$GENE) ) {
  		text(genes.in.locus.snp[i,]$START + (genes.in.locus.snp[i,]$SIZE / 2), -offset + i/nrow(genes.in.locus.snp) + ( big.range / 20 ), labels=genes.in.locus.snp[i,]$GENE, cex=0.8)
  	}
  }
}
