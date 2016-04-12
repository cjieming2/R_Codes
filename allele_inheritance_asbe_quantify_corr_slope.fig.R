setwd("C:/Users/Jieming/Documents/thesis/mark_work/allele_specificity/inheritance")
# setwd("C:/Users/JM/thesis/mark_work/allele_specificity/compare trios ratios inheritance/test")

## parameters
qnumDataPointsCutOff = 10;
colratio1=6
colratio2=8

title = 'combined'
ctcf.na12878.91.file <- "inherit.combined_CTCF_NA12878_intHets.bed.alleleBed.inherit.combined_CTCF_NA12891_intHets.bed.alleleBed.asb.bed"
ctcf.na12878.92.file <- "inherit.combined_CTCF_NA12878_intHets.bed.alleleBed.inherit.combined_CTCF_NA12892_intHets.bed.alleleBed.asb.bed"
ctcf.na12891.92.file <- "inherit.combined_CTCF_NA12891_intHets.bed.alleleBed.inherit.combined_CTCF_NA12892_intHets.bed.alleleBed.asb.bed";
pu1.na12878.91.file <- "inherit.combined_PU1_NA12878_intHets.bed.alleleBed.inherit.combined_PU1_NA12891_intHets.bed.alleleBed.asb.bed"
pu1.na12878.92.file <- "inherit.combined_PU1_NA12878_intHets.bed.alleleBed.inherit.combined_PU1_NA12892_intHets.bed.alleleBed.asb.bed"
pu1.na12891.92.file <- "inherit.combined_PU1_NA12891_intHets.bed.alleleBed.inherit.combined_PU1_NA12892_intHets.bed.alleleBed.asb.bed";
ase.na12878.91.file <- "inherit.rnaseq.NA12878.old.intHets.bed.alleleBed.inherit.rnaseq.NA12891.old.intHets.bed.alleleBed.asb.bed";
ase.na12878.92.file <- "inherit.rnaseq.NA12878.old.intHets.bed.alleleBed.inherit.rnaseq.NA12892.intHets.bed.alleleBed.asb.bed";
ase.na12891.92.file <- "inherit.rnaseq.NA12891.old.intHets.bed.alleleBed.inherit.rnaseq.NA12892.intHets.bed.alleleBed.asb.bed";


## hets snp.calls chr-pos chr pos refallele mat pat child phase
ctcf.na12878.91.data <- read.table(ctcf.na12878.91.file, header=F, sep = "\t", stringsAsFactors = F)
ctcf.na12878.92.data <- read.table(ctcf.na12878.92.file, header=F, sep = "\t", stringsAsFactors = F)
pu1.na12878.91.data <- read.table(pu1.na12878.91.file, header=F, sep = "\t", stringsAsFactors = F)
pu1.na12878.92.data <- read.table(pu1.na12878.92.file, header=F, sep = "\t", stringsAsFactors = F)
ase.na12878.91.data <- read.table(ase.na12878.91.file, header=F, sep = "\t", stringsAsFactors = F)
ase.na12878.92.data <- read.table(ase.na12878.92.file, header=F, sep = "\t", stringsAsFactors = F)

## calculate correlation
ctcf.na12878.91.cor = cor(ctcf.na12878.91.data[,colratio1],ctcf.na12878.91.data[,colratio2],method="spearman");
ctcf.na12878.92.cor = cor(ctcf.na12878.92.data[,colratio1],ctcf.na12878.92.data[,colratio2],method="spearman");
pu1.na12878.91.cor = cor(pu1.na12878.91.data[,colratio1],pu1.na12878.91.data[,colratio2],method="spearman");
pu1.na12878.92.cor = cor(pu1.na12878.92.data[,colratio1],pu1.na12878.92.data[,colratio2],method="spearman");
ase.na12878.91.cor = cor(ase.na12878.91.data[,colratio1],ase.na12878.91.data[,colratio2],method="spearman");
ase.na12878.92.cor = cor(ase.na12878.92.data[,colratio1],ase.na12878.92.data[,colratio2],method="spearman");

p2c.r = c(ctcf.na12878.91.cor, ctcf.na12878.92.cor,
          pu1.na12878.91.cor, pu1.na12878.92.cor,
          ase.na12878.91.cor, ase.na12878.92.cor)

## plotting
x11()
par(mar=c(5,10,4,6),xpd=TRUE)
barplot(p2c.r, beside=T, col=c("blue","red"), cex.axis=2, 
        cex.names=1.5, cex.lab=2, xlab="TFs", ylab="correlation",
        names.arg=c("CTCF","CTCF","PU.1","PU.1","ASE","ASE"), ylim=c(0,1))
legend(15,1,c('Child v Father','Child v Mother'), 
              col=c('blue','red'), 
              text.col = "black", bg = 'white',cex=2, pch=15)

###########################################################################
## calculate slope
ctcf.na12891.92.data <- read.table(ctcf.na12891.92.file, header=F, sep = "\t", stringsAsFactors = F)
pu1.na12891.92.data <- read.table(pu1.na12891.92.file, header=F, sep = "\t", stringsAsFactors = F)
ase.na12891.92.data <- read.table(ase.na12891.92.file, header=F, sep = "\t", stringsAsFactors = F)

ctcf.na12878.91.slope = lm(ctcf.na12878.91.data[,colratio1]~ctcf.na12878.91.data[,colratio2])
ctcf.na12878.92.slope = lm(ctcf.na12878.92.data[,colratio1]~ctcf.na12878.92.data[,colratio2])
ctcf.na12891.92.slope = lm(ctcf.na12891.92.data[,colratio1]~ctcf.na12891.92.data[,colratio2])
pu1.na12878.91.slope = lm(pu1.na12878.91.data[,colratio1]~pu1.na12878.91.data[,colratio2])
pu1.na12878.92.slope = lm(pu1.na12878.92.data[,colratio1]~pu1.na12878.92.data[,colratio2])
pu1.na12891.92.slope = lm(pu1.na12891.92.data[,colratio1]~pu1.na12891.92.data[,colratio2])
ase.na12878.91.slope = lm(ase.na12878.91.data[,colratio1]~ase.na12878.91.data[,colratio2])
ase.na12878.92.slope = lm(ase.na12878.92.data[,colratio1]~ase.na12878.92.data[,colratio2])
ase.na12891.92.slope = lm(ase.na12891.92.data[,colratio1]~ase.na12891.92.data[,colratio2])
p2c.slope = c(round(summary(ctcf.na12878.91.slope)$coefficients[[2]],2),round(summary(ctcf.na12878.92.slope)$coefficients[[2]],2),
              round(summary(ctcf.na12891.92.slope)$coefficients[[2]],2),
              round(summary(pu1.na12878.91.slope)$coefficients[[2]],2),round(summary(pu1.na12878.92.slope)$coefficients[[2]],2),
              round(summary(pu1.na12891.92.slope)$coefficients[[2]],2),
              round(summary(ase.na12878.91.slope)$coefficients[[2]],2),round(summary(ase.na12878.92.slope)$coefficients[[2]],2),
              round(summary(ase.na12891.92.slope)$coefficients[[2]],2))