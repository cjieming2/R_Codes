# setwd("C:/Users/JM/thesis/mark_work/allele_specificity/compare trios ratios inheritance")
# setwd("C:/Users/JM/thesis/mark_work/allele_specificity/compare trios ratios inheritance/ase/allHets")
setwd("C:/Users/Jieming/Documents/thesis/mark_work/allele_specificity/inheritance")

## parameters
numDataPointsCutOff = 5;
title = '3_2Hets'

# TF = 'CTCF'
# na12878.91.file <- "inherit.combined_CTCF_NA12878_intHets.bed.alleleBed.inherit.combined_CTCF_NA12891_intHets.bed.alleleBed.asb.bed"
# na12878.92.file <- "inherit.combined_CTCF_NA12878_intHets.bed.alleleBed.inherit.combined_CTCF_NA12892_intHets.bed.alleleBed.asb.bed"
# na12891.92.file <- "inherit.combined_CTCF_NA12891_intHets.bed.alleleBed.inherit.combined_CTCF_NA12892_intHets.bed.alleleBed.asb.bed"
# TF = 'POL2'
# na12878.91.file <- "inherit.combined_Pol2_NA12878_intHets.bed.alleleBed.inherit.combined_POL2_NA12891_intHets.bed.alleleBed.asb.bed";
# na12878.92.file <- "inherit.combined_Pol2_NA12878_intHets.bed.alleleBed.inherit.combined_POL2_NA12892_intHets.bed.alleleBed.asb.bed";
# na12891.92.file <- "inherit.combined_POL2_NA12891_intHets.bed.alleleBed.inherit.combined_POL2_NA12892_intHets.bed.alleleBed.asb.bed";
# TF = 'PU1'
# na12878.91.file <- "inherit.combined_PU1_NA12878_intHets.bed.alleleBed.inherit.combined_PU1_NA12891_intHets.bed.alleleBed.asb.bed"
# na12878.92.file <- "inherit.combined_PU1_NA12878_intHets.bed.alleleBed.inherit.combined_PU1_NA12892_intHets.bed.alleleBed.asb.bed"
# na12891.92.file <- "inherit.combined_PU1_NA12891_intHets.bed.alleleBed.inherit.combined_PU1_NA12892_intHets.bed.alleleBed.asb.bed"
TF = 'rnaseq'
na12878.91.file <- "inherit.rnaseq.NA12878.old.intHets.bed.alleleBed.inherit.rnaseq.NA12891.old.intHets.bed.alleleBed.asb.bed"
na12878.92.file <- "inherit.rnaseq.NA12878.old.intHets.bed.alleleBed.inherit.rnaseq.NA12892.intHets.bed.alleleBed.asb.bed"
na12891.92.file <- "inherit.rnaseq.NA12891.old.intHets.bed.alleleBed.inherit.rnaseq.NA12892.intHets.bed.alleleBed.asb.bed"

## hets snp.calls chr-pos chr pos refallele mat pat child phase
na12878.91.data <- read.table(na12878.91.file, header=F, sep = "\t", stringsAsFactors = F)
na12878.92.data <- read.table(na12878.92.file, header=F, sep = "\t", stringsAsFactors = F)
na12891.92.data <- read.table(na12891.92.file, header=F, sep = "\t", stringsAsFactors = F)

## plotting e.g. 78.91for inheritance with genotype of third individual
## this function takes the input and outputs a figure
## 1) calc counts,correl,cortest of the input data an
## 2) 
plot_inherit <- function(data,colRatio1,colRatio2,name1,name2,colNumIndiv,colThirdIndiv)
{
  data.ratio <- cbind(data[,colRatio1], data[,colRatio2])
  name = paste(name1,"_",name2,sep='')
  
  counts   = nrow(data);
  correl   = cor(data.ratio[,1],data.ratio[,2],method="spearman");
  cortest  = cor.test(data.ratio[,1],data.ratio[,2], alternative="two.sided", method="spearman")
  
  pdf(paste(title,"_",name,".pdf",sep=''));
  plot(data.ratio, main=name, 
       xlab=paste('ref/total for ',name1),
       ylab=paste('ref/total for ',name2),
       ylim=c(0,1),xlim=c(0,1));
  text(0.2,0.45,paste("p val=",signif(cortest$p.value,3)))
  text(0.2,0.5,paste("Spearman's cor=",signif(correl,3)))
  text(0.2,0.55,paste("numPoints=",counts))
  
  
  ##### color by homo 00 and 11 #####
  ## hets snp.calls chr-pos chr pos refallele mat pat child phase
  my.00.data <- data[(data[,colNumIndiv] == 2) & (data[,colThirdIndiv] == 0),]
  my.11.data <- data[(data[,colNumIndiv] == 2) & (data[,colThirdIndiv] == 11),]
  
  my.00.data.ratio <- cbind(my.00.data[,colRatio1], my.00.data[,colRatio2])
  my.11.data.ratio <- cbind(my.11.data[,colRatio1], my.11.data[,colRatio2])
  
  points(my.00.data.ratio,col='red');
  points(my.11.data.ratio,col='blue');
  
  # text(0.3,0.75,paste("p val=",signif(cortest$p.value,3)))
  # text(0.3,0.8,paste("Spearman's cor=",signif(correl,3)))
  # text(0.3,0.85,paste("numPoints=",counts))
  dev.off();
}

plot_inherit(na12878.91.data,6,8,paste("NA12878",TF),paste("NA12891",TF),5,7)
plot_inherit(na12878.92.data,6,8,paste("NA12878",TF),paste("NA12892",TF),5,7)
plot_inherit(na12891.92.data,6,8,paste("NA12891",TF),paste("NA12892",TF),5,7)
