# setwd("C:/Users/JM/thesis/mark_work/allele_specificity/compare trio venns total")
setwd("C:/Users/JM/thesis/mark_work/allele_specificity/compare trio venns asb ase")

## filenames
# title   = "trio_by_total"
# na12878.file <- "NA12878.total.snps.gtypes";
# na12891.file <- "NA12891.total.snps.gtypes";
# na12892.file <- "NA12892.total.snps.gtypes";

# title   = "trio_by_hets"
# na12878.file <- "NA12878.het.snps.gtypes";
# na12891.file <- "NA12891.het.snps.gtypes";
# na12892.file <- "NA12892.het.snps.gtypes";

# title = 'trio_asb_intHets_all'
# # title = 'trio_asb_intHets_common'
# na12878.file <- "na12878.asb.combined.all.intHets.sorted.mergedstrict.bed";
# # na12878.file <- "na12878.asb.combined.common.intHets.sorted.mergedstrict.bed";
# na12891.file <- "na12891.asb.combined.all.intHets.sorted.mergedstrict.bed";
# # na12891.file <- "na12891.asb.combined.common.intHets.sorted.mergedstrict.bed";
# na12892.file <- "na12892.asb.combined.all.intHets.sorted.mergedstrict.bed";

# title = 'trio_asb_counts_all'
# title = 'trio_asb_counts_common'
# # na12878.file <- "na12878.asb.combined.all.counts.sorted.mergedStrict.bed";
# na12878.file <- "na12878.asb.combined.common.counts.sorted.mergedStrict.bed";
# # na12891.file <- "na12891.asb.combined.all.counts.sorted.mergedStrict.bed";
# na12891.file <- "na12891.asb.combined.common.counts.sorted.mergedStrict.bed";
# na12892.file <- "na12892.asb.combined.all.counts.sorted.mergedStrict.bed";

# title = 'trio_ase_counts'
# na12878.file <- "na12878.ase.combined.PolyA.counts.sorted.mergedStrict.bed";
# na12891.file <- "na12891.ase.combined.PolyA.counts.sorted.mergedStrict.bed";
# na12892.file <- "na12892.ase.combined.PolyA.counts.sorted.mergedStrict.bed";

title = 'trio_ase_intHets'
na12878.file <- "na12878.ase.combined.PolyA.intHets.sorted.mergedstrict.bed";
na12891.file <- "na12891.ase.combined.PolyA.intHets.sorted.mergedstrict.bed";
na12892.file <- "na12892.ase.combined.PolyA.intHets.sorted.mergedstrict.bed";


## hets snp.calls chr-pos chr pos refallele mat pat child phase
na12878.hets <- read.table(na12878.file, header=F, sep = "\t", stringsAsFactors = F)
na12891.hets <- read.table(na12891.file, header=F, sep = "\t", stringsAsFactors = F)
na12892.hets <- read.table(na12892.file, header=F, sep = "\t", stringsAsFactors = F)
pos     = list(paste(na12891.hets[,1],na12891.hets[,2],na12891.hets[,3],sep="-"),
               paste(na12892.hets[,1],na12892.hets[,2],na12892.hets[,3],sep="-"),
               paste(na12878.hets[,1],na12878.hets[,2],na12878.hets[,3],sep="-"))


## venn by position
library(VennDiagram)
pdf(paste(title,".pdf"))
posVenn = venn.diagram(pos, NULL, fill=c("green", "blue", "red"),
                        cex = 2, cat.fontface=4, main.fontfamily ="Helvetica",
                        alpha = c(0.5,0.5,0.5),
                        category.names=c("NA12891/P", "NA12892/M","NA12878/D"))
grid.draw(posVenn)
dev.off()

## venn by position
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("limma")
# biocLite("statmod")
# library(limma);
# len = length(pos);
# 
# # union list
# all.uniq = unique(unlist(pos));
# # count for each unique label which NA12878 NA12891 NA12892 has
# tmp = matrix(0, ncol=len, nrow=length(all.uniq));
# colnames(tmp) = c("NA12878","NA12891","NA12892");
# rownames(tmp) = all.uniq;
# for(ith in 1:len){
#   tmp[pos[[ith]],ith] <- 1;
# }
# # sum up
# count.sum = apply(tmp,1,sum)
# # count.sum
# # 3 (4-1) way venn
# # if(len<4){
# #   windows();
# #   vennDiagram(vennCounts(tmp), main=title,cex=1);
# # }
# x11(); vennDiagram(vennCounts(tmp), main=title,cex=1);
# 
# # print files
# f1 = c(0);
# ctr= 0;
# for(jth in 1:length(all.uniq)){
#   if(tmp[jth,"NA12878"]==1 && tmp[jth,"NA12891"]==1 && tmp[jth,"NA12891"]==0)
#   {
#     f1[ctr,1] = c(f1[ctr,1],rownames(tmp)[jth]);
#     ctr = ctr + 1;
#   }
# } 