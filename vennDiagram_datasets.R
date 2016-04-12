setwd("C:/Users/JM/thesis/mark_work/allele_specificity/datasets/venn")

## filenames
title   = "datasets_venn"
geuv.f <- "geuvadis2013.txt";
mont.f <- "montgomery2010.txt";
kaso.f <- "kasowski2013.txt";
kilp.f <- "kilpinen2013.txt";
mcvi.f <- "mcvicker2013.txt";
pick.f <- "pickrell2010.txt";
enco.f <- "encode2012.txt";


## read files
geuv <- read.table(geuv.f, header=F, sep = "\t", stringsAsFactors = F)
mont <- read.table(mont.f, header=F, sep = "\t", stringsAsFactors = F)
kaso <- read.table(kaso.f, header=F, sep = "\t", stringsAsFactors = F)
kilp <- read.table(kilp.f, header=F, sep = "\t", stringsAsFactors = F)
mcvi <- read.table(mcvi.f, header=F, sep = "\t", stringsAsFactors = F)
pick <- read.table(pick.f, header=F, sep = "\t", stringsAsFactors = F)
enco <- read.table(enco.f, header=F, sep = "\t", stringsAsFactors = F)


# pos     = list(geuv,mont,kaso,kilp,mcvi,pick,enco)
pos     = list(geuv,mont,kaso);


## venn by position
library(VennDiagram)
pdf(paste(title,".pdf"))
# posVenn = venn.diagram(pos, NULL, col="transparent",
#                        fill=c("orange", "green", "red", 
#                               "cyan", "violetred", "blue",
#                               "yellow"),
#                         cex = 2, cat.fontface=4, main.fontfamily ="Helvetica",
#                         alpha = 0.5,
#                         category.names=c("geuvadis2013", "montgomery2010","kasowski2013",
#                                          "kilpinen2013","mcvicker2013","pickrell2010","encode2012"))
posVenn = venn.diagram(pos, 
                       NULL, 
                       col="transparent",
                       fill=c("blue", "green", "red"),
                       alpha = c(0.5,0.5,0.5),
                       cex = 2, 
                       cat.fontface=4, 
                       main.fontfamily ="Helvetica",
                       category.names=c("geuvadis2013", "montgomery2010","kasowski2013"))

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