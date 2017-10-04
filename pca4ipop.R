setwd("/Users/jiemingchen/Documents/varimed/ipop/pca")

## my own library
source("/Users/jiemingchen/R_codes/jmRlib.R")
library(ggplot2)
library(RColorBrewer)
library(dplyr)

## input
eigenvec = read.delim("plink2.eigenvec.idcorr", header = T, sep = "\t", stringsAsFactors = FALSE, na.strings = "")

eigenval = read.delim("plink2.eigenval", header = F, sep = "\t", stringsAsFactors = FALSE, na.strings = "")
names(eigenval) = "eigenvalue"

## sample 
## sample-id	population	superpopulation	dataset
## HG00096	GBR	EUR	1KG-p3
pops = read.delim("samples.1KGp3.ipop.list", header=T, sep="\t", stringsAsFactors=FALSE, na.strings = "")

## eigenvec and pop info
eigenvec.pop = merge(pops, eigenvec, by.x = "sample.id", by.y = "X.FID")

## choosing colors ####
colorCount = length(unique(eigenvec.pop$population))

nc="" # AFR - red ; EAS - green ; EUR - blue ; AMR - purple ; SAS - orange
nc[1] = "magenta" #ACB
nc[2] = "NA" #######AFR
nc[3] = "NA" #######AMR
nc[4] = "tomato"  #ASW
nc[5] = "orange"  #BEB
nc[6] = "yellowgreen" #CDX
nc[7] = "cyan"        #CEU
nc[8] = "green"       #CHB
nc[9] = "turquoise"   #CHS
nc[10] = "purple"      #CLM
nc[11] = "NA" #######EAS
nc[12] = "violetred1"  #ESN
nc[13] = "NA" #######EUR
nc[14] = "skyblue1"   #FIN
nc[15] = "slateblue"  #GBR
nc[16] = "tan1"       #GIH
nc[17] = "maroon1"     #GWD
nc[18] = "slategray2"  #IBS
# nc[15] = "black"          #ICGC_TGCA
nc[19] = "lightsalmon" #ITU
nc[20] = "springgreen" #JPT
nc[21] = "olivedrab"   #KHV
nc[22] = "mediumvioletred" #LWK
nc[23] = "sienna"          #MSL
nc[24] = "mediumorchid"    #MXL
nc[25] = "magenta4" #PEL
nc[26] = "orange3" #PJL
nc[27] = "violetred4" #PUR
nc[28] = "NA" ########SAS
nc[29] = "lightsalmon3" #STU
nc[30] = "steelblue1" #TSI
nc[31] = "yellow" #YRI

## plot population ####
x11(type="cairo")
pmain = ggplot(eigenvec.pop, aes(PC1,PC2))
phisto = geom_point(aes(color=population), position=position_jitter(), alpha=0.5, size = 3)
plabels = labs(x=paste("PC1 (", round(eigenval[1,1]/sum(eigenval)*100,2), "%)", sep=""),
               y=paste("PC2 (", round(eigenval[2,1]/sum(eigenval)*100,2), "%)", sep=""))
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))

## based on only 2583 snps
p=pmain + phisto + plabels + paxes + ptitle + scale_color_manual(values = nc) +
  geom_point(data=eigenvec.pop %>% filter(dataset=="IPOP"), aes(x=PC1, y=PC2, shape = population),size = 3, position=position_jitter()) + 
  guides(shape=guide_legend(title="IPOP population"))
print(p)
# ggsave("1kgp3-72ind-ipop-2583snps-full-pop.pdf")

## 2 EUR in Asian clusters, 3 in SAS cluster
eigenvec.pop %>% filter(PC2 < -0.0125, dataset == "IPOP", population == "EUR") 
eigenvec.pop %>% filter(PC2 < (0.0125) & PC2 > 0, dataset == "IPOP", population == "EUR") 


############################################
## superpopulation plot
## colors super populations
nc.super = ""
nc.super[1] = "red" ## AFR
nc.super[2] = "purple" ## AMR
nc.super[3] = "green" ## EAS
nc.super[4] = "cyan" ## EUR
nc.super[5] = "orange" ## SAS

## plot superpopulation ####
x11(type="cairo")
pmain = ggplot(eigenvec.pop, aes(PC1,PC2))
phisto = geom_point(aes(color=superpopulation), position=position_jitter(), alpha=0.5, size=3)
plabels = labs(x=paste("PC1 (", round(eigenval[1,1]/sum(eigenval)*100,2), "%)", sep=""),
               y=paste("PC2 (", round(eigenval[2,1]/sum(eigenval)*100,2), "%)", sep=""))
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 20),
              axis.title.y = element_text(face = "bold",colour = "black", size = 20),
              axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))

p=pmain + phisto + plabels + paxes + ptitle + scale_color_manual(values = nc.super) +
  geom_point(data=eigenvec.pop %>% filter(dataset=="IPOP"), aes(x=PC1, y=PC2, shape = population),size=3, position=position_jitter()) +
  guides(shape=guide_legend(title="IPOP population"))
print(p)
# ggsave("1kgp3-72ind-ipop-2583snps-superpop.pdf")

## 2 EUR in Asian clusters, 3 in SAS cluster
eigenvec.pop %>% filter(PC2 < -0.0125, dataset == "IPOP", population == "EUR") 
eigenvec.pop %>% filter(PC2 < (0.0125) & PC2 > 0, dataset == "IPOP", population == "EUR") 

