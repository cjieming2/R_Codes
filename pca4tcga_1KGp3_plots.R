## eigenvec and pop info
eigenvec.pop = merge(pops, eigenvec, by.x = "sample.id", by.y = "X.FID")
eigenvec.pop[eigenvec.pop$dataset == "ICGC" | eigenvec.pop$dataset == "TCGA",]$population = "TCGA_ICGC"

eigenvec.pop[eigenvec.pop$dataset == "ICGC" | eigenvec.pop$dataset == "TCGA",]$superpopulation = paste("TCGA_ICGC_", eigenvec.pop[eigenvec.pop$dataset == "ICGC" | eigenvec.pop$dataset == "TCGA",]$superpopulation, sep="")

## colors super populations
nc.supern = ""
nc.supern[1] = "red" ## AFR
nc.supern[2] = "purple" ## AMR
nc.supern[3] = "green" ## EAS
nc.supern[4] = "cyan" ## EUR
nc.supern[5] = "orange" ## SAS
nc.supern[6] = "black" ## TCGA_ICGC_ADM
nc.supern[7] = "pink" ## TCGA_ICGC_AFR
nc.supern[8] = "magenta" ## TCGA_ICGC_AMR
nc.supern[9] = "yellowgreen" ## TCGA_ICGC_ASN
nc.supern[10] = "blue" ## TCGA_ICGC_EUR
nc.supern[11] = "brown" ## TCGA_ICGC_SAN

## plot superpopulation ####
pmain = ggplot(eigenvec.pop, aes(PC1,PC2))
phisto = list(geom_point(aes(color=superpopulation, shape=superpopulation, size=superpopulation, alpha=superpopulation), position=position_jitter()), geom_point(data=subset(eigenvec.pop, superpopulation == "TCGA_ICGC_SAN"), color="brown", shape=16, size=4, alpha=1))

plabels = labs(x=paste("PC1 (", round(eigenval[1,1]/sum(eigenval)*100,2), "%)", sep=""),
               y=paste("PC2 (", round(eigenval[2,1]/sum(eigenval)*100,2), "%)", sep=""))
paxes = theme(axis.title.x = element_text(face = "bold",colour = "black", size = 32),
              axis.title.y = element_text(face = "bold",colour = "black", size = 32),
              axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 30))
ptitle = theme(plot.title = element_text(lineheight = 50, face = "bold", size = 30))

## segregate TCGA_SAN people

pmain + phisto + plabels + paxes + ptitle + 
  scale_color_manual(values = nc.supern) + 
  scale_shape_manual(values=c(rep(3,5),rep(16,6))) +
  scale_size_manual(values=c(rep(3,5),rep(4,6))) + 
  scale_alpha_manual(values=c(rep(1,5),0.5,rep(1,5))) + 
  theme(legend.title = element_text(size=15, face="bold"), legend.text=element_text(size=15)) 
# ggsave("pca_superpop.pdf", device = "pdf")
