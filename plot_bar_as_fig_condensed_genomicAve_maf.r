# setwd("C:/Users/JM/thesis/mark_work/FIG/eQTLsnps")
setwd("C:/Users/JM/thesis/mark_work/FIG/allelicSNPs")

## data input
# filename1 = "fishersresults_enrichment-tfm-sf-yaoorder.txt"
# filename1 = "fishersresults_enrichment-tfp-sf-yaoorder.txt"
filename1 = "proportions_AS_rare_lt005_maf_noCDS_highcov.txt";

data1 = read.table(filename1,sep="\t",header=FALSE,row.names=1);

## for CDS RNAseq only
# x11(20,7) # adjust this for aspect ratio! device
# par(mar=c(5,15,5,5), lty=1) # lty=0 invis border
# colors = c("darkblue","blue","lightblue") # no0bins darkgreen=chipseq, darkblue=rnaseq
# barplot(data2$V2,xlab="Fraction of rare SNPs", names.arg=c('Non-peak','AS','Non-AS'), horiz=TRUE,
#         cex.lab=1.5,cex.axis=1.5,cex.names=1.5, col=colors,xpd=FALSE, xlim=c(0.01,0.032), width=0.5,
#         las=1, # set hori y axis label
#         space=c(0,0,0)) # space between bars
# # abline(v=0.01261034,col='red',lty=2) # pseudogenes
# # abline(v=0.012406049,col='red',lty=2) # pseudogenes no CDS
# # abline(v=0.013386563,col='red',lty=2) # genomic average no CDS pass1kgmask
# abline(v=0.013171341,col='red',lty=2) # genomic ave no CDS
# tiff(filename = "Rplot%03d.tif", width = 480, height = 480,
#      units = "px", pointsize = 12,
#      bg = "white", res = 300,
#      restoreConsole = TRUE)
# # xlim=c(0.01,0.016),


## barplots
x11(20,7) # adjust this for aspect ratio! device
par(mar=c(5,10,5,5), lty=1) # this adjust the margins, btm,left,top,right; lty=0 invis border

## for peaks
colors3 = c("darkblue","darkblue","darkblue","darkgreen","darkgreen","darkgreen") # 0bins darkgreen=chipseq, darkblue=rnaseq
colors2 = c("darkblue","darkblue","darkgreen","darkgreen") # 0bins darkgreen=chipseq, darkblue=rnaseq

ordered3 <- c(data1["chipseq.nonpeak","V4"],data1["chipseq.nonas","V4"],data1["chipseq.as","V4"],
              data1["rnaseq.nonpeak","V4"],data1["rnaseq.nonas","V4"],data1["rnaseq.as","V4"])

ordered2 <- c(data1["chipseq.nonas","V4"],data1["chipseq.as","V4"],
              data1["rnaseq.nonas","V4"],data1["rnaseq.as","V4"])

barsy3=barplot(ordered3,xlab="Fraction of rare SNPs in a Personal Genome", names.arg=c('r.np','r.nas','r.as','c.np','c.nas','c.as'), horiz=TRUE,
        cex.lab=1.5,cex.axis=1.5,cex.names=1.5, col=colors3,xpd=FALSE, xlim=c(0.01,0.09), width=0.5,
        las=1, # set hori y axis label
        space=c(0,0,0,0.2,0,0)) # space between bars

# barsy2=barplot(ordered2,xlab="Fraction of rare SNPs in a Personal Genome", names.arg=c('r.nas','r.as','c.nas','c.as'), horiz=TRUE,
#               cex.lab=1.5,cex.axis=1.5,cex.names=1.5, col=colors2,xpd=FALSE, xlim=c(0.01,0.07), width=0.5,
#               las=1, # set hori y axis label
#               space=c(0,0,0.2,0)) # space between bars

abline(v=data1["background","V4"],col='red',lty=2) # genomic ave no CDS
# abline(v=data1["pseudogenes","V4"],col='pink',lty=2) # pseudogenes no CDS

## 95% CI
# barsx_lower=rbind(0.01509992,0.01291343,0.01290524,0.01220467)
# barsx_upper=rbind(0.01639706,0.01518349,0.01469161,0.01412036)
# arrows(barsx_lower,barsy,barsx_upper,barsy,length=0.1,angle=90,code=3)