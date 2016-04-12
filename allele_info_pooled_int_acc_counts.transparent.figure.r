setwd("C:/Users/Jieming/Documents/thesis/mark_work/allele_specificity/datasets_accN/datasets_pooled/ase")
source("C:/Users/Jieming/Documents/Shared/scripts-R_perl_shell_macros/R codes/jmRlib.R")

## data
# chipseq = read.table('info.pooled.383samples.chipseq.txt',header=T,stringsAsFactors = F)
# # chipseq = chipseq[ order(chipseq$TF),  ]
# chipseq = chipseq[ order(chipseq$population),  ]
rnaseq  = read.table('rnaseq-info-pop-transparent-plot.txt',header=T,stringsAsFactors=F)
rnaseq = rnaseq[ order(rnaseq$population),  ]

#################### chipseq
# color by TF
# tf = as.data.frame(rainbow(length(unique(chipseq$TF))),unique(chipseq$TF))
# cols = data.frame(stringsAsFactors=F)
# for(i in 1:nrow(chipseq))
# {
#   for(j in 1:nrow(tf))
#   {
#     if(chipseq$TF[i] == row.names(tf)[j])
#     {
#       cols[i,1] = row.names(tf)[j]
#       cols[i,2] = tf[j,1]
#     }
#   }
# }

cols = as.character(nrow(chipseq))
for(i in 1:nrow(chipseq))
{
  if(chipseq$population[i] == "GBR"){cols[i] = 'red'} else
    if(chipseq$population[i] == "FIN"){cols[i] = 'purple'} else
      if(chipseq$population[i] == "CEU"){cols[i] = 'blue'} else
        if(chipseq$population[i] == "YRI"){cols[i] = 'green'} else
          if(chipseq$population[i] == "CHB"){cols[i] = 'orange'} else
            if(chipseq$population[i] == "JPT"){cols[i] = 'yellow'} else
              if(chipseq$population[i] == "TSI"){cols[i] = 'grey50'}
}

# acc counts
x11()
ylimit= c(0,max(chipseq$nCounts.accN))
# ylimit=c(0,20000)
par(mar=c(8,7,4,1),cex.axis=1, cex.lab=2, cex.main=2)
# barplot(chipseq$nCounts.accN, names.arg=paste(chipseq$sample,"_",chipseq$TF), col=add.alpha(cols[,2],alpha=0.2), border=add.alpha(cols[,2],alpha=0.5),xlab="individuals grouped by TF",ylab="number of SNVs", 
#         ylim=ylimit, las=2)
barplot(chipseq$nCounts.accN, names.arg=paste(chipseq$sample,"_",chipseq$TF), col=add.alpha(cols,alpha=0.2), border=add.alpha(cols,alpha=0.5),xlab="individuals grouped by TF",ylab="number of SNVs", 
        ylim=ylimit, las=2)
# axis(1, at=1:nrow(cols), labels =cols[,1],las=2)

legend(2,40000,unique(chipseq$population),col=c("blue","orange","yellow","green"),pch=15)

## acc counts limit
# ylimit = c(0,2500)
# x11()
# par(mar=c(5,5,4,1),cex.axis=1, cex.lab=2, cex.main=2)
# barplot(chipseq$nCounts.accN, col=add.alpha(cols, alpha=0.3), xlab="population",ylab="number of SNVs", ylim=ylimit)
# axis(1, at=c(48,92,135,217,259,310,415), labels =c("CEU","CHB","FIN","GBR","JPT","TSI","YRI"),las=2)

## intHets binomial counts
# ylimit=c(0,2000)
# x11()
# par(mar=c(5,5,4,1),cex.axis=1, cex.lab=2, cex.main=2)
# barplot(chipseq$intHets.bin, col=add.alpha(cols, alpha=0.3), xlab="population",ylab="number of SNVs", ylim=ylimit)
# axis(1, at=c(48,92,135,217,259,310,415), labels =c("CEU","CHB","FIN","GBR","JPT","TSI","YRI"),las=2)


## intHets betabinomial counts
# x11()
# par(mar=c(5,5,4,1),cex.axis=1, cex.lab=2, cex.main=2)
barplot(chipseq$intHets.betabin, col=add.alpha(cols,alpha=0) , border="black", xlab="",ylab="", xaxt="n", yaxt="n", add=T, ylim=ylimit)
# barplot(chipseq$intHets.betabin, col=cols, xlab="population",ylab="number of SNVs")
# axis(1, at=c(48,92,135,217,259,310,415), labels =c("CEU","CHB","FIN","GBR","JPT","TSI","YRI"),las=2)


#################### rnaseq
# color by pop
cols = as.character(nrow(rnaseq))
for(i in 1:nrow(rnaseq))
{
  if(rnaseq$population[i] == "GBR"){cols[i] = 'red'} else
    if(rnaseq$population[i] == "FIN"){cols[i] = 'purple'} else
      if(rnaseq$population[i] == "CEU"){cols[i] = 'blue'} else
        if(rnaseq$population[i] == "YRI"){cols[i] = 'green'} else
          if(rnaseq$population[i] == "CHB"){cols[i] = 'orange'} else
            if(rnaseq$population[i] == "JPT"){cols[i] = 'yellow'} else
              if(rnaseq$population[i] == "TSI"){cols[i] = 'grey50'}
}

# acc counts
x11()
ylimit= c(0,max(rnaseq$ncounts.bb.min6))
# ylimit=c(0,20000)
par(mar=c(5,5,4,1),cex.axis=1, cex.lab=2, cex.main=2)
barplot(rnaseq$ncounts.bb.min6, col=add.alpha(cols,alpha=0.2), border=add.alpha(cols,alpha=0.5),xlab="population",ylab="number of SNVs", ylim=ylimit)
axis(1, at=c(48,92,135,217,259,310,415), labels =c("CEU","CHB","FIN","GBR","JPT","TSI","YRI"),las=2)

## acc counts limit
# ylimit = c(0,2500)
# x11()
# par(mar=c(5,5,4,1),cex.axis=1, cex.lab=2, cex.main=2)
# barplot(rnaseq$nCounts.accN, col=add.alpha(cols, alpha=0.3), xlab="population",ylab="number of SNVs", ylim=ylimit)
# axis(1, at=c(48,92,135,217,259,310,415), labels =c("CEU","CHB","FIN","GBR","JPT","TSI","YRI"),las=2)

## intHets binomial counts
# ylimit=c(0,2000)
# x11()
# par(mar=c(5,5,4,1),cex.axis=1, cex.lab=2, cex.main=2)
# barplot(rnaseq$intHets.bin, col=add.alpha(cols, alpha=0.3), xlab="population",ylab="number of SNVs", ylim=ylimit)
# axis(1, at=c(48,92,135,217,259,310,415), labels =c("CEU","CHB","FIN","GBR","JPT","TSI","YRI"),las=2)


## intHets betabinomial counts
# x11()
# par(mar=c(5,5,4,1),cex.axis=1, cex.lab=2, cex.main=2)
barplot(rnaseq$intHets.bb.min6, col=cols, xlab="",ylab="", add=T, ylim=ylimit)
# barplot(rnaseq$intHets.betabin, col=cols, xlab="population",ylab="number of SNVs")
axis(1, at=c(48,92,135,217,259,310,415), labels =c("CEU","CHB","FIN","GBR","JPT","TSI","YRI"),las=2)
