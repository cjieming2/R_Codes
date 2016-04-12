setwd("C:/Users/Jieming/Documents/thesis/mark_work/allele_specificity/alleleseq_in_R")
library(VGAM)
### data
filename1 = "counts.txt"
data1 = read.table(filename1, header=T, stringsAsFactors=F)
filename2 = "betabinomial/b_chosen.grad.txt"
data2 = read.table(filename2, header=T, stringsAsFactors=F)

## parameters
p=0.5 ## binomial null p
FDR.thresh = 0.05
iter = 10

## finding the second highest value for the SNP
cACGT = data.frame(cA=data1$cA,cC=data1$cC,cG=data1$cG,cT=data1$cT)
lower = apply(cACGT,1,function(x) sort(x, partial=3)[3])

## find total num = highest + second highest
higher = apply(cACGT,1,max)
total = higher+lower


## empirical tests
b = data2$b.choice
p.bin = apply(data.frame(2 * mapply(pbinom,lower,total,p)),1,function(x) min(x,1))
p.betabin = apply(data.frame(2 * mapply(pbetabinom,lower,total,p,b)),1,function(x) min(x,1))
data1$p.betabin = p.betabin

## simulations
p.thresh = data.frame(seq(0,1,by=0.001))
cutoff <- function(x,y) sum(y<=x)
fp.sim.bin = matrix(0,nrow(p.thresh),1)
fp.sim.betabin = matrix(0,nrow(p.thresh),1)

#p.sim.bin = matrix(0,length(p.bin),iter)
for (j in 1:iter)
{
  choose.sim = apply(data.frame(total),1,function(x) sum(sample(c(0,1),x,replace=TRUE)))
  p.sim.bin = apply(data.frame(2 * mapply(pbinom,choose.sim,total,p)),1,function(x) min(x,1))
  p.sim.betabin = apply(data.frame(2 * mapply(pbetabinom,choose.sim,total,p,b)),1,function(x) min(x,1))
  fp.sim.bin = fp.sim.bin + apply(p.thresh,1,cutoff,y=p.sim.bin)
  fp.sim.betabin = fp.sim.betabin + apply(p.thresh,1,cutoff,y=p.sim.betabin)
}



## FDR.txt
tp.bin = apply(p.thresh,1,cutoff,y=p.bin)+1
tp.betabin = apply(p.thresh,1,cutoff,y=p.betabin)+1
fp.sim.bin = fp.sim.bin / iter
fp.sim.betabin = fp.sim.betabin / iter 
fdr.bin = fp.sim.bin / tp.bin
fdr.betabin = fp.sim.betabin / tp.betabin
p.choice.bin = max(p.thresh[,1][fdr.bin<=FDR.thresh])
p.choice.betabin = max(p.thresh[,1][fdr.betabin<=FDR.thresh])

FDR.txt = data.frame(cbind(p.thresh,tp.bin,fp.sim.bin,fdr.bin,
                           tp.betabin,fp.sim.betabin,fdr.betabin))
colnames(FDR.txt) <- c("pval","P.bin","FP.bin","FDR.bin",
                       "P.betabin","FP.betabin","FDR.betabin")
FDR.txt[is.na(FDR.txt)] <- 0
FDR.txt[FDR.txt == "Inf"] <- 0 

## take in counts.txt and filter
interestingHets.betabinom = data1[data1$p.betabin<=p.choice.betabin,]
write.table(interestingHets.betabinom,file="interestingHets.betabinom.txt", sep="\t",
            row.names=FALSE,quote=FALSE)
write.table(FDR.txt,file="FDR.betabinomial.txt",sep="\t",row.names=FALSE,quote=FALSE)
write(rbind(paste("p.choice.bin=",p.choice.bin),paste("p.choice.betabin=",p.choice.betabin)),
      file="FDR.betabinomial.txt",append=TRUE)
