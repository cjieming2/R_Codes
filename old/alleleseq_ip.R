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
step = 0.01
p.thresh = data.frame(seq(0,1,by=step))
cutoff <- function(x,y) sum(y<=x)
fp.sim.bin = matrix(0,nrow(p.thresh),1)
fp.sim.betabin = matrix(0,nrow(p.thresh),1)
p.sim.bin = matrix(0,nrow(data1),iter)
p.sim.betabin = matrix(0,nrow(data1),iter)

#p.sim.bin = matrix(0,length(p.bin),iter)
for (j in 1:iter)
{
  choose.sim = apply(data.frame(total),1,function(x) sum(sample(c(0,1),x,replace=TRUE)))
  p.sim.bin[,j] = apply(data.frame(2 * mapply(pbinom,choose.sim,total,p)),1,function(x) min(x,1))
  p.sim.betabin[,j] = apply(data.frame(2 * mapply(pbetabinom,choose.sim,total,p,b)),1,function(x) min(x,1))
}

fp.sim.bin = apply(p.thresh,1,cutoff,y=p.sim.bin) / iter
fp.sim.betabin = apply(p.thresh,1,cutoff,y=p.sim.betabin) / iter

## FDR.txt
tp.bin = apply(p.thresh,1,cutoff,y=p.bin)+1
tp.betabin = apply(p.thresh,1,cutoff,y=p.betabin)+1
fdr.bin = fp.sim.bin / tp.bin
fdr.betabin = fp.sim.betabin / tp.betabin
p.choice.bin = max(p.thresh[,1][fdr.bin<=FDR.thresh])
p.choice.betabin = max(p.thresh[,1][fdr.betabin<=FDR.thresh])

fdr.choice.bin = max(fdr.bin[fdr.bin<=FDR.thresh])
fdr.choice.betabin = max(fdr.betabin[fdr.betabin<=FDR.thresh])

## bisection method to find p value
bisect <- function(p,p.sim,iter,p.choice,fdr,fdr.threshold,by)
{
  p.fdr.e = matrix(0,100,3)
  e.prev = 10
  flag = 3
  ctr = 1
  p.fdr.e[ctr,1] = p.choice
  p.fdr.e[ctr,2] = fdr
  p.fdr.e[ctr,3] = e.prev
  
  while(flag)
  {
    start = max(0,(p.choice - by/2))
    end = p.choice + by/2
    by = by/4
    range = seq(start,end,by)
    
    for (i in range)
    {
      tp = cutoff(i,p)
      fp = cutoff(i,p.sim) / iter
      fdr.ind = fp/tp
      e.curr = fdr.threshold - fdr.ind
      ctr = ctr + 1
      
      p.fdr.e[ctr,1] = i
      p.fdr.e[ctr,2] = fdr.ind
      p.fdr.e[ctr,3] = e.curr
      e.prev = p.fdr.e[(ctr-1),3]
      p.choice = i
        
      if(e.curr < 0){ break }
      
#       print(paste(start,"|",end,"|",i,"|",ctr,"|",by)) ##debug
#       print(paste("fdr.thresh=",fdr.threshold,"|fdr=",fdr.ind,"|fdrmatrix=",p.fdr.e[(ctr-1),2],
#                   "e.curr=",p.fdr.e[ctr,3],"|e.prev=",p.fdr.e[ctr-1,3])) ##debug
    }
#     print(paste(start,"|",end,"|",i,"|",ctr,"|",by)) ##debug
#     print(paste("fdr.thresh=",fdr.threshold,"|fdr=",p.fdr.e[ctr,2],"|fdrprev=",p.fdr.e[ctr-1,2])) 
#     break##debug
#     print(paste("tp=",tp,"|fp=",fp)) ##debug
#     print(paste("e.curr=",e.curr,"|e.prev=",e.prev)) ##debug
    
    if(signif(p.fdr.e[ctr-1,3],3) == signif(p.fdr.e[ctr,3],3)){ flag = 0 }
  }  
  return(p.fdr.e)
}
p.choice.bin.1 = bisect(p.bin,p.sim.bin,iter,p.choice.bin,fdr.choice.bin,FDR.thresh,step)
p.choice.betabin.1 = bisect(p.betabin,p.sim.betabin,iter,p.choice.betabin,fdr.choice.betabin,FDR.thresh,step)

p.choice.bin.1 = p.choice.bin.1[p.choice.bin.1[,3]>0,]
p.choice.bin.2 = p.choice.bin.1[nrow(p.choice.bin.1),1]
p.choice.betabin.1 = p.choice.betabin.1[p.choice.betabin.1[,3]>0,]
p.choice.betabin.2 = p.choice.betabin.1[nrow(p.choice.betabin.1),1]

## formatting FDR.txt
FDR.txt = data.frame(cbind(p.thresh,tp.bin,fp.sim.bin,fdr.bin,
                           tp.betabin,fp.sim.betabin,fdr.betabin))
colnames(FDR.txt) <- c("pval","P.bin","FP.bin","FDR.bin",
                       "P.betabin","FP.betabin","FDR.betabin")
FDR.txt[is.na(FDR.txt)] <- 0
FDR.txt[FDR.txt == "Inf"] <- 0 


## take in counts.txt and filter
interestingHets.betabinom = data1[data1$p.betabin<=p.choice.betabin,]

## printing files
write.table(interestingHets.betabinom,file="interestingHets.betabinom.txt", sep="\t",
            row.names=FALSE,quote=FALSE)
write.table(FDR.txt,file="FDR.betabinomial.txt",sep="\t",row.names=FALSE,quote=FALSE)
write(rbind(paste("p.choice.bin.old=",p.choice.bin),paste("p.choice.betabin.old=",p.choice.betabin)),
      file="FDR.betabinomial.txt",append=TRUE)
write(rbind(paste("p.choice.bin=",p.choice.bin.2),paste("p.choice.betabin=",p.choice.betabin.2)),
      file="FDR.betabinomial.txt",append=TRUE)