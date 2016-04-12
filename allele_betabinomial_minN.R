setwd('C:/Users/Jieming/Documents/thesis/mark_work/allele_specificity/datasets_accN/datasets_pooled')
# dbinom(5,5,0.5) *2
# dbinom(4,5,0.5) *2
# dbinom(3,5,0.5) *2
# dbinom(2,5,0.5) *2
# dbinom(1,5,0.5) *2
# dbinom(0,5,0.5) *2
########### betabinomial
library(VGAM)
## note that the counts_min5 and intHets rowcount include the header since they are .txt
# filename = "zstatistics_chipseq_pooled.b.sse.p.txt"
filename = "zstatistics_rnaseq_382samples_pooled.b.sse.p.txt"
data = read.table(filename,header=T, stringsAsFactors=F)

## parameters
maxN = 300
a = 0.5
b = data$b.grad
#b = 0.046875 ## solo
p = data$p.choice.bb ## for each dataset, an FDR of 5% (e.g) correspond to this p value
#p = 0.004666015625 ## solo

## for each b do a betabinomial distribution, for maxN
betabin <- function(x)
{
  ## 2 tail
  j = dbetabinom(0,seq(0,maxN),0.5,x)*2
  j[j>1] = 1 
  
#   print (j) ##debug
  return(j)
}

out.b = lapply(as.matrix(b),betabin)

findN <- function(x,y)
{
  k = which.min(x[x>y])
  return(k)
}

out.p = mapply(findN, out.b, p)

## comment these out for solos
data$minN = out.p
write.table(data,file=paste(filename,".minN.new",sep=''), sep="\t",
           row.names=FALSE,quote=FALSE)