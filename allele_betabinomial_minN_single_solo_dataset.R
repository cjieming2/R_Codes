# dbinom(5,5,0.5) *2
# dbinom(4,5,0.5) *2
# dbinom(3,5,0.5) *2
# dbinom(2,5,0.5) *2
# dbinom(1,5,0.5) *2
# dbinom(0,5,0.5) *2
########### betabinomial
library(VGAM)

## parameters
maxN = 479960
a = 0.5
b = 0.041796875 ## solo
p = 0.00366650390625 ## solo; for each dataset, an FDR of 5% (e.g) correspond to this p value

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

minN = mapply(findN, out.b, p)