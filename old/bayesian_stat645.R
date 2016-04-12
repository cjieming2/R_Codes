setwd("C:/Users/JM/yale/courses/courses spring 2011/stat 645 stats mtds in genetics and bioinfo/assignments")

########################################################################################
# this function calculates the probability of the data, P(y)
loglikelihood<-function(data)
{
  counts = table(data)
  a = array(1/length(counts), dim = length(counts))
  
  # log gamma hence hence use sum, + and -. sum over all different "amino acids" present
  dprob = sum(lgamma(a + counts) - lgamma(a)) + lgamma(1) - lgamma(1+sum(counts))
  
  return (dprob)
}

########################################################################################
# calculate log likelihood of I
calcloglikelihoodI<-function(untreated,treated,combined)
{
  indicator = matrix(0,ncol(combined),4)
  loglikelihoodI = 0
  for (i in 1:ncol(combined))
  {
     logl_unt = loglikelihood(untreated[,i])
     logl_t = loglikelihood(treated[,i]) 
     logl_combined = loglikelihood(combined[,i])
     
     indicator[i,1] = logl_combined      #H1
     indicator[i,2] = logl_unt + logl_t  #H2
     if(indicator[i,1]>indicator[i,2]) indicator[i,3]=1 else indicator[i,3]=2
     
     indicator[i,4] = indicator[i,indicator[i,3]]
  }
  return(indicator)
}

########################################################################################
########################################################################################
########################################################################################
# read in files
untreated <- read.table("Untreated_4146.txt",colClasses="character")
treated   <- read.table("IDV_949.txt",colClasses="character")
combined <- rbind(treated,untreated)
labels <- rbind(c(rep("untreated",nrow(untreated)), rep("treated",nrow(treated))))

## (a)
# calculate how many different amino acids at a certain position
desiredPos = 82

# log likelihood
# H1 - from same pmf
# H2 - from different pmf
logP_dataUntreated = loglikelihood(untreated[,desiredPos])
logP_datatreated = loglikelihood(treated[,desiredPos])
logP_combined = loglikelihood(combined[,desiredPos])
summm = logP_dataUntreated + logP_datatreated
print (paste("H1_loglikelihood:",logP_combined))
print (paste("H2_loglikelihood:",summm))

## (b)
# calculate posterior of indicator variable 
# calculate first likelihood of indicators; H1: combined, H2: separate sum
# prior P(H1) = P(H2) = 0.5, P(I) = 0.5*0.5= 0.25
priorI = 0.25
indicator = calcloglikelihoodI(untreated,treated,combined)
loglikelihoodI = sum(indicator[,4])

# likelihoodI = exp(loglikelihoodI)
# posterior = priorI * likelihoodI

# metropolis-hasting
niter = 10000
mhresults <- matrix(0,nrow=niter,ncol=ncol(combined))

# first randomly assign I for each position
mhresults[1,] <- sample(1:2,99,replace=TRUE)
rsampling <- mhresults[1,]