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
calcloglikelihoodI<-function(untreated,treated,combined,states)
{
  indicator = c(rep(0,length(states)))
  loglikelihoodI = 0
  
  for (i in 1:length(states))
  {
    #### H1
    if(states[i]==1)
    {
      logl_combined = loglikelihood(combined[,i])
      indicator[i] = logl_combined      #H1
    }
    else
    {
      #### H2
      if(states[i]==2)
      {
        logl_unt = loglikelihood(untreated[,i])
        logl_t = loglikelihood(treated[,i]) 
        indicator[i] = logl_unt + logl_t  #H2
      }
    }
#     else
#     {
#       #### H3
#       if(states[i]==3)
#       {
#         hyp3.combined = rbind(hyp3.combined,combined[,i])
#       }
#     }
#     else
#     {
#       #### H4
#       if(states[i]==4)
#       {
#         hyp4.untreated = rbind(hyp4.untreated,untreated[,i])
#         hyp4.treated = rbind(hyp4.treated,treated[,i])
#       }
#     }
  }
  
  
  return(indicator)
}

########################################################################################
# gives the state of the indicator H1,H2,H3,H4
generate<-function(x,prob)
{ 
  prob <- cumsum(prob)
  for (k in 1:length(prob))
  {
    if(x>prob[k])
    {
      next
    }
    return (k)
  }
}
########################################################################################
########################################################################################
########################################################################################
# MAIN
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
### metropolis-hasting H1H2
# variables
niter = 10
burnin = 3
mh <- matrix(0,nrow=niter,ncol=ncol(combined))
plotData = array(0,dim=niter)
logprob1 <- array(0,dim=ncol(combined))
logprob2 <- array(0,dim=ncol(combined))

# first randomly assign I for each position
prob = c(0.5,0.5)
indicators = runif(ncol(combined),0,1)
mh[1,] = sapply(indicators,generate,prob)

# run through N steps of markov chain
for (i in 2:niter)
{
  # choose a random i
  mh[i,] = mh[(i-1),]
  
  j = sample(1:ncol(combined),1)
  x = mh[i,] # a vector
  y = mh[i,]
  y[j] = 3-x[j] # complement of state x, ie if x=1,y=2, vice versa
  
  Px = sum(calcloglikelihoodI(untreated,treated,combined,x))
  Py = sum(calcloglikelihoodI(untreated,treated,combined,y))
  Pratio = Py-Px
  
  alpha = min(0,Pratio)
  u = runif(1,0,1)
  if(log(u) <= alpha){x.next = y}
  else {x.next = x}
  
  mh[i,] = x.next
  plotData[i+1] = plotData[i] - Px + sum(calcloglikelihoodI(untreated,treated,combined,x.next))
  
  # if past burnin, update posterior
  if(i>burnin){logprob1 = logprob1 + (mh[i,]==1)}
  if(i>burnin){logprob2 = logprob2 + (mh[i,]==2)}
}

print (paste("P(I=1|Data):"))
print (round((logprob1/(niter-burnin)),2))
print (paste("P(I=2|Data):"))
print (round((logprob2/(niter-burnin)),2))
x11()
plot(plotData)

## (c)
# for positions 82 and 54
desiredPos2 = 54
logP_dataUntreated2 = loglikelihood(untreated[,desiredPos2])
logP_datatreated2 = loglikelihood(treated[,desiredPos2])
logP_combined54 = loglikelihood(combined[,desiredPos2])
summm1 = logP_combined + logP_combined54
summm2 = logP_dataUntreated + logP_datatreated + logP_dataUntreated2 + logP_datatreated2

# log likelihood
# H1 - independent but from same pmf
# H2 - independent but from different pmf
print (paste("H1_loglikelihood82.54:",summm1))
print (paste("H2_loglikelihood82.54:",summm2))

# H3 - dependent but from same pmf
combined82_54 = rbind(combined[,desiredPos],combined[,desiredPos2])
logP_combined82_54 = loglikelihood(combined82_54)

# H4 - dependent ubt from different pmf
dataUntreated82_54 = rbind(untreated[,desiredPos],untreated[,desiredPos2])
datatreated82_54 = rbind(treated[,desiredPos],treated[,desiredPos2])
logP_dataUntreated82_54 = loglikelihood(dataUntreated82_54)
logP_datatreated82_54 = loglikelihood(datatreated82_54)
summm82_54 = logP_dataUntreated82_54 + logP_datatreated82_54

print (paste("H3_loglikelihood82.54:",logP_combined82_54))
print (paste("H4_loglikelihood82.54:",summm82_54))

## (d)
### metropolis-hasting H1H2H3H4
# variables
mhd <- matrix(0,nrow=niter,ncol=ncol(combined))
plotDatad = array(0,dim=niter)
logprob1d <- array(0,dim=ncol(combined))
logprob2d <- array(0,dim=ncol(combined))

# first randomly assign I for each position
probd = c(0.25,0.25,0.25,0.25)
indicatorsd = runif(ncol(combined),0,1)
mhd[1,] = sapply(indicators,generate,probd)

# run through N steps of markov chain
for (i in 2:niter)
{
  # choose a random i
  mh[i,] = mh[(i-1),]
  
  # H1, H2
  j = sample(1:ncol(combined),1)
  x = mh[i,] # a vector
  y = mh[i,]
  y[j] = 3-x[j] # complement of state x, ie if x=1,y=2, vice versa
  
  Px = sum(calcloglikelihoodI(untreated,treated,combined,x))
  Py = sum(calcloglikelihoodI(untreated,treated,combined,y))
  Pratio = Py-Px
  
  alpha = min(0,Pratio)
  u = runif(1,0,1)
  if(log(u) <= alpha){x.next = y}
  else {x.next = x}
  
  mh[i,] = x.next
  plotData[i+1] = plotData[i] - Px + sum(calcloglikelihoodI(untreated,treated,combined,x.next))
  
  # H3, H4
#   hyp3
  
  # if past burnin, update posterior
  if(i>burnin){logprob1 = logprob1 + (mh[i,]==1)}
  if(i>burnin){logprob2 = logprob2 + (mh[i,]==2)}
}

print (paste("P(I=1|Data):"))
print (round((logprob1/(niter-burnin)),2))
print (paste("P(I=2|Data):"))
print (round((logprob2/(niter-burnin)),2))
plot(plotData)