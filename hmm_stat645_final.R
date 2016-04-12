setwd("C:/Users/JM/yale/courses/courses spring 2011/stat 645 stats mtds in genetics and bioinfo/assignments")

###############################################################################
# observation generator
# input matrix should be a vector not a matrix
generate<-function(emis_mat)
{ 
  emis_mat <- cumsum(emis_mat)
  obsT = 0
  probT <- runif(1)
  for (k in 1:length(emis_mat))
  {
    if(probT>emis_mat[k])
    {
      next
    }
    obsT = k
    break
  }
  
  # store the probability for debugging and the status (state or obs)
  combined <- cbind(probT,obsT)
  return(combined)
}

###############################################################################3
# function to generate a HMM state (and pipe out the observation)
# returns a list with 2 cols: state and observation
hmmgenerate<-function(nstates,initial_prob,trans_mat,emis_mat)
{
  # first state and observation; 1=fair, 2=loaded with supposed initial prob  
  curr_state <- generate(initial_prob)
  curr_obs <- generate(emis_mat[,curr_state])
  
  # matrices that store all the states and obs
  states <- rbind(curr_state)
  obs <- rbind(curr_obs)
  
  # second state and on
  for (i in 2:nstates)
  {
    # generate the status
    # for state, use row of transition matrix
    # for obs, use column of emission matrix
    curr_state <- generate(trans_mat[curr_state,])
    curr_obs <- generate(emis_mat[,curr_state])
    
    states <- rbind(states,curr_state)
    obs <- rbind(obs,curr_obs)
  }
  
  states[,2][states[,2]==1] <- "fair"
  states[,2][states[,2]==2] <- "loaded"
  finalmat <- data.frame(stateProb=states[,1],states=states[,2],obsProb=obs[,1],observations=obs[,2])
  
  return (finalmat)
}

###################################################################################
# calculate total probability of O, P(O|lambda) using the forward algorithm
# outputs a matrix of fair and loaded probabilities and the sum at each time point
# output gets filled in ascending order (normal)
# obs input is a vector not a matrix
fwd<-function(obs,initial_prob,trans_mat,emis_mat)
{
  # initialize the first observation
  P_ot <- matrix(0,length(obs),(length(initial_prob)+1))
  P_sum = 0
  
  # 1=fair, 2=loaded
  # for the first
  for (i in 1:length(initial_prob))
  {
    P_ot[1,i] = (emis_mat[obs[1],i] * initial_prob[i])
    P_sum = P_sum + P_ot[1,i]
  }

  P_ot[1,(length(initial_prob)+1)] = P_sum
  
  # loop 2 onwards
  for (t in 2:length(obs))
  {
    # initialize sum on each round
    P_sum = 0
    
    # 1=fair 2=loaded
    for (currstate in 1:length(initial_prob))
    {
      P_transum = 0;
      for (prevstate in 1:length(initial_prob))
      {
        P_ot1qt1qt = P_ot[(t-1),prevstate] * trans_mat[prevstate,currstate]
        P_transum = P_transum + P_ot1qt1qt
      }
      P_ot[t,currstate] = P_transum * emis_mat[obs[t],currstate]
      P_sum = P_sum + P_ot[t,currstate]
    }
    
    P_ot[t,(length(initial_prob)+1)] = P_sum
  }
  
  return (P_ot)
}

####################################################################################
# calculate total probability of O P(O|lambda) using the backward algorithm
# outputs a matrix of fair and loaded probabilities and the sum at each time point
# note that since this is backward, the matrix gets filled in descending order
# obs input is a vector not a matrix
bwd<-function(obs,initial_prob,trans_mat,emis_mat)
{    
  # initialize; 1=fair, 2=loaded
  P_ot <- matrix(0,length(obs),(length(initial_prob)+1))
  tempsum = 0
  P_sum = 0
  
  for (i in 1:length(initial_prob))
  {
    for (j in 1:length(initial_prob))
    {
      tempsum = tempsum + emis_mat[obs[length(obs)],i] * trans_mat[i,j]
    }
    P_ot[length(obs),i] = tempsum
    P_sum = P_sum + P_ot[length(obs),i]
  }
  P_ot[length(obs),(length(initial_prob)+1)] = P_sum

  
  # loop
  for (b in (length(obs)-1):2)
  {
    P_sum = 0
    for (cstate in 1:length(initial_prob))
    {
      P_otot1qt2 = 0
      for (estate in 1:length(initial_prob))
      {
        P_otot1qt2 = P_otot1qt2 + (P_ot[(b+1),estate] * emis_mat[obs[b],estate] * trans_mat[cstate,estate])
      }
      P_ot[b,cstate] = P_otot1qt2
      P_sum = P_sum + P_ot[b,cstate]
    }
    P_ot[b,(length(initial_prob)+1)] = P_sum
  }

  P1F = P_ot[2,1] * emis_mat[obs[1],1] * initial_prob[1] 
  P1L = P_ot[2,2] * emis_mat[obs[1],2] * initial_prob[2]
  P1sum = P1F + P1L
  P_ot[1,] = c(P1F,P1L,P1sum)
  
  return (P_ot)
}
  
###############################################################################
# calculate forward-backward algo
fb<-function(obs,initial_prob,trans_mat,emis_mat)
{
  fwdResults<-fwd(obs,initial_prob,trans_mat,emis_mat)
  colnames(fwdResults) <- c("fwdfair","fwdloaded","fwdsumProb")
  bwdResults<-bwd(obs,initial_prob,trans_mat,emis_mat)
  colnames(bwdResults) <- c("bwdfair","bwdloaded","bwdsum")
  
  fwdResults = rbind(fwdResults,c(rep(0,ncol(fwdResults))))
  bwdResults = rbind(c(rep(0,ncol(bwdResults))),bwdResults)
  fbstates <- c()
  fairfb <- (fwdResults[,"fwdfair"] + bwdResults[,"bwdfair"])/(fwdResults[,"fwdsumProb"])
  loadfb <- (fwdResults[,"fwdloaded"] + bwdResults[,"bwdloaded"])/(fwdResults[,"fwdsumProb"])
  
  fbstates <- fairfb>loadfb
  fbstates[fbstates=="TRUE"] <- "fair"
  fbstates[fbstates=="FALSE"] <- "loaded"
  fbstates <- cbind(fairfb,loadfb,fbstates)
  fbstates <- fbstates[1:(nrow(fbstates)-1),1:ncol(fbstates)]
  return(fbstates)
}

###############################################################################
# viterbi algo on HMM
viterbi<-function(obs,initial_prob,trans_mat,emis_mat)
{
  v <- matrix(0,length(obs),length(initial_prob))
  p <- matrix(0,length(obs),length(initial_prob))
  
  # initialize: the first state
  for (i in 1:length(initial_prob))
  {
    v[1,i] = initial_prob[i] * emis_mat[obs[1],i]
  }
  
  for (t in 2:length(obs))
  {    
    # 1=fair, 2=loaded
    for (currstate in 1:length(initial_prob)) #qt
    {
      temp2 = 0
      maxi = 0
      maxistate = 0
      for(prevstate in 1:length(initial_prob)) #qt-1
      {
        # prev viterbi probability * transition prob
        temp2 = v[(t-1),prevstate] * trans_mat[prevstate,currstate]
        
        if(temp2 > maxi)
        {
          maxi = temp2
          maxistate = prevstate
        }
      }
      v[t,currstate] = maxi * emis_mat[obs[t],maxistate]
      p[(t-1),currstate] = maxistate
    }
  }
  
  lastmaxi = 0
  lastmaxistate = 0
  for(laststate in 1:length(initial_prob))
  {
    if(v[length(obs),laststate] > lastmaxi)
    {
      lastmaxi = v[length(obs),laststate]
      lastmaxistate = laststate
    }
  }
  p[length(obs),lastmaxistate] = lastmaxistate
  finalp = p[,lastmaxistate]
  
  vp <- cbind(v,finalp)
  return (vp)
}
###############################################################################
###############################################################################
###############################################################################
# MAIN FUNCTION
## generate hidden states and observations for a fair and loaded dice
## fair die - 6-sided, loaded die - 1/10 for 5 sides and 1/2 for 6th

## (a)
# define variables
nstates = 200
initial_prob = c(0.5,0.5)

# transition matrix
trans_mat = matrix(c(0.99,0.2,0.01,0.8),ncol=2)
row.names(trans_mat) <- c("fair","loaded")
colnames(trans_mat) <- c("fair","loaded")

# emission matrix
emis_mat = matrix(c(rep(1/6,6),rep(1/10,5),0.5),ncol=2)
row.names(emis_mat) <- c(1:6)
colnames(emis_mat) <- c("fair","loaded")

## (b)
# generate hmm model
hmm <- hmmgenerate(nstates,initial_prob,trans_mat,emis_mat)

## (c)
# use forward and backward algo to calc P(O|lambda)
# use f-b algo to estimate states
fwdResults <- fwd(hmm$observations,initial_prob,trans_mat,emis_mat)
colnames(fwdResults) <- c("fwdfair","fwdloaded","fwdsumProb")

## (d)
bwdResults <- bwd(hmm$observations,initial_prob,trans_mat,emis_mat)
colnames(bwdResults) <- c("bwdfair","bwdloaded","bwdsum")

## (e)
fbResults <- fb(hmm$observation,initial_prob,trans_mat,emis_mat)

# calculate percentage accuracy of state estimation by f-b algo
truepos = length(hmm$states[hmm$states == fbResults[,"fbstates"]])
percentAcc = truepos / nrow(hmm) * 100

## (f)
# viterbi algo
viterbiResults <- viterbi(hmm$observations,initial_prob,trans_mat,emis_mat)
colnames(viterbiResults) <- c("vfair","vloaded","vpath")
viterbiResults[,"vpath"][viterbiResults[,"vpath"]==1] <- "fair"
viterbiResults[,"vpath"][viterbiResults[,"vpath"]==2] <- "loaded"

# calculate percentage accuracy of state estimation by Viterbi
trueposVit = length(hmm$states[hmm$states == viterbiResults[,"vpath"]])
percentAccVit = trueposVit / nrow(hmm) * 100

# compare Viterbi with f-b algo
trueposVit_fb = length(viterbiResults[,"vpath"][fbResults[,"fbstates"] == viterbiResults[,"vpath"]])
percentAccVit_fb = trueposVit_fb / nrow(viterbiResults) * 100

## (g)
# EM algorithm
# start off with the following initial "arbitrary" states and parameters (as per (a)):
# number of observations/states = 200
# initial_prob = c(0.5,0.5)
# trans_mat = matrix(c(0.99,0.2,0.01,0.8),ncol=2)
# emis_mat = matrix(c(rep(1/6,6),rep(1/10,5),0.5),ncol=2)
# only the observations are known and fixed
em_initial_prob = initial_prob
em_trans_mat = trans_mat
new_em_trans_mat = matrix(0,nrow(em_trans_mat),ncol(em_trans_mat))
em_emis_mat = emis_mat
new_em_emis_mat = matrix(0,nrow(em_emis_mat),ncol(em_emis_mat))
ctr = 0
tratio_trans = c() 
tratio_emis = c()

# variables
n_obs = matrix(0,nrow(em_emis_mat),nstates)
n_states = matrix(0,nrow(em_trans_mat),nstates)
while(1)
{  
  # counter
  ctr = ctr + 1
    
  # use forward-backward algo to find path
  em_fbResults <- fb(hmm$observations,em_initial_prob,em_trans_mat,em_emis_mat)
  old_path = em_fbResults[,3]
  
  ## reconstruct new emission matrix
  ## find MLE
  # boolean matrices of each observable data in set of observations
  for (z in 1:nrow(em_emis_mat))
  {
    n_obs[z,] = (hmm$observations==z)
  }
  
  # number of fair and loaded states
  for (y in 1:nrow(em_trans_mat))
  {
    n_states[y,] = (em_fbResults[,3]==rownames(em_trans_mat)[y])
  }
  
  for (a in 1:nrow(em_trans_mat))
  {
    for (b in 1:nrow(em_emis_mat))
    {
      # overlay matrices: number of conditional (obs|state); both need to be TRUE 
      os = n_states[a,] + n_obs[b,]
      
      # find number of both-TRUEs
      # b 1|, 2|, 3|, 4|, 5| ,6|  ----a |fair |loaded
      new_em_emis_mat[b,a] = length(n_states[a,][os==2])/length(n_states[a,][n_states[a,]==1])
    }
  }
  
  row.names(new_em_emis_mat) <- c(1:6)
  colnames(new_em_emis_mat) <- c("fair","loaded")
  
  ## reconstruct new transition matrix
  ff_ll = c(rep(0,(length(old_path)-1)))
  for (st in 1:(length(old_path)-1))
  {
    ff_ll[st] = paste(old_path[st],old_path[(st+1)],sep=",")
  }
  
  totaldenom1 = length(ff_ll[ff_ll=="fair,fair"]) + length(ff_ll[ff_ll=="fair,loaded"])
  totaldenom2 = length(ff_ll[ff_ll=="loaded,loaded"]) + length(ff_ll[ff_ll=="loaded,fair"])
  new_em_trans_mat[1,1] = length(P_ff_ll[P_ff_ll=="fair,fair"])/totaldenom1
  new_em_trans_mat[2,2] = length(P_ff_ll[P_ff_ll=="loaded,loaded"])/totaldenom2
  new_em_trans_mat[1,2] = 1-new_em_trans_mat[1,1]
  new_em_trans_mat[2,1] = 1-new_em_trans_mat[2,2]
  row.names(new_em_trans_mat) <- c("fair","loaded")
  colnames(new_em_trans_mat) <- c("fair","loaded")
  
  # threshold
  tratio_trans[ctr] = mean(new_em_trans_mat / em_trans_mat)
  tratio_emis[ctr] = mean(new_em_emis_mat / em_emis_mat)
  
  if((0.99<=tratio_trans[ctr] && tratio_trans[ctr]<=1) && (0<=tratio_emis[ctr] && tratio_emis[ctr]<=1))
  {
    break
  }
  else
  {
    em_emis_mat = new_em_emis_mat
    em_trans_mat = new_em_trans_mat
  }
  
#   break #debug
}