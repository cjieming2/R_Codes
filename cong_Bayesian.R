# load the data

treated = read.table("IDV_949.txt", stringsAsFactors = F, colClasses = "character")
untreated = read.table("Untreated_4146.txt", stringsAsFactors = F, colClasses = "character")

#prepare the data
label = factor(c(rep('treated', dim(treated)[1]), rep('untreated', dim(untreated)[1])))
dat = rbind(treated, untreated)

##############################################################
### Define some functions that will be used throughout this problem


#calculated probability of data (single column) given data from H1
H1Prob = function(dat, label)
{
	counts = table(dat)
	a = array(1/length(counts), dim = length(counts))
	
	lprob = sum(lgamma(a + counts) - lgamma(a)) + lgamma(1) - lgamma(1+sum(counts))
	return(lprob) 
}

#calculated probability of data (multiple columns) given data from H1
multiH1Prob = function(dat, label)
{
	dat = as.matrix(dat, nrow = length(label))
	if(dim(dat)[2] == 0)
	{
		return(0)	
	}	
	probs = array(dim = dim(dat)[2])
	for(i in 1:length(probs))
	{
		probs[i] = H1Prob(dat[,i], label)
	}
	
	return(sum(probs))
}

#calculated probability of data (single column) given data from H2
H2Prob = function(dat, label)
{
	levs = levels(label)
	types = unique(dat)
	
	counts1 = getCounts(dat[label == levs[1]], types)
	counts2 = getCounts(dat[label == levs[2]], types)	
	a = array(1/length(types), dim = length(types))
	
	lprob1 = sum(lgamma(a + counts1) - lgamma(a)) + lgamma(1) - lgamma(1+sum(counts1))
	lprob2 = sum(lgamma(a + counts2) - lgamma(a)) + lgamma(1) - lgamma(1+sum(counts2))
	
	return(lprob1 + lprob2)
}


#calculated probability of data (multiple columns) given data from H2
multiH2Prob = function(dat, label)
{
	dat = as.matrix(dat, nrow = length(label))
	if(dim(dat)[2] == 0)
	{
		return(0)
	}
	if(is.null(dim(dat)))
	{
		dim(dat) = c(length(dat), 1)	
	}	
	probs = array(dim = dim(dat)[2])
	
	for(i in 1:length(probs))
	{
		probs[i] = H2Prob(dat[,i], label)
	}
	
	return(sum(probs))		
}

#calculated probability of data given data from H3
H3Prob = function(dat, label)
{
	dat = as.matrix(dat, nrow = length(label))
	if(dim(dat)[2] == 0)
	{
		return(0)
	}
	cDat = matrix(nrow = length(label), ncol = 0)

	for(i in 1:dim(dat)[2])
	{
		cDat = paste(cDat, dat[,i], sep = '')
	}
		
	counts = table(cDat)
	types = unique(cDat)
	counts = counts[types]
	
	dat = data.frame(dat)
	uniqueChars = lapply(dat, unique)
	numPossibleStates = prod(sapply(uniqueChars, length))
	
	a = array(1/numPossibleStates, dim = length(types))

	
	lprob = sum(lgamma(a + counts) - lgamma(a)) + lgamma(1) - lgamma(1+sum(counts))
		
	return(lprob) 	
}

#calculated probability of data given data from H4
H4Prob = function(dat, label)
{
	dat = as.matrix(dat, nrow = length(label))
	if(dim(dat)[2] == 0)
	{
		return(0)
	}
	cDat = matrix(nrow = length(label), ncol = 0)
	for(i in 1:dim(dat)[2])
	{
		cDat = paste(cDat, dat[,i], sep = '')
	}
		
	levs = levels(label)
	types = unique(cDat)
	counts = table(cDat)
	counts = counts[types]
	
	counts1 = getCounts(cDat[label == levs[1]], types)
	counts2 = getCounts(cDat[label == levs[2]], types)

	dat = data.frame(dat)
	uniqueChars = lapply(dat, unique)
	numPossibleStates = prod(sapply(uniqueChars, length))
	

	a = array(1/numPossibleStates, dim = length(types))
	
	lprob1 = sum(lgamma(a + counts1) - lgamma(a)) + lgamma(1) - lgamma(1+sum(counts1))
	lprob2 = sum(lgamma(a + counts2) - lgamma(a)) + lgamma(1) - lgamma(1+sum(counts2))
		
	return(lprob1 + lprob2)	
		
}


#a helper function to obtain the count of each character
getCounts = function(s, states)
{
	counts = array(dim = length(states), dimnames = states)
	for(i in 1:length(states))
	{
		counts[i] = sum(s == states[i])
	}
	return(counts)
}

#do Metropolis-Hastings sampling in the case where only H1 and H2 are considered (used in b)
MH12 = function(dat, label, prior, nIter)
{
	states = c(1, 2)
	trace = numeric(nIter)
	numPos = dim(dat)[2]
	
	preComputedProbs = matrix(nrow = 2, ncol = numPos)
	for(i in 1:numPos)
	{
		preComputedProbs[1,i] = H1Prob(dat[,i], label)
		preComputedProbs[2,i] = H2Prob(dat[,i], label)
	}
	
	samples = matrix(nrow = nIter, ncol = numPos)
	initial = sample(x = states, size = numPos, replace = T)
	s = initial
	curProb = sum(preComputedProbs[1, s == 1]) + sum(preComputedProbs[2, s == 2])
	curProb = curProb + sum(getCounts(s, states)*log(prior))
	
	for(i in 1:nIter)
	{
		print(i)
		proposePos = sample(x = numPos, size = 1)
		restStates = states[states != s[proposePos]]
		propose = s
		propose[proposePos] = sample(x = restStates, size = 1)
		proposeProb = sum(preComputedProbs[1, propose == 1]) + sum(preComputedProbs[2, propose == 2])
		proposeProb = proposeProb + sum(getCounts(propose, states)*log(prior))
		
		if(runif(1) < exp(proposeProb - curProb))
		{
			s = propose
			curProb = proposeProb	
		}
		samples[i,] = s	
		trace[i] = curProb
	}
	return(list(samples = samples, initial = initial, trace = trace))
}


#do Metropolis-Hastings sampling in the case where H1, H2, H3 and H4 are considered (used in d and e)
MH1234 = function(dat, label, prior, nIter)
{	
	states = c(1:4)
	trace = numeric(nIter)
	numPos = dim(dat)[2]
	
	preComputedProbs = matrix(nrow = 2, ncol = numPos)
	for(i in 1:numPos)
	{
		preComputedProbs[1,i] = H1Prob(dat[,i], label)
		preComputedProbs[2,i] = H2Prob(dat[,i], label)
	}

	samples = matrix(nrow = nIter, ncol = numPos)
	
	initial = sample(4, size = numPos, replace = T)
	s = initial
	curProb = sum(preComputedProbs[1, s == 1]) + sum(preComputedProbs[2, s == 2]) + H3Prob(dat[, s==3], label) + H4Prob(dat[, s==4], label)
	curProb = curProb + sum(getCounts(s, states)*log(prior))
	
	for(i in 1:nIter)
	{
		print(i)
		proposePos = sample(x = numPos, size = 1)
		restStates = states[states != s[proposePos]]
		propose = s
		propose[proposePos] = sample(x = restStates, size = 1)
		
		proposeProb = sum(preComputedProbs[1, propose==1]) + sum(preComputedProbs[2, propose==2]) + H3Prob(dat[, propose==3], label) + H4Prob(dat[, propose==4], label)
		proposeProb = proposeProb + sum(getCounts(propose, states)*log(prior))
		
		if(runif(1) < exp(proposeProb - curProb))
		{
			s = propose
			curProb = proposeProb	
		}
		
		samples[i,] = s	
		
		if(sum(samples[i,] == 3) == 1)
		{
			propose[samples[i,] == 3] = 1	
		}
		if(sum(samples[i,] == 4) == 1)
		{
			propose[samples[i,] == 4] = 2	
		}
		
		trace[i] = curProb

	}
	return(list(samples = samples, initial = initial, trace = trace))
		
}
##############################################################


# part a)
aa = tapply(dat[,82], label, table)
sapply(aa, length)

log.p.82.H1 = H1Prob(dat[,82], label)
log.p.82.H2 = H2Prob(dat[,82], label)

#P(H1|data at 82)
p.H1.82 = exp(log.p.82.H1 - log.p.82.H2)/(exp(log.p.82.H1 - log.p.82.H2) + 1)
#P(H2|data at 82)
p.H2.82 = 1/(exp(log.p.82.H1 - log.p.82.H2) + 1)

# part b)
nIter = 20000
samples = MH12(dat = dat, label = label, nIter = nIter, prior = c(0.5, 0.5))

marginal.p = matrix(nrow = 99, ncol = 2)
for(i in 1:2)
{
	marginal.p[,i] = apply(samples[[1]][seq(10001, 20000, 10),] == i, MARGIN = 2, mean)
}

# part c)


log.p.82.54.H1 = multiH1Prob(dat[,c(82,54)], label)
log.p.82.54.H2 = multiH2Prob(dat[,c(82,54)], label)
log.p.82.54.H3 = H3Prob(dat[,c(82,54)], label)
log.p.82.54.H4 = H4Prob(dat[,c(82,54)], label)

#P(H1|data at 82 and 54)
p.H1.82.54 = exp(log.p.82.54.H1 - log.p.82.54.H4)/(exp(log.p.82.54.H1 - log.p.82.54.H4) + exp(log.p.82.54.H2 - log.p.82.54.H4) + exp(log.p.82.54.H3 - log.p.82.54.H4) + 1)

#P(H2|data at 82 and 54)
p.H2.82.54 = exp(log.p.82.54.H2 - log.p.82.54.H4)/(exp(log.p.82.54.H1 - log.p.82.54.H4) + exp(log.p.82.54.H2 - log.p.82.54.H4) + exp(log.p.82.54.H3 - log.p.82.54.H4) + 1)

#P(H3|data at 82 and 54)
p.H3.82.54 = exp(log.p.82.54.H3 - log.p.82.54.H4)/(exp(log.p.82.54.H1 - log.p.82.54.H4) + exp(log.p.82.54.H2 - log.p.82.54.H4) + exp(log.p.82.54.H3 - log.p.82.54.H4) + 1)

#P(H4|data at 82 and 54)
p.H4.82.54 = exp(log.p.82.54.H4 - log.p.82.54.H4)/(exp(log.p.82.54.H1 - log.p.82.54.H4) + exp(log.p.82.54.H2 - log.p.82.54.H4) + exp(log.p.82.54.H3 - log.p.82.54.H4) + 1)


# part d)

nIter = 20000
samples = MH1234(dat = dat, label = label, nIter = nIter, prior = c(0.5, 0.5, 0.5, 0.5))

marginal.p = matrix(nrow = 99, ncol = 4)
for(i in 1:4)
{
	marginal.p[,i] = apply(samples[[1]][seq(1001, 5000, 10),] == i, MARGIN = 2, mean)
}

# part e)

nIter = 20000
samples = MH1234(dat = dat, label = label, nIter = nIter, prior = c(0.48, 0.02, 0.48, 0.02))

marginal.p = matrix(99, 4)
for(i in 1:99)
{
	for(j in 1:4)
	{
		marginal.p[i,j] = sum(samples[[1]][-c(1:10000),i] == k)/(nrow(samples[[1]]) - 10000)	
	}	
}
