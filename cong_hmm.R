#set up the model
transition = matrix(c(0.99, 0.01, 0.2, 0.8), ncol = 2, nrow = 2, byrow = T)

modelSize = array(c(2, 6))
rownames(modelSize) = c("states", "symbols")
initial = c(0.5, 0.5)
emission = matrix(c(rep(1/6, 6), rep(1/10, 5), 1/2), ncol = 6, nrow = 2, byrow = T)
diceModel = list(modelSize = modelSize, initial = initial, transition = transition, emission = emission)

##############################################################
# some functions that are used throughout the problem

# generate samples(both hidden states and observations)
generateSample = function(model, nStep)
{
	samples = matrix(ncol = 2, nrow = nStep)
	colnames(samples) = c("hiddenStates", "observations")
	samples[1, 1] = sample(2, 1, prob = model$initial)
	for(i in 2:nStep)
	{
		samples[i, 1] = sample(2, 1, prob = model$transition[samples[i-1, 1],])
	}	
	for(i in 1:nStep)
	{
		samples[i, 2] = sample(6, 1, prob = model$emission[samples[i, 1],])
	}
	samples
}


#forward summation at each step
forwardSum = function(model, obs)
{
	summations = matrix(ncol = model$modelSize['states'], nrow = length(obs))

	summations[1, ] = model$emission[, obs[1]]*model$initial

	for(i in 2:dim(summations)[1])
	{
		summations[i,] = model$emission[,obs[i]]*(summations[i-1,]%*%(model$transition))
	}
	summations
}

#calculate the likelihood of all observations using forward summation
forwardSumTotal = function(summations)
{
	sum(summations[dim(summations)[1],])	
}

#backward summation at each step
backwardSum = function(model, obs)
{
	summations = matrix(ncol = model$modelSize['states'], nrow = length(obs))

	summations[dim(summations)[1], ] = model$emission[, obs[dim(summations)[1]]]

	for(i in (dim(summations)[1]-1):1)
	{
		summations[i,] = model$emission[,obs[i]]*(summations[i+1,]%*%t(model$transition))
	}
	summations
}

#calculate the likelihood of all observations using backward summation
backwardSumTotal = function(summations, model)
{
	sum(summations[1,]*model$initial)	
}


#calculate the marginal distribution of hidden states at each step using the results from forward summation and backward summation
getHidStMarginDistribution = function(fSummations, bSummations, model)
{
	sumTotal = forwardSumTotal(fSummations)
	newBSummations = bSummations
	newBSummations[dim(bSummations)[1],] = 1
	for(i in 1:(dim(bSummations)[1]-1))
	{
		newBSummations[i,] = bSummations[i+1,]%*%t(model$transition)
			
	}
	marginDistribution = fSummations*newBSummations/sumTotal
		
}

#use viterbi algorithm to find the most probably path
viterbi = function(model, obs)
{
	maxMatrix = matrix(nrow = length(obs)-1, ncol = model$modelSize['states'])
	traceMatrix = matrix(nrow = length(obs)-1, ncol = model$modelSize['states'])
	maxHidPath = array(dim=length(obs))
	
	probs = log10(model$emission[,obs[1]]) + log10(model$initial) + log10(model$transition)
	maxMatrix[1,] = apply(probs, MARGIN = 2, max)
	traceMatrix[1,] = apply(probs, MARGIN = 2, which.max)
	
	for(i in 2:(length(obs)-1))
	{
		probs = log10(model$emission[,obs[i]]) + maxMatrix[i-1,] + log10(model$transition)
		maxMatrix[i,] = apply(probs, MARGIN = 2, max)
		traceMatrix[i,] = apply(probs, MARGIN = 2, which.max)
	}
	probs = log10(model$emission[,obs[i+1]]) + maxMatrix[i,]
	maxHidPath[length(obs)] = which.max(probs)
	for(i in (length(obs)-1):1)
	{
		maxHidPath[i] = traceMatrix[i,maxHidPath[i+1]]
				
	}
	
	maxHidPath
}

#calculate the MLE of model parameters given observations and marginal hidden states distributions
mle = function(modelSize, initial, obs, currentPar)
{
	fSum = forwardSum(currentPar, obs)
	bSum = backwardSum(currentPar, obs)
	marginDistribution = getHidStMarginDistribution(fSum, bSum, currentPar)

	transition = matrix(0, nrow = modelSize['states'], ncol = modelSize['states'])
	emission = matrix(0, nrow = modelSize['states'], ncol = modelSize['symbols'])
	
	
	# the expected number of transitions
	transition[1, 1] = 	sum(marginDistribution[1:(length(obs)-1),1]* marginDistribution[2:length(obs),1])
	transition[1, 2] = 	sum(marginDistribution[1:(length(obs)-1),1]* marginDistribution[2:length(obs),2])
	transition[2, 1] = 	sum(marginDistribution[1:(length(obs)-1),2]* marginDistribution[2:length(obs),1])
	transition[2, 2] = 	sum(marginDistribution[1:(length(obs)-1),2]* marginDistribution[2:length(obs),2])

	# the MLE of transitions probabilities
	transition = transition + 1e-50 # add a very small number to avoid dividing by 0
	transition = transition/apply(transition, MARGIN=1, sum)
	
	# the expected number of emissions
	emission[1,] = tapply(marginDistribution[,1], INDEX = factor(obs, levels = c(1:modelSize['symbols'])), sum)
	emission[2,] = tapply(marginDistribution[,2], INDEX = factor(obs, levels = c(1:modelSize['symbols'])), sum)

	# the MLE of emission probabilities
	emission[is.na(emission)] = 0
	emission = emission + 1e-50
	emission = emission/apply(emission, MARGIN=1, sum)
	
	model = list(modelSize = modelSize, initial = initial, emission = emission, transition = transition)
	print(model)
	
	return(model)	
}

#use EM algorithm to estimate model parameters(emission and transition probabilities)
em = function(modelSize, obs)
{
	transition = matrix(c(0.5, 0.5, 0.5, 0.5), ncol = 2, nrow = 2, byrow = T)
	emission = matrix(c(3/8, rep(1/8, 5), rep(1/8, 5), 3/8), ncol = 6, nrow = 2, byrow = T)

	initial = c(0.5, 0.5)
	
	currentPar = list(modelSize = modelSize, transition = transition, emission = emission, initial = initial)
	nextPar = mle(modelSize, initial, obs, currentPar)
		
	nIter = 0
	while(!cmpPar(currentPar, nextPar, 1e-6))
	{
		print(currentPar)
		currentPar = nextPar
		nextPar = mle(modelSize, initial, obs, currentPar)
		nIter = nIter + 1
		if(nIter >= 10000)
		{
			print("Failed to converge.\n")
			break	
		}
	}
	
	model = mle(modelSize, initial, obs, nextPar)
}
##########################################################


#compare two different sets of parameters, return TRUE if the maximum difference is smaller than the given precision
cmpPar = function(par1, par2, prec)
{
	maxDif = max(abs(par1$emission - par2$emission), abs(par1$transition - par2$transition))
	if(maxDif < prec)
	{
		return (TRUE)	
	} else
	{
		return (FALSE)	
	}
}
##########################################################


# part b)
samples = generateSample(diceModel, 200)

# part c)
fSum = forwardSum(diceModel, samples[,'observations'])
fSumTotal = forwardSumTotal(fSum)

# part d)
bSum = backwardSum(diceModel, samples[,'observations'])
bSumTotal = backwardSumTotal(bSum, diceModel)

# part e)
marginalHidStDistribution = getHidStMarginDistribution(fSum, bSum, diceModel)
marginalMaxHidStPath = apply(marginalHidStDistribution, MARGIN=1, which.max)

# part f)
maxHidPath = viterbi(diceModel, samples[,2])

# part g)
emModel = em(modelSize, samples[,2])

