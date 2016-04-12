#Get poverty data

poverty=read.csv("http://www.reuningscherer.net/stat230/datasets/poverty.csv")

#transform select variables and reassemble into a new dataset

newpov=cbind(poverty[,1:2],sqrt(poverty[,3]),sqrt(poverty[,4]), log(poverty[,5:7]),
  log(poverty[,8]-min(poverty[,8])*1.1), sqrt(abs(poverty[,9]))*sign(poverty[,9]),
  poverty[,10:14])

names(newpov)=c("NAME","agricult","sqrtfemilit","sqrtmalilit","logimport","logexport","loggni",
            "loginflat","newforest","urban","fertilit","mortalit","popgrow","lifexp")

#Keep only complete cases
newpov=newpov[complete.cases(newpov),]

write.csv(newpov,"C:/Documents and Settings/Jonathan/My Documents/Classes/Datasets/Poverty/povtransR.csv")

##########################################################
##  Perform Factor Analysis using Maximum Likelihood
##  with Varimax Rotation
##########################################################


fact1=factanal(newpov[-1],factors=3,rotation="varimax")
fact1

#get loading plot for first two factors
plot(fact1$loadings, pch=18, col='red')
abline(h=0)
abline(v=0)
text(fact1$loadings, labels=names(newpov),cex=0.5)

#get loading plot for factors 1 and 3
plot(fact1$loadings[,c(1,3)], pch=18, col='red')
abline(h=0)
abline(v=0)
text(fact1$loadings[,c(1,3)], labels=names(newpov),cex=0.5)

#get reproduced correlation matrix
repro1=fact1$loadings%*%t(fact1$loadings)
#residual correlation matrix
resid1=fact1$cor-repro1
round(resid1,2)

#get root-mean squared residuals
len=length(resid1[upper.tri(resid1)])
RMSR1=sqrt(sum(resid1[upper.tri(resid1)]^2)/len)
RMSR1

#get proportion of residuals greater than 0.05 in absolute value
sum(rep(1,len)[abs(resid1[upper.tri(resid1)])>0.05])/len




##########################################################
##  Perform Factor Analysis using PAF
##  with Varimax Rotation
##########################################################

library(psych)

fact2=factor.pa(newpov[-1],nfactors=3,rotate="varimax")
fact2

#get loading plot for first two factors
plot(fact2$loadings, pch=18, col='red')
abline(h=0)
abline(v=0)
text(fact2$loadings, labels=names(newpov),cex=0.5)

#get loading plot for factors 1 and 3
plot(fact2$loadings[,c(1,3)], pch=18, col='red')
abline(h=0)
abline(v=0)
text(fact2$loadings[,c(1,3)], labels=names(newpov),cex=0.5)

#get reproduced correlation matrix
repro2=fact2$loadings%*%t(fact2$loadings)
#residual correlation matrix
resid2=cor(newpov[,-1])-repro2
round(resid2,2)

#get root-mean squared residuals
len=length(resid2[upper.tri(resid2)])
RMSR2=sqrt(sum(resid2[upper.tri(resid2)]^2)/len)
RMSR2

#get proportion of residuals greater than 0.05 in absolute value
sum(rep(1,len)[abs(resid2[upper.tri(resid2)])>0.05])/len





##########################################################
##  Perform Factor Analysis using iterative PCA
##  with Varimax Rotation
##########################################################

library(psych)

fact3=factor.pa(newpov[-1],nfactors=3,rotate="varimax", SMC=FALSE)
fact3

#get loading plot for first two factors
plot(fact3$loadings, pch=18, col='red')
abline(h=0)
abline(v=0)
text(fact3$loadings, labels=names(newpov),cex=0.5)

#get loading plot for factors 1 and 3
plot(fact3$loadings[,c(1,3)], pch=18, col='red')
abline(h=0)
abline(v=0)
text(fact3$loadings[,c(1,3)], labels=names(newpov),cex=0.5)

#get reproduced correlation matrix
repro3=fact3$loadings%*%t(fact3$loadings)
#residual correlation matrix
resid3=cor(newpov[,-1])-repro3
round(resid3,2)

#get root-mean squared residuals
len=length(resid3[upper.tri(resid3)])
RMSR3=sqrt(sum(resid2[upper.tri(resid3)]^2)/len)
RMSR3

#get proportion of residuals greater than 0.05 in absolute value
sum(rep(1,len)[abs(resid3[upper.tri(resid3)])>0.05])/len


##########################################################
##  Get KMO and other measurements (and also PAF again)
##########################################################

library(rela)
fact4=paf(as.matrix(newpov[,-1]))
fact4  #LOTS of output here!






