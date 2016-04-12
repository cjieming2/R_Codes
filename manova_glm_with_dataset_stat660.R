setwd('C:/Users/JM/yale/courses/courses spring 2011/stats 660 multivariate stats in social and env sciences/our dataset')
# setwd('C:/Users/JM/yale/courses/courses spring 2011/stats 660 multivariate stats in social and env sciences/homework/HW5-manova-glm')

#A function to make chi-square quantile plots 
#to test for multivariate normality of data or residuals
reschisqplot<-function(vars,label){
  #usually, vars is xxx$residuals or data from one group and label is for plot
  x<-cov(scale(vars),use="pairwise.complete.obs")
  squares<-sort(diag(as.matrix(scale(vars))%*%solve(x)%*%as.matrix(t(scale(vars)))))
  quantiles<-quantile(squares)
  hspr<-quantiles[4]-quantiles[2]
  cumprob<-c(1:length(vars[,1]))/length(vars[,1])-1/(2*length(vars[,1]))
  degf<-dim(x)[1]
  quants<-qchisq(cumprob,df=degf)
  gval<-(quants**(-1+degf/2))/(exp(quants/2)*gamma(degf/2)*(sqrt(2)**degf))
  scale<-hspr / (qchisq(.75,degf)-qchisq(.25,degf))
  se<-(scale/gval)*sqrt(cumprob*(1-cumprob)/length(squares))
  lower<-quants-2*se
  upper<-quants+2*se
  
  plot(quants,squares,col='red',pch=19,cex=1.2,xlab="Chi-Square Quantiles",
       ylab=label,main=paste("Chi-Square Quantiles for",label),ylim=range(upper,lower),xlim=range(c(0,quants)))
  lines(c(0,100),c(0,100),col=1)
  lines(quants,upper,col="blue",lty=2,lwd=2)
  lines(quants,lower,col="blue",lty=2,lwd=2)
  legend(0,range(upper,lower)[2]*.9,c("Data","95% Conf Limits"),lty=c(0,2),col=c("red","blue"),lwd=c(2,2),
         pch=c(19,NA))
}

## input data
# mean.rm.na <- function(x) mean(x,na.rm=T)
# file = "ohiocrimehm.txt"
file = "E-GEOD-8536-processed-data-1631482821-t-abridged.txt"
data <- read.table(file,header=T,sep="\t")

# dataname = "Ohio Crime Survey 1995"
dataname = "Yeast Expression Data"

## normal and covariance equal
# x11()
# par(mfrow=c(2,2))
# reschisqplot(data$education,label=paste(dataname,"education"))
# reschisqplot(data$gender,label=paste(dataname,"gender"))

## make interaction plots to see how the factors interact
x11()
par(mfrow=c(2,2))
replicate=data$rep
interaction.plot(data$hr,replicate,data$X10205_at,
                 lwd=3,col=c("red","blue","black"),ylim=c(4,7),
                 xlab="time (hr)",ylab="mean intensity",fixed=TRUE,xtick=TRUE,
                 main=paste("Interaction between time and replicate in",dataname, "at 10205_at"),
                 leg.bty="o")

interaction.plot(data$hr,replicate,data$X10222_at,
                 lwd=3,col=c("red","blue","black"),ylim=c(4,7),
                 xlab="time (hr)",ylab="mean intensity",fixed=TRUE,xtick=TRUE,
                 main=paste("Interaction between time and replicate in",dataname, "at 10222_at"),
                 leg.bty="o")

interaction.plot(data$hr,replicate,data$X10337_at,
                 lwd=3,col=c("red","blue","black"),ylim=c(4,7),
                 xlab="time (hr)",ylab="mean intensity",fixed=TRUE,xtick=TRUE,
                 main=paste("Interaction between time and replicate in",dataname, "at 10337_at"),
                 leg.bty="o")

interaction.plot(data$hr,replicate,data$X5733_at,
                 lwd=3,col=c("red","blue","black"),ylim=c(4,7),
                 xlab="time (hr)",ylab="mean intensity",fixed=TRUE,xtick=TRUE,
                 main=paste("Interaction between time and replicate in",dataname, "at 5733_at"),
                 leg.bty="o")

## MANOVA
# fit linear model
lin.mod = as.matrix(data[,1:4]) ~ data$gender + data$education + 
  data$gender*data$education
mod1 = manova(lin.mod)

# get univariate results
summary.aov(mod1)

# get multivariate results
summary.manova(mod1)
summary.manova(mod1, test="Wilks")

## contrast
require(contrast)
repclass <- factor(data$rep)
hrclass <- factor(data$hr)
ex1.lin.mod <- lm(data[,3] ~ repclass + hrclass + repclass*hrclass, data=data)
contrast(ex1.lin.mod, 
         a=list(repclass="1",hrclass=levels(hrclass)),
         b=list(repclass="2",hrclass=levels(hrclass)))

## make plots to see linaer rships between covariates and responses
x11()
pairs(data[,1:6])

## fit linear model
fit.lin.mod = lm(as.matrix(data[,3:6]) ~ data$rep + data$hr + data$rep*data$hr)
summary.lm(fit.lin.mod)

## plot chisquare plot of residuals
x11()
reschisqplot(fit.lin.mod$residuals,label=paste(dataname,"MANOVA Residuals"))