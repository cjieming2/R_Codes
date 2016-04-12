setwd('C:/Users/JM/yale/courses/courses spring 2011/stats 660 multivariate stats in social and env sciences/our dataset')

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
     ylab=label,main=paste("Chi-Square Quantiles for",label),ylim=range(upper,lower))
    lines(c(0,100),c(0,100),col=1)
    lines(quants,upper,col="blue",lty=2,lwd=2)
    lines(quants,lower,col="blue",lty=2,lwd=2)
    legend(locator(1),c("Data","95% Conf Limits"),lty=c(0,2),col=c("red","blue"),lwd=c(2,2),
      pch=c(19,NA))
}

## input data
file = "E-GEOD-8536-processed-data-1631482821.txt"
data <- read.table(file, header=TRUE, row.names=1, sep="\t")
dataname = "Yeast expression data"

## transform data
data <- sqrt(data)

## check univariate normality (for each column) using boxplot
x11()
boxplot(data)
title(paste(("Boxplot for each col of "),dataname))

## check univariate normality (for each column) using qq plot
x11()
op <- par(mar=rep(2,4)) # remove margins in the figure
plot.new()
subplotparam = ceiling(sqrt(length(data))) # setting the size of subplots
par(mfrow=c(subplotparam-1,subplotparam)) # plot subplots

for (i in 2:length(data))
{
  qqnorm(data[,i],col='red',pch=19,cex=1.2,
         main=paste("QQ Plot ",dataname," col",i), 
         plot.it=1, xlab="theoretical normal quantiles", ylab="data quantiles")
  
  qqline(data[,i],lwd=2)
}

## perform chisq quantile plot to test multivariate normality
x11()
reschisqplot(data,label=dataname)

## compute correlation matrix
corrmatrix <- cor(data,method="spearman")

## perform PCA using correlation matrix
pcaresults <- princomp(data,cor=TRUE,scores=TRUE)

## decide which PC to keep - screeplot
screeplot(pcaresults,type="lines",col="red",lwd=2,pch=19,cex=1.2,
          main=paste("Scree Plot of PCA results of ",dataname))

## decide which PC to keep - eigenvalues (R gives sd, so need to sq)
pcaeigenval <- pcaresults$sdev^2

## looking at the loadings
pcaloadings <- pcaresults$loadings

## plot pairwise PCs
colors <- c("green", "red", "blue","purple","magenta","pink","cyan") ##Set of colours to use for all plots
cats <- matrix(c("12h","rep1","340h","rep2","48h","rep2","60h","rep1","24h","rep1","48h","rep3",
                 "1h","rep2","48h","rep1","24h","rep2","340h","rep3",
                 "120h","rep1","1h","rep3","12h","rep3","24h","rep3","120h","rep3",
                 "1h","rep1","60h","rep3","340h","rep1","12h","rep2","60h","rep2",
                 "120h","rep2"),ncol=2,byrow=TRUE,
               dimnames=list(c("X12.hour.rep.1","X340.hour.rep.2","X48.hour.rep.2",
                               "X60.hour.rep.1","X24.hour.rep.1","X48.hour.rep.3",
                               "X1.hour.rep.2","X48.hour.rep.1","X24.hour.rep.2",
                               "X340.hour.rep.3","X120.hour.rep.1","X1.hour.rep.3",
                               "X12.hour.rep.3","X24.hour.rep.3","X120.hour.rep.3",
                               "X1.hour.rep.1","X60.hour.rep.3","X340.hour.rep.1",
                               "X12.hour.rep.2","X60.hour.rep.2","X120.hour.rep.2"),
                             c("time","rep")))
pcamerged <- merge(cats,pcaloadings[,1:dim(pcaloadings)[2]],by=0)
pcascores <- pcaresults$scores
x11()
pairs(pcascores[,1:5], main="", pch = 20)
x11()
categories <- pcamerged$time
pairs(pcamerged[,4:8],col = colors[unclass(factor(categories))],pch=20,cex=1.5)
x11()
xaxis <- pcamerged$Comp.1
yaxis <- pcamerged$Comp.2
plot(xaxis,yaxis, bg = colors[unclass(factor(categories))],main="Plot of PCA loadings labels by time",
     xlab="PC1loadings",ylab="PC2loadings",pch=21)
legend(locator(1), sort(matrix((unique(categories)))), 
       pt.bg = colors,text.col = "black", pch = 21, bg = 'white')
identify(pcamerged$Comp.1,pcamerged$Comp.2,labels=as.character(pcamerged$Row.names))
