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

#An example of usage of this function
#reschisqplot(danielaaov$residuals,label="Daniela MANOVA Residuals")







