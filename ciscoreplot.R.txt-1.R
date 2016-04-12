
#make score plot with confidence ellipse.
#arguments are output from PCA, vector with components for plotting (usually c(1,2) or c(1,3)
#and a vector of names for the points

ciscoreplot<-function(x,comps,namevec){
  y1<-sqrt(5.99*(x$sdev[comps[1]]^2))
  ymod<-y1-y1%%.05
  y1vec<-c(-y1,seq(-ymod,ymod,by=0.05),y1)
  y2vecpos<-sqrt((5.99-(y1vec^2)/x$sdev[comps[1]]^2)*x$sdev[comps[2]]^2)
  y2vecneg<--sqrt((5.99-(y1vec^2)/x$sdev[comps[1]]^2)*x$sdev[comps[2]]^2)
  y2vecpos[1]<-0
  y2vecneg[1]<-0
  y2vecpos[length(y2vecpos)]<-0
  y2vecneg[length(y2vecneg)]<-0
  plot(x$scores[,comps[1]],x$scores[,comps[2]],pch=19,cex=1.2,ylim=c(min(y2vecneg,x$scores[,comps[2]]),max(y2vecpos,x$scores[,comps[2]])),
    main="PC Score Plot with 95% CI Ellipse", xlab=paste("Scores for PC",comps[1],sep=" "), ylab=paste("Scores for PC",comps[2],sep=" "),
    xlim=c(min(y1vec,x$scores[,comps[1]]),max(y1vec,x$scores[,comps[1]])))
  lines(y1vec,y2vecpos,col="Red",lwd=2)
  lines(y1vec,y2vecneg,col="Red",lwd=2)
  outliers<-((x$scores[,comps[1]]^2)/(x$sdev[comps[1]]^2)+(x$scores[,comps[2]]^2)/(x$sdev[comps[2]]^2))>5.99
  points(x$scores[outliers,comps[1]],x$scores[outliers,comps[2]],pch=19,cex=1.2,col="Blue")
  text(x$scores[outliers,comps[1]],x$scores[outliers,comps[2]],col="Blue",lab=namevec[outliers])
}

#here is an example of how to call the function for the poverty dataset
ciscoreplot(pc1,c(1,2),povtrans[,1])

