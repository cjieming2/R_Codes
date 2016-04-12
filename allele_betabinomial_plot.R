library(VGAM)

colors <- c("red","blue","green","orange","cyan","pink","purple",
            "brown","black","slategray1","violetred","tan","deeppink","darkgreen", 
            "orchid","darksalmon","antiquewhite3","magenta","darkblue","peru","slateblue",
            "thistle","tomato","rosybrown1","royalblue","olivedrab") ##Set of 26 colours (no red) to use for all plots

k = seq(0,1,by=0.1)
n = 100
a = apply(as.data.frame(k),1,function(x) dbetabinom(seq(0,n),n,0.5,x))

leg = matrix(0,length(k))

x11()
for(i in 1:length(k))
{
  if(i==1)
  {
    plot(a[,i],ylim=c(0,0.25),type='b',col=colors[i])
    leg[i] = paste("b=",k[i])
  }
  else
  {
    par(new=TRUE)
    plot(a[,i], col=colors[i],ylim=c(0.,0.25), type='b')
    leg[i] = paste("b=",k[i])
  }
}

legend(2,0.25,leg,colors[1:length(k)])
