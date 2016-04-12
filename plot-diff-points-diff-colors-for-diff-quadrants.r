x <- rnorm(50)
y <- rnorm(50)

color <- rep("blue",50)
color[(x>0)&(y<0)] <-"red"
color[(x<0)&(y>0)] <-"green"
color[(x<0)&(y<0)] <-"yellow"


plot(x,y,col=color,pch=19)
abline(h=0) #just to check
abline(v=0) #just to check 