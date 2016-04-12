### this clipplot requires the TeachingDemos package
### clipplot(fun, xlim = par("usr")[1:2], ylim = par("usr")[3:4])
### threshold = 0


x <- rnorm(1:100)
x11()
plot(x, type='l',col='red')
clipplot( lines(x, col='blue'), ylim = c(min(x), 0) )

#################################################################
## same as above

x <- seq(1,100)
y <- rnorm(100)
plot(x,y, type='b', col='blue')
clipplot( lines(x,y, type='b', col='red'), ylim=c(par('usr')[3],0))

#################################################################

attach(iris)

tmp <- c('red','green','blue')
names(tmp) <- levels(Species)
x11()
plot(Petal.Width,Petal.Length, col=tmp[Species])
for(s in levels(Species)){
  clipplot( abline(
    lm(Petal.Length~Petal.Width, data=iris, subset=Species==s),
    col=tmp[s]),
    xlim=range(Petal.Width[Species==s]))
}

detach(iris)





