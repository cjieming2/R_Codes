
x11()
x = c(9.375,18.75,37.5,75,150,300,600)
x1 = c(31.25,62.5,125,250,500,1000,2000)
y = c(0.022,0.039,0.079,0.148,0.287,0.4625,0.7515)

plot(x,y,col="black",pch=19,xlim=c(0,600),ylim=c(0,1))
par(new=T)
plot(x1,y,col="red",pch=19,xlim=c(0,600),ylim=c(0,1))

## fit1 linear ####
fit = lm(y~x)
summary(fit)

xx = seq(0,600,50)
par(new=T)
# plot(xx,predict(fit,data.frame(x=xx)), col="blue",type="b",lty=2,xlim=c(0,600),ylim=c(0,1))
yy = fit$coefficient[1] + fit$coefficient[2]*xx
plot(xx,yy,col="blue",type="b",lty=2,xlim=c(0,600),ylim=c(0,1))

## fit2 poly deg 2 ####
fit.poly = lm(y ~ poly(x, 2, raw=TRUE))
summary(fit.poly)

par(new=T)
# plot(xx,predict(fit.poly,data.frame(x=xx)), col="orange",type="b",lty=2,xlim=c(0,600),ylim=c(0,1))
yy.poly = fit.poly$coefficient[1] + fit.poly$coefficient[2]*xx + fit.poly$coefficient[3]*xx^2
plot(xx,yy.poly,col="orange", type="b",lty=2,xlim=c(0,600),ylim=c(0,1))

## fit3 poly deg 3 ####
fit.poly3 = lm(y ~ poly(x, 3, raw=TRUE))
summary(fit.poly3)

par(new=T)
# plot(xx,predict(fit.poly3,data.frame(x=xx)), col="orange",type="b",lty=2,xlim=c(0,600),ylim=c(0,1))
yy.poly3 = fit.poly3$coefficient[1] + fit.poly3$coefficient[2]*xx + fit.poly3$coefficient[3]*xx^2
plot(xx,yy.poly3,col="green", type="b",lty=2,xlim=c(0,600),ylim=c(0,1))
