
setwd("C:/Users/JM/thesis/mark_work/ekta_netsnps/correlations_hgnc_0deg")

## input
data <- read.table("mppVnonmpp.txt", header=TRUE, sep="\t")
d1 <- na.omit(data$mppNSsnpdensity)
d2 <- na.omit(data$nonmppNSsnpdensity)

d1name = "mpp"
d2name = "non-mpp"
xname = "NS SNP density"

colors = c("blue","red")

## plots
x11()
b1 <- boxplot(data,ylab=xname)

x11()
plot(density(d2),main="",xlab=xname,ylab="frequency",col=colors[1],lwd=3,cex.lab=1.5,cex.axis=1.2) #lwd=length width
lines(density(d1),main="",xlab=xname,ylab="frequency",col=colors[2],lwd=3,cex.lab=1.5,cex.axis=1.2)

legend(locator(1),c(d2name,d1name),col=colors,lwd=c(2,2))

## kstest
ks.test(d1,d2)