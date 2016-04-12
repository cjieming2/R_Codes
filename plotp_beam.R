setwd("C:/Users/JM/Desktop")

# datafile3 = "lnp-combined8.txt"
# datafile2 = "lnp.txt"
# datafile1 = "lnp.txt_"
datafile1 = "combined8-lnp.txt"
# data3 <- read.table(datafile3, header=TRUE)
# data2 <- read.table(datafile2, header=TRUE)
data1 <- read.table(datafile1, header=TRUE)

x11()
# par(mfrow=c(3,1))
# plot(data1[,1],data1[,2], 
#      xlab="iterations",ylab="ln p", main="ln p for 1 MCMC chain (diagnostic)", cex.lab=1.5, cex.main=1.5,
#      ylim=c(-4000000,-4400000))

# plot(data2[,1],data2[,2], 
#      xlab="iterations",ylab="ln p", main="ln p for 250 MCMC chains", cex.lab=1.5, cex.main=1.5,
#      ylim=c(-3959950,-3959450))
# 
plot(data1[,1],data1[,2], 
     xlab="iterations",ylab="ln p", main="ln p for 2000 MCMC chains", cex.lab=1.5, cex.main=1.5,
     ylim=c(-4000000,-4400000))

# hist(data$SymPval)