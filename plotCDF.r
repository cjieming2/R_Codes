setwd("C:/Documents and Settings/chenjm/Desktop")

# a.txt file should have only the column for call rate 
x <- read.table("snps.txt", header=TRUE)
x <- sort(x[,])
#x <- trunc(x/1000)

Fn <- ecdf(x)

plot(x, Fn(x), lwd = 2, col = "red", type = 'l', xlab="snp", ylab="snp-call-rate", main="Cumulative Distribution of reclustered-SNP-call-rate")
#plot(x, Fn(x), lwd = 2, col = "red", type = 'p', xlab="snp", ylab="snp-call-rate", main="Cumulative Distribution of reclustered-SNP-call-rate")


# this gives the percentiles of the values in x
# Fn(x)

median(x)
quantile(x, probs=c(10,25,50,75)/100, type=1, names = TRUE)
max(x)
min(x)

100*Fn(1000) 
100*Fn(5000)
100*Fn(10000)
100*Fn(20000) 
100*Fn(50000) 
100*Fn(100000) 
100*Fn(250000)
100*Fn(500000)