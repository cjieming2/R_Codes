setwd('C:/Users/Jieming/Documents/thesis/mark_work/allele_specificity/datasets')
# dbinom(5,5,0.5) *2
# dbinom(4,5,0.5) *2
# dbinom(3,5,0.5) *2
# dbinom(2,5,0.5) *2
# dbinom(1,5,0.5) *2
# dbinom(0,5,0.5) *2
x11();par(mar=c(5,5,4,4),xpd=TRUE)
plot(seq(0,5),dbinom(seq(0,5),5,0.5),type="l",xlab='#success',ylab='pdf',cex.lab=2, cex.axis=2)


n = 100
p = 0.5
k = seq(0,n)
m = numeric(length(k))

for(i in 0:n) ## variable n
{
  m[i+1] = min(dbinom(i,i,p)*2,1)
}

plot(k,m, type="l")
write.table(cbind(seq(0,n),m),"binomialChart.txt",
            sep='\t',col.names=F,row.names=F)

