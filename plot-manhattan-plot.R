setwd("C:/Users/JM/Desktop")

library(gap)

data = read.table("manplot",header=T) ## includes snp-id, chromosome, position
color=c("red2","green","orange1","royalblue","yellow3","darkslategrey","purple3","turquoise","hotpink","lightgreen","salmon","skyblue1","goldenrod","slategrey","purple1","maroon","darkgreen","orange4","darkblue","brown","gray","mediumpurple")
x11()
# change parameters of x axis
par(las=2, xpd=TRUE, cex.axis=1.8, cex=1)

# apparently mhtplot function only changes x- axis parameters
mhtplot(data,control=mht.control(logscale=FALSE,colors=color,cex=1,
                                 xline=1.5,yline=1.5),
        pch=19,ylab="average posterior")

# change y-axis
axis(2,pos=-10,cex.axis=1)

# draw hline
abline(h=0.99,col='black',lty=3,intercept=FALSE)

# title
title("Manhattan plot for average marginal posterior from 1437 GWAS SNPs")





