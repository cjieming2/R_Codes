setwd('C:/Users/Jieming/Documents/thesis/dissertation/thesis_defense/figures')

library(ggplot2)
library(gridExtra)

# file = 'sequencing_costs_apr2015.jmedited.txt'
file = 'sequencing_costs_apr2015.jmedited.moore.txt'
data = read.table(file,header=T,stringsAsFactors=F)

## from 2001
# miny=3
miny=3
maxy=8
x11()
par(mar=c(5,5,1,1),cex.axis=1, cex.lab=2, cex.main=2)
plot(seq(1,nrow(data)),log10(data[,3]), xaxt="n",xlab="", ylab ="", col=c("black"),yaxt="n",
     ylim=c(miny,maxy), type='l',lwd=6)
grid(nx="",ny=nrow(seq(miny,maxy)),col="grey", lty=1)
par(new=T)
plot(seq(1,nrow(data)),log10(data[,3]), xaxt="n",xlab="", ylab ="", col=c("black"),yaxt="n",
     ylim=c(miny,maxy), type='l',lwd=6)
axis(1,at=seq(1,nrow(data)),label=data[,1],las=2)
# axis(2,at=seq(miny,maxy),label=c("$1K","$10K","$100K","$1M","$10M","$100M","$1B","$10B"),las=2)
axis(2,at=seq(miny,maxy),label=c("$1K","$10K","$100K","$1M","$10M","$100M"),las=2)
mtext("Cost per genome", side=2,line=3.5,cex=2)

## moore's law
par(new=T)
# plot(seq(1,nrow(data)),na.approx(log10(data[,4]),na.rm=F), xaxt="n",xlab="", ylab ="", col=c("grey"),yaxt="n",
#      ylim=c(miny,maxy), type='l',lwd=6)

plot(seq(1,nrow(data)),log10(data[,4]), xaxt="n",xlab="", ylab ="", col=c("grey"),yaxt="n",
     ylim=c(miny,maxy),lwd=6)

## until 2003 (zoomed)
row=33
# x11()
par(mar=c(5,5,1,1),cex.axis=1, cex.lab=2, cex.main=2)
plot(seq(1,row),log10(data[1:33,3]), xaxt="n",xlab="", ylab ="", col=c("black"),yaxt="n",
     ylim=c(miny,maxy), type='l',lwd=6)
grid(nx="",ny=nrow(seq(miny,maxy)),col="grey", lty=1)
par(new=T)
plot(seq(1,row),log10(data[1:33,3]), xaxt="n",xlab="", ylab ="", col=c("black"),yaxt="n",
     ylim=c(miny,maxy), type='l',lwd=6)
axis(1,at=seq(1,nrow(data)),label=data[,1],las=2)
axis(2,at=seq(miny,maxy),label=c("$1K","$10K","$100K","$1M","$10M","$100M"),las=2)
mtext("Cost per genome", side=2,line=3.5,cex=2)


## showing the increase in number of exomes/genomes = individual's DNA sequence
year=c("2007","2008","2009","2010","2011","2012","2013","2014","2015")
x11()
