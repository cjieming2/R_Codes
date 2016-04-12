#load libraries for analysis
library(vegan)
library(akima)
library(mgcv)

#get the data.
alaska <- read.delim("http://reuningscherer.net/stat660/data/streamdata.txt",sep="\t")

#1
#several sites had NO species. Delete these.
alaska<-alaska[apply(alaska[,4:14],1,sum)>0,]
#assign names to rows
rownames(alaska)<-alaska$StreamID
#make species count matrix
alaskacount<-alaska[,4:14]
#make environmental data matrix
alaskaenv<-alaska[,15:20]
#perform correspondence analysis
alaskaca<-cca(alaskacount)
x11()
plot(alaskaca,main="Correspondence Analysis for Alaska Stream Data")
#add environmental variables
fit<-envfit(alaskaca, alaskaenv, permutations=1000)
fit
#take out insignificant env vars.
fit<-envfit(alaskaca, alaskaenv [,-c(2,3,5,6)], permutations=1000)
plot(fit) # plots arrow on current plot
fit

#2
#perform detrended correspondence analysis
alaskadca <-decorana(alaskacount)
x11()
plot(alaskadca,main="Detrended CA for Alaska Stream Data")
#add environmental variables
fit<-envfit(alaskadca, alaskaenv, permutations=1000)
fit
#take out insignificant env vars.
fit<-envfit(alaskadca, alaskaenv [,-c(2)], permutations=1000)
plot(fit)
fit

#3
#perform nonmetric multidimensional scaling (NMDS)
x11()
alaskanmds<-metaMDS(alaskacount)
fig<-ordiplot(alaskanmds,type="none",cex=1.1,main="NMDS for Alaska Stream Data")
text(fig,"species",col="red",cex=1.1)
text(fig,"sites",col="black",cex=0.8)
#add environmental variables
fit<-envfit(alaskanmds, alaskaenv, permutations=1000)
fit
#take out insignificant env vars.
fit<-envfit(alaskadca, alaskaenv [,-c(2,3,5,6)], permutations=1000)
plot(fit)
fit

#4
#perform canonical correspondence analysis
alaskacca<-cca(alaskacount, alaskaenv[,-c(2,3,5,6)], main="Canonical CA for Alaska Stream Data")
x11()
plot(alaskacca)
alaskacca

#5
#make interpolated surface on nmds
x11()
alaskanmds<-metaMDS(alaskacount)
fit1<-envfit(alaskanmds ~temp+age, alaskaenv)
#first plot with interpolated surface
fig<-ordiplot(alaskanmds,type="none",cex=1.1,main="NMDS for Alaska Stream Data")
text(fig,"species",col="red",cex=0.7)
text(fig,"sites",col="black",cex=0.7)
plot(fit1)
tmp1 <- with(alaskaenv, ordisurf(alaskanmds, age, add = TRUE))
tmp2<-with(alaskaenv, ordisurf(alaskanmds, temp, add = TRUE, col = "green4"))