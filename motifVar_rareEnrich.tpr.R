setwd('C:/Users/Jieming/Documents/thesis/lynne_work/motifVar/tpr/exac')
# setwd('C:/Users/Jieming/Documents/thesis/lynne_work/motifVar/tpr/1kg.esp6500')
# setwd('C:/Users/Jieming/Documents/thesis/lynne_work/motifVar/tpr/3tpr')
# setwd('C:/Users/Jieming/Documents/thesis/lynne_work/motifVar/tpr/exac')
# setwd('C:/Users/Jieming/Documents/thesis/lynne_work/motifVar/tpr/single proteins/TTC21B')
# setwd('C:/Users/Jieming/Documents/thesis/lynne_work/motifVar/tpr/1kg.esp6500')

domain = 'ExAC TPRs'

## data
# data = read.table('combined.1KG.ESP6500.snps.nonmono.smartDomain2gPos.TPR.34aa.noS.enrich', sep="\t", header=T, stringsAsFactors = F)
# data = read.table('1KG.snps.nonmono.smartDomain2gPos.TPR.34aa.noS.enrich', sep="\t", header=T, stringsAsFactors = F)
# data = read.table('final.merged.ExAC.ens73.TPR.34aa.vepoutput.coding.canonical.sorted.auto.571.uniqseq.1tpr.AC.enrich', sep="\t", header=T, stringsAsFactors = F)
data = read.table('final.merged.ExAC.ens73.TPR.34aa.vepoutput.coding.canonical.sorted.auto.enrich', sep="\t", header=T, stringsAsFactors = F)
data[data==0.001] <- 0
data.mnoS = matrix(nrow = nrow(data), ncol = 2)
data.mS = matrix(nrow = nrow(data), ncol = 2)
order = read.table('order_relentropy-tpr.txt', header=T, stringsAsFactors = F, comment.char="", sep="\t")
data$V20 = factor(data$resNum, levels=levels(data$resNum)<-order$pos)

## fisher's test
## pval and foldchange in data.m
## S and noS are with singletons (AC=1) and without singletons
## each row in the matrix is the position
for(i in 1:nrow(data))
{
  d.noS <- matrix(c(data$numRareNS.noS[i],data$numCommNS.noS[i],data$numRareS.noS[i],
                          data$numCommS.noS[i]), 2,2, 
                   dimnames = list(Pathways = c("rare","comm"),SNPs = c("NS","S")))  
  
  d.S <- matrix(c(data$numRareNS[i],data$numCommNS[i],data$numRareS[i],
                data$numCommS[i]), 2,2, 
              dimnames = list(Pathways = c("rare","comm"),SNPs = c("NS","S")))  
  
  x.noS = fisher.test(d.noS,alternative="two.sided")
  x.S = fisher.test(d.S,alternative="two.sided")
  
  data.mnoS[i,1] = x.noS$p.value
  data.mnoS[i,2] = x.noS$estimate
  
  data.mS[i,1] = x.S$p.value
  data.mS[i,2] = x.S$estimate
}


############################################ plot bar plots 
###### rare.ns vs comm.ns ######
#### noS
rare.comm.ns.noS = data$numRareNS.noS + data$numCommNS.noS
yuplimit = max(rare.comm.ns.noS)+5
x11()
par(cex.axis=1, cex.lab=2, cex.main=2)
# barplot(rare.comm.ns,names.arg=data$resNum,main=paste('variants by ', domain, ' motif pos'),col='blue')
# barplot(data$numRareNS.noS,add=T,xaxt="n",col='red') # tab.var.noSingle
barplot(rare.comm.ns.noS,names.arg=data$resNum,main=paste('variants by ', domain, ' motif pos'),
        col='blue', ylim=c(0,yuplimit))
barplot(data$numRareNS.noS,add=T,xaxt="n",col='red',ylim=c(0,yuplimit)) # tab.var.noSingle

legend(3,yuplimit,c("rare.ns.noS","comm.ns.noS"),
       col=c("red","blue"), cex=2, pt.cex=2,
       text.col = "black", pch = 15, bg = 'white')

#### S
# rare.comm.ns = data$numRareNS + data$numCommNS
# yuplimit = 250
# x11()
# par(cex.axis=1, cex.lab=2, cex.main=2)
# barplot(rare.comm.ns,names.arg=data$resNum,main=paste('variants by ', domain, ' motif pos'),
#         col='blue', ylim=c(0,yuplimit))
# barplot(data$numRareNS,add=T,xaxt="n",col='red',ylim=c(0,yuplimit)) # tab.var
# 
# legend(3,yuplimit,c("rare.ns","comm.ns"),
#        col=c("red","blue"), cex=2, pt.cex=2,
#        text.col = "black", pch = 15, bg = 'white')


############################################ 
###### rare.s vs comm.s ####
#### noS
rare.comm.s.noS = data$numRareS.noS + data$numCommS.noS
yuplimit = max(rare.comm.s.noS)+5
x11()
par(cex.axis=1, cex.lab=2, cex.main=2)
barplot(rare.comm.s.noS,names.arg=data$resNum,main=paste('variants by ',domain,' motif pos'),
        col='cyan',ylim=c(0,yuplimit))
barplot(data$numRareS.noS,add=T,xaxt="n",col='magenta',ylim=c(0,yuplimit)) # tab.var.noSingle


legend(3,yuplimit,c("rare.s.noS","comm.s.noS"),
       col=c("magenta","cyan"), cex=2, pt.cex=2,
       text.col = "black", pch = 15, bg = 'white')

#### S
# rare.comm.s = data$numRareS + data$numCommS
# yuplimit = 150
# x11()
# par(cex.axis=1, cex.lab=2, cex.main=2)
# barplot(rare.comm.s,names.arg=data$resNum,main=paste('variants by ',domain,' motif pos'),
#         col='cyan',ylim=c(0,yuplimit))
# barplot(data$numRareS,add=T,xaxt="n",col='magenta',ylim=c(0,yuplimit)) # tab.var
# 
# legend(3,yuplimit,c("rare.s","comm.s"),
#        col=c("magenta","cyan"), cex=2, pt.cex=2,
#        text.col = "black", pch = 15, bg = 'white')

#############################################################
### common v rare
## noS
num.Comm.noS = data$totVar.noS  - data$numRare.noS
a = log(data$numRare.noS/num.Comm.noS)
a.finite = a[!is.infinite(a)]
a[a == "Inf"] <- max(a.finite) +1
a[a == "-Inf"] <- -max(a.finite) -1
x11()
barplot(a,names.arg=data$resNum,
        main=paste('rare:common variants noS by ',domain,'motif pos'), 
        col='red')

## S
# x11()
# num.Comm = data$totVar  - data$numRare
# a = log(data$numRare/num.Comm)
# a[a == "Inf"] <- 2.5
# a[a == "-Inf"] <- -2.5
# barplot(a,names.arg=data$resNum,
#         main=paste('rare:common variants by ',domain,'motif pos'), 
#         col='red')

#############################################################
### ratio of ns-s
## noS
num.S.noS = data$numRareS.noS + data$numCommS.noS
x11()
par(cex.axis=1, cex.lab=2, cex.main=2)
b=log(data$num.NS.noS/num.S.noS)
b.finite = b[!is.infinite(a)]
b[b == "Inf"] <- max(a.finite) +1
b[b == "-Inf"] <- -max(a.finite) -1
barplot(b, names.arg=data$resNum,
        main=paste('ns:s variants noS by ',domain,'motif pos'), 
        col='orange',ylim=c(-2,2))

# ordered
data.ordered = data[order(data$V20),]
num.S.noS.ordered = data.ordered$numRareS.noS + data.ordered$numCommS.noS
x11()
b1=log(data.ordered$num.NS.noS/num.S.noS.ordered)
b1.finite = b1[!is.infinite(b1)]
b1[b1 == "Inf"] <- max(b1.finite) +1
b1[b1 == "-Inf"] <- -max(b1.finite) -1
barplot(b1, names.arg=data.ordered$resNum,
        main=paste('ns:s variants noS by ',domain,'motif pos'), 
        col='orange',ylim=c(-2,2))

## S
# num.S = data$numRareS + data$numCommS
# x11()
# par(cex.axis=1, cex.lab=2, cex.main=2)
# b=log(data$num.NS/num.S)
# b[b == "-Inf"] <- -2.5
# b[b == "Inf"] <- 2.5
# barplot(b, names.arg=data$resNum,
#         main=paste('ns:s variants by ',domain,'motif pos'), 
#         col='orange')

#############################################################
### ratio of ns.rare v s.rare
## noS
x11()
par(cex.axis=1, cex.lab=2, cex.main=2)
c=log(data$numRareNS.noS/data$numRareS.noS)
c.finite = c[!is.infinite(c)]
c[c == "Inf"] <- max(c.finite) +1
c[c == "-Inf"] <- -max(c.finite) -1
barplot(c, names.arg=data$resNum,
        main=paste('ns.rare:s.rare variants noS by ',domain,'motif pos'), 
        col='red')

## S
# x11()
# par(cex.axis=1, cex.lab=2, cex.main=2)
# c=log(data$numRareNS/data$numRareS)
# c[c == "-Inf"] <- -2.5
# c[c == "Inf"] <- 2.5
# barplot(c, names.arg=data$resNum,
#         main=paste('ns.rare:s.rare variants by ',domain,'motif pos'), 
#         col='red')

#############################################################
### ratio of ns.comm v s.comm
## noS
x11()
par(cex.axis=1, cex.lab=2, cex.main=2)
d=log(data$numCommNS.noS/data$numCommS.noS)
d.finite = d[!is.infinite(d)]
d[d == "Inf"] <- max(d.finite) +1
d[d == "-Inf"] <- -max(d.finite) -1
# d[d == "NaN"] <- 0
barplot(d, names.arg=data$resNum,
        main=paste('ns.comm:s.comm variants noS by ',domain,'motif pos'), 
        col='black')

#############################################################
#### combining ns.rarecomm:s.rarecomm
### noS
x11()
barplot(rbind(c,d), names.arg=data$resNum,
        main=paste('ns.rarecomm:s.rarecomm variants noS by ',domain,'motif pos'), 
        col=c('red','blue'),beside=TRUE)

legend(25,2.5,c("nsrare:srare","nscomm:scomm"),c('red','blue'))


#############################################################
#### nsrare:nscommon
### noS
x11()
par(cex.axis=1, cex.lab=2, cex.main=2)
z=log(data$numRareNS.noS/data$numCommNS.noS)
z.finite = z[!is.infinite(z)]
z[z == "Inf"] <- max(z.finite) +1
z[z == "-Inf"] <- -max(z.finite) -1
yuplimit = max(z.finite)+1
barplot(z, names.arg=data$resNum,ylim=c(0,yuplimit),
        main=paste('ns.rare:ns.comm variants noS by ',domain,'motif pos'), 
        col='purple')

#############################################################
#### odds ratio = nsrare:nscommon/srare:scommon
### noS
x11()
par(cex.axis=1, cex.lab=2, cex.main=2)
numer = data$numRareNS.noS/data$numCommNS.noS
denom = data$numRareS.noS/data$numCommS.noS
y=log(numer/denom)
y.finite = y[!is.infinite(y)]
y[y == "Inf"] <- max(y.finite) +1
y[y == "-Inf"] <- -max(y.finite) -1
yuplimit = max(y.finite) +1
barplot(y, names.arg=data$resNum,ylim=c(-yuplimit,yuplimit),
        main=paste('ns.rare:ns.comm/s.rare:s.comm variants noS by ',domain,'motif pos'), 
        col='purple')

#############################################################
#### relative risk = nsrare:ns/srare:s
### noS
x11()
par(cex.axis=1, cex.lab=2, cex.main=2)
numer1 = data$numRareNS.noS/data$num.NS.noS
denom1 = data$numRareS.noS/(data$totVar.noS-data$num.NS.noS)
j=log(numer1/denom1)
barplot(j, names.arg=data$resNum,
        main=paste('ns.rare:ns/s.rare:s variants noS by ',domain,'motif pos'), 
        col='pink')
