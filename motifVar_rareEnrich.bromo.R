setwd('C:/Users/Jieming/Documents/thesis/lynne_work/motifVar/bromo')

domain = 'bromo'

## data
data = read.table('1KG.snps.nonmono.smartDomain2gPos.BROMO.109aa.enrich', header=T, stringsAsFactors = F)
data[data==0.001] <- 0
data.mnoS = matrix(nrow = nrow(data), ncol = 2)
data.mS = matrix(nrow = nrow(data), ncol = 2)

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

### plot bar plots
# rare.ns vs comm.ns
rare.comm.ns = data$numRareNS.noS + data$numCommNS.noS
x11()
par(cex.axis=1, cex.lab=2, cex.main=2)
barplot(rare.comm.ns,names.arg=data$resNum,main=paste('variants by ', domain, ' motif pos'),col='blue')
barplot(data$numRareNS.noS,add=T,xaxt="n",col='red') # tab.var.noSingle


legend(3,0.5,c("rare.ns","comm.ns"),
       col=c("red","blue"), cex=2, pt.cex=2,
       text.col = "black", pch = 15, bg = 'white')

# rare.s vs comm.s
rare.comm.s = data$numRareS.noS + data$numCommS.noS
x11()
par(cex.axis=1, cex.lab=2, cex.main=2)
barplot(rare.comm.s,names.arg=data$resNum,main=paste('variants by ',domain,' motif pos'),col='cyan')
barplot(data$numRareNS.noS,add=T,xaxt="n",col='magenta') # tab.var.noSingle


legend(3,0.5,c("rare.s","comm.s"),
       col=c("magenta","cyan"), cex=2, pt.cex=2,
       text.col = "black", pch = 15, bg = 'white')

#############################################################
### common v rare
x11()
num.Comm.noS = data$totVar.noS  - data$numRare.noS
a = log(data$numRare.noS/num.Comm.noS)
a[a == "Inf"] <- 2.5
a[a == "-Inf"] <- -2.5
barplot(a,names.arg=data$resNum,
        main=paste('rare:common variants by ',domain,'motif pos'), 
        col='red')

### ratio of ns-s
num.S.noS = data$numRareS.noS + data$numCommS.noS
x11()
par(cex.axis=1, cex.lab=2, cex.main=2)
b=log(data$num.NS.noS/num.S.noS)
b[b == "-Inf"] <- -2.5
b[b == "Inf"] <- 2.5
barplot(b, names.arg=data$resNum,
        main=paste('ns:s variants by ',domain,'motif pos'), 
        col='orange')

### ratio of ns.rare v s.rare
x11()
par(cex.axis=1, cex.lab=2, cex.main=2)
c=log(data$numRareNS.noS/data$numRareS.noS)
c[c == "-Inf"] <- -2.5
c[c == "Inf"] <- 2.5
barplot(c, names.arg=data$resNum,
        main=paste('ns.rare:s.rare variants by ',domain,'motif pos'), 
        col='green')

### ratio of ns.comm v s.comm
x11()
par(cex.axis=1, cex.lab=2, cex.main=2)
d=log(data$numCommNS.noS/data$numCommS.noS)
d[d == "-Inf"] <- -2.5
d[d == "Inf"] <- 2.5
d[d == "NaN"] <- 0
barplot(d, names.arg=data$resNum,
        main=paste('ns.comm:s.comm variants by ',domain,'motif pos'), 
        col='black')

x11()
barplot(rbind(c,d), names.arg=data$resNum,
        main=paste('ns.rarecomm:s.rarecomm variants by ',domain,'motif pos'), 
        col=c('green','black'),ylim=c(-2.5,2.5),beside=TRUE)

legend(20,2.5,c("nsrare:srare","nscomm:scomm"),c('green','black'))
