setwd("C:/Users/JM/thesis/mark_work/allele_specificity/sharing")

## data
# asb.int.AFR = read.table('mt2AS.asb.maf.pop.AFR.bed',stringsAsFactors = F)
# asb.int.ASI = read.table('mt2AS.asb.maf.pop.ASI.bed',stringsAsFactors = F)
# asb.int.CEU = read.table('mt2AS.asb.maf.pop.CEU.bed',stringsAsFactors = F)
# asb.int.CHB = read.table('mt2AS.asb.maf.pop.CHB.bed',stringsAsFactors = F)
# asb.int.EUR = read.table('mt2AS.asb.maf.pop.EUR.bed',stringsAsFactors = F)
# asb.int.FIN = read.table('mt2AS.asb.maf.pop.FIN.bed',stringsAsFactors = F)
# asb.int.GBR = read.table('mt2AS.asb.maf.pop.GBR.bed',stringsAsFactors = F)
# asb.int.JPT = read.table('mt2AS.asb.maf.pop.JPT.bed',stringsAsFactors = F)
# asb.int.TSI = read.table('mt2AS.asb.maf.pop.TSI.bed',stringsAsFactors = F)
# asb.int.YRI = read.table('mt2AS.asb.maf.pop.YRI.bed',stringsAsFactors = F)
# ase.int.AFR = read.table('mt2AS.ase.maf.pop.AFR.bed',stringsAsFactors = F)
# ase.int.ASI = read.table('mt2AS.ase.maf.pop.ASI.bed',stringsAsFactors = F)
# ase.int.CEU = read.table('mt2AS.ase.maf.pop.CEU.bed',stringsAsFactors = F)
# ase.int.CHB = read.table('mt2AS.ase.maf.pop.CHB.bed',stringsAsFactors = F)
# ase.int.EUR = read.table('mt2AS.ase.maf.pop.EUR.bed',stringsAsFactors = F)
# ase.int.FIN = read.table('mt2AS.ase.maf.pop.FIN.bed',stringsAsFactors = F)
# ase.int.GBR = read.table('mt2AS.ase.maf.pop.GBR.bed',stringsAsFactors = F)
# ase.int.JPT = read.table('mt2AS.ase.maf.pop.JPT.bed',stringsAsFactors = F)
# ase.int.TSI = read.table('mt2AS.ase.maf.pop.TSI.bed',stringsAsFactors = F)
# ase.int.YRI = read.table('mt2AS.ase.maf.pop.YRI.bed',stringsAsFactors = F)
asb.int.AFR = read.table('mt2AS.asb.maf.cou.pop.AFR.bed',stringsAsFactors = F)
asb.int.ASI = read.table('mt2AS.asb.maf.cou.pop.ASI.bed',stringsAsFactors = F)
asb.int.CEU = read.table('mt2AS.asb.maf.cou.pop.CEU.bed',stringsAsFactors = F)
asb.int.CHB = read.table('mt2AS.asb.maf.cou.pop.CHB.bed',stringsAsFactors = F)
asb.int.EUR = read.table('mt2AS.asb.maf.cou.pop.EUR.bed',stringsAsFactors = F)
asb.int.FIN = read.table('mt2AS.asb.maf.cou.pop.FIN.bed',stringsAsFactors = F)
asb.int.GBR = read.table('mt2AS.asb.maf.cou.pop.GBR.bed',stringsAsFactors = F)
asb.int.JPT = read.table('mt2AS.asb.maf.cou.pop.JPT.bed',stringsAsFactors = F)
asb.int.TSI = read.table('mt2AS.asb.maf.cou.pop.TSI.bed',stringsAsFactors = F)
asb.int.YRI = read.table('mt2AS.asb.maf.cou.pop.YRI.bed',stringsAsFactors = F)
ase.int.AFR = read.table('mt2AS.ase.maf.cou.pop.AFR.bed',stringsAsFactors = F)
ase.int.ASI = read.table('mt2AS.ase.maf.cou.pop.ASI.bed',stringsAsFactors = F)
ase.int.CEU = read.table('mt2AS.ase.maf.cou.pop.CEU.bed',stringsAsFactors = F)
ase.int.CHB = read.table('mt2AS.ase.maf.cou.pop.CHB.bed',stringsAsFactors = F)
ase.int.EUR = read.table('mt2AS.ase.maf.cou.pop.EUR.bed',stringsAsFactors = F)
ase.int.FIN = read.table('mt2AS.ase.maf.cou.pop.FIN.bed',stringsAsFactors = F)
ase.int.GBR = read.table('mt2AS.ase.maf.cou.pop.GBR.bed',stringsAsFactors = F)
ase.int.JPT = read.table('mt2AS.ase.maf.cou.pop.JPT.bed',stringsAsFactors = F)
ase.int.TSI = read.table('mt2AS.ase.maf.cou.pop.TSI.bed',stringsAsFactors = F)
ase.int.YRI = read.table('mt2AS.ase.maf.cou.pop.YRI.bed',stringsAsFactors = F)


## data with MAF<0.01 (rare)
maf=0.01
asb.int.hist.AFR <- table(asb.int.AFR$V11[which(asb.int.AFR$V5<=maf)]); asb.int.hist.AFR <- prop.table(asb.int.hist.AFR)
asb.int.hist.ASI <- table(asb.int.ASI$V11[which(asb.int.ASI$V5<=maf)]); asb.int.hist.ASI <- prop.table(asb.int.hist.ASI)
asb.int.hist.EUR <- table(asb.int.EUR$V11[which(asb.int.EUR$V5<=maf)]); asb.int.hist.EUR <- prop.table(asb.int.hist.EUR)

asb.int.hist.CEU <- table(asb.int.CEU$V9[which(asb.int.CEU$V5<=maf)]); asb.int.hist.CEU <- prop.table(asb.int.hist.CEU)
asb.int.hist.CHB <- table(asb.int.CHB$V9[which(asb.int.CHB$V5<=maf)]); asb.int.hist.CHB <- prop.table(asb.int.hist.CHB)
asb.int.hist.FIN <- table(asb.int.FIN$V9[which(asb.int.FIN$V5<=maf)]); asb.int.hist.FIN <- prop.table(asb.int.hist.FIN)
asb.int.hist.GBR <- table(asb.int.GBR$V9[which(asb.int.GBR$V5<=maf)]); asb.int.hist.GBR <- prop.table(asb.int.hist.GBR)
asb.int.hist.JPT <- table(asb.int.JPT$V9[which(asb.int.JPT$V5<=maf)]); asb.int.hist.JPT <- prop.table(asb.int.hist.JPT)
asb.int.hist.TSI <- table(asb.int.TSI$V9[which(asb.int.TSI$V5<=maf)]); asb.int.hist.TSI <- prop.table(asb.int.hist.TSI)
asb.int.hist.YRI <- table(asb.int.YRI$V9[which(asb.int.YRI$V5<=maf)]); asb.int.hist.YRI <- prop.table(asb.int.hist.YRI)

ase.int.hist.AFR <- table(ase.int.AFR$V11[which(ase.int.AFR$V5<=maf)]); ase.int.hist.AFR <- prop.table(ase.int.hist.AFR)
ase.int.hist.ASI <- table(ase.int.ASI$V11[which(ase.int.ASI$V5<=maf)]); ase.int.hist.ASI <- prop.table(ase.int.hist.ASI)
ase.int.hist.EUR <- table(ase.int.EUR$V11[which(ase.int.EUR$V5<=maf)]); ase.int.hist.EUR <- prop.table(ase.int.hist.EUR)

ase.int.hist.CEU <- table(ase.int.CEU$V9[which(ase.int.CEU$V5<=maf)]); ase.int.hist.CEU <- prop.table(ase.int.hist.CEU)
ase.int.hist.CHB <- table(ase.int.CHB$V9[which(ase.int.CHB$V5<=maf)]); ase.int.hist.CHB <- prop.table(ase.int.hist.CHB)
ase.int.hist.FIN <- table(ase.int.FIN$V9[which(ase.int.FIN$V5<=maf)]); ase.int.hist.FIN <- prop.table(ase.int.hist.FIN)
ase.int.hist.GBR <- table(ase.int.GBR$V9[which(ase.int.GBR$V5<=maf)]); ase.int.hist.GBR <- prop.table(ase.int.hist.GBR)
ase.int.hist.JPT <- table(ase.int.JPT$V9[which(ase.int.JPT$V5<=maf)]); ase.int.hist.JPT <- prop.table(ase.int.hist.JPT)
ase.int.hist.TSI <- table(ase.int.TSI$V9[which(ase.int.TSI$V5<=maf)]); ase.int.hist.TSI <- prop.table(ase.int.hist.TSI)
ase.int.hist.YRI <- table(ase.int.YRI$V9[which(ase.int.YRI$V5<=maf)]); ase.int.hist.YRI <- prop.table(ase.int.hist.YRI) ## works





## plot col9 #uniqpop 7 pop CEU CHB JPT YRI TSI GBR FIN
# ASE
x11()
par(mar=c(5,8,4,4),xpd=TRUE)
xlimit=5;
ylimit=0.9;

# par(cex.axis=2, cex.lab=1, cex.main=1.2, cex.sub=1)
par(cex.axis=2, cex.lab=2)
plot(ase.int.hist.GBR, col="lightseagreen", xlim=c(1,xlimit), ylim=c(0,ylimit), type="l", lty=1, xlab='number of unique populations', ylab='proportion of rare SNVs that are shared', main='ASE')
par(new=TRUE); plot(ase.int.hist.FIN, type="l", lty=1, col="skyblue", xlim=c(1,xlimit), ylim=c(0,ylimit), axes=FALSE, xlab='',ylab='')
par(new=TRUE); plot(ase.int.hist.CEU, type="l", lty=1, col="blue", xlim=c(1,xlimit), ylim=c(0,ylimit), axes=FALSE, xlab='',ylab='')
par(new=TRUE); plot(ase.int.hist.TSI, type="l", lty=1, col="cyan", xlim=c(1,xlimit), ylim=c(0,ylimit), axes=FALSE, xlab='',ylab='')
par(new=TRUE); plot(ase.int.hist.CHB, type="l", lty=1, col="green", xlim=c(1,xlimit), ylim=c(0,ylimit), axes=FALSE, xlab='',ylab='')
par(new=TRUE); plot(ase.int.hist.JPT, type="l", lty=1, col="darkgreen", xlim=c(1,xlimit), ylim=c(0,ylimit), axes=FALSE, xlab='',ylab='')
par(new=TRUE); plot(ase.int.hist.YRI, type="l", lty=1, col="red", xlim=c(1,xlimit), ylim=c(0,ylimit), axes=FALSE, xlab='',ylab='')

## legend
legend(4,0.8,c('CEU','GBR','FIN','TSI','CHB','JPT','YRI'), 
       col=c('blue','lightseagreen','skyblue','cyan','green','darkgreen','red'), 
       text.col = "black", pch = c(16,16,16,16,16,16,16), bg = 'white',cex=2)

# ASB
x11()
par(mar=c(5,8,4,4),xpd=TRUE)
xlimit=2;
ylimit=1;

par(cex.axis=2, cex.lab=2);
               plot(asb.int.hist.CEU, col="blue", xlim=c(1,xlimit), ylim=c(0,ylimit), type="l", lty=1, xlab='number of unique populations', ylab='proportion of rare SNVs that are shared',main='ASB')
par(new=TRUE); plot(asb.int.hist.CHB, type="l", lty=1, col="green", xlim=c(1,xlimit), ylim=c(0,ylimit), axes=FALSE, xlab='',ylab='')
# par(new=TRUE); plot(asb.int.hist.FIN, type="l", lty=1, col="skyblue", xlim=c(1,xlimit), ylim=c(0,ylimit), axes=FALSE, xlab='',ylab='')
# par(new=TRUE); plot(asb.int.hist.GBR, type="l", lty=1, col="lightseagreen", xlim=c(1,xlimit), ylim=c(0,ylimit), axes=FALSE, xlab='',ylab='')
par(new=TRUE); plot(asb.int.hist.JPT, type="l", lty=1, col="darkgreen", xlim=c(1,xlimit), ylim=c(0,ylimit), axes=FALSE, xlab='',ylab='')
# par(new=TRUE); plot(asb.int.hist.TSI, type="l", lty=1, col="cyan", xlim=c(1,xlimit), ylim=c(0,ylimit), axes=FALSE, xlab='',ylab='')
par(new=TRUE); plot(asb.int.hist.YRI, type="l", lty=1, col="red", xlim=c(1,xlimit), ylim=c(0,ylimit), axes=FALSE, xlab='',ylab='')

## legend
legend(1.5,0.8,c('CEU','CHB','JPT','YRI'), 
       col=c('blue','green','darkgreen','red'), 
       text.col = "black", pch = c(16,16,16,16,16,16,16), bg = 'white',cex=2)


##################################################################
## world pop
##################################################################

## plot col11 EUR ASI AFR
# ASE
x11()
par(mar=c(5,8,4,4),xpd=TRUE)
xlimit=3;
ylimit=1;

# par(cex.axis=2, cex.lab=1, cex.main=1.2, cex.sub=1)
par(cex.axis=2, cex.lab=2)
               plot(ase.int.hist.AFR, col="red", xlim=c(1,xlimit), ylim=c(0,ylimit), type="l", lty=1, xlab='number of unique populations', ylab='proportion of rare SNVs that are shared',main='ASE')
par(new=TRUE); plot(ase.int.hist.ASI, type="l", lty=1, col="green", xlim=c(1,xlimit), ylim=c(0,ylimit), axes=FALSE, xlab='',ylab='')
par(new=TRUE); plot(ase.int.hist.EUR, type="l", lty=1, col="blue", xlim=c(1,xlimit), ylim=c(0,ylimit), axes=FALSE, xlab='',ylab='')

## legend
legend(2.5,0.8,c('EUR','ASI','AFR'), 
       col=c('blue','green','red'), 
       text.col = "black", pch = c(16,16,16,16,16,16,16), bg = 'white',cex=2)

# ASB
x11()
par(mar=c(5,8,4,4),xpd=TRUE)
xlimit=2;
ylimit=1;

# par(cex.axis=2, cex.lab=1, cex.main=1.2, cex.sub=1)
par(cex.axis=2, cex.lab=2)
               plot(asb.int.hist.AFR, col="red", xlim=c(1,xlimit), ylim=c(0,ylimit), type="l", lty=1, xlab='number of unique populations', ylab='proportion of rare SNVs that are shared', main='ASB')
par(new=TRUE); plot(asb.int.hist.ASI, type="l", lty=1, col="green", xlim=c(1,xlimit), ylim=c(0,ylimit), axes=FALSE, xlab='',ylab='')
par(new=TRUE); plot(asb.int.hist.EUR, type="l", lty=1, col="blue", xlim=c(1,xlimit), ylim=c(0,ylimit), axes=FALSE, xlab='',ylab='')

## legend
legend(1.5,0.8,c('EUR','ASI','AFR'), 
       col=c('red','green','blue'), 
       text.col = "black", pch = c(16,16,16,16,16,16,16), bg = 'white',cex=2)