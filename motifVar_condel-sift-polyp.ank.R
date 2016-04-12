setwd("C:/Users/Jieming/Documents/thesis/lynne_work/motifVar/ank/blosum_sift_polyphen2_gerp")
source("C:/Users/Jieming/Documents/Shared/scripts-R_perl_shell_macros/R codes/jmRlib.R")
library(ggplot2)
library(matrixcalc)
library(vioplot)

boxbubbleplot <- function(title,resNumCol,scoreCol,lastCol,data, cexAxisFontSize=2)
{
  title = title
  resNum = data[,resNumCol]
  blosum = data[,scoreCol]
  m=max(blosum)
  tt.blosum=table(resNum, blosum)
  range01 = function(x){(x-min(x))/(max(x)-min(x))}
  # #   x11()
  #   par(mar=c(5,6,4,4),xpd=TRUE,mar=c(1, 4, 1, 1) + 1)
  #   boxplot(blosum~resNum,xlab="position",ylab=paste(title,sep=" "),cex.lab=2, cex.axis=cexAxisFontSize, las=2)
  #   
  #   ## convert table counts of SNVs at each position to dataframe
  #   a.blosum=as.data.frame(table(resNum))[,2]
  #   text(0,m+0.5,paste(a.blosum,collapse=" "),adj = c(0,0))
  #   
  #   ###### bubble plot
  #   df = expand.grid(sort(unique(resNum)),sort(unique(blosum)))
  #   aa = as.vector(table(resNum, blosum))
  #   df$value = aa
  #   
  #   # x11()
  #   area = pi * (df$value/150)^2
  #   par(mar=c(5,6,4,4),xpd=TRUE)
  #   area = (range01(df$value)/2)^2 * pi
  #   symbols(df$Var1,df$Var2,circles=area, inches=F, bg=rgb(1,0,0,0.5), fg="black", yaxt="n", 
  #           xaxt="n",add=T)
  #   
  #   # axis(1,at=seq(1,34), las=2)
  #   # axis(2,at=seq(-4,9))
  #   # text(0,3,paste(a.blosum,collapse=" "),adj = c(0,0))
  
  ## ordered
#     x11()
  par(xpd=TRUE,mar=c(3, 4, 1, 1) + 1)
  boxplot(blosum~data[,lastcol],xlab="position",ylab=paste(title,sep=" "),cex.lab=2, cex.axis=cexAxisFontSize, las=2)
  a.blosumordered=as.data.frame(table(data[,lastcol]))[,2]
  text(0,m+0.5,paste(a.blosumordered,collapse=" "),adj = c(0,0))
  
  #### bubble
#   x11()
  par(xpd=TRUE)
  df.ordered = expand.grid(sort(unique(data[,lastcol])),sort(unique(blosum)))
  aa.ordered = as.vector(table(data[,lastcol], blosum))
  df.ordered$value = aa.ordered
  area.ordered = (range01(df.ordered$value)/2)^2 * pi
  symbols(df.ordered$Var1,df.ordered$Var2,circles=area.ordered, inches=F, bg=rgb(1,0,0,0.5), fg="black", yaxt="n", 
          xaxt="n",xlab='position',ylab=paste(title, "scores",sep=" "),add=T)
  # axis(1,at=seq(1,34),labels=as.vector(sort(unique(data[,lastcol]))), las=2)
  # axis(2,at=seq(-4,9))
  
  #   return(blosum);
}

####################################################################

filename = 'final.merged.ExAC.ens73.ANK.30aa.vepoutput.coding.canonical.noORF.sorted.auto.blosum-30-35-40-45-50-62.gerp.sift.bed'
data1 = read.table(filename, header=F, stringsAsFactors = F, comment.char="", sep="\t")

## remove AC=1 V8
data = data1[data1$V8 > 1,]

order = read.table('order_relentropy-ank.txt', header=T, stringsAsFactors = F, comment.char="", sep="\t")
lastcol=ncol(data)+1
resNumCol=45

## we insert ordered info into lastCol so order info is lastCol
data[,lastcol] = factor(data[,resNumCol], levels=levels(data[,resNumCol])<-order$pos)

## remove synonymous
data.ns.only <- data[(data$V56 == "NON_SYNONYMOUS_CODING" | data$V56 == "NON_SYNONYMOUS_CODING,SPLICE_SITE"),]




## condel ####
## conseq "neutral"     NA            "deleterious"
## resNum, condel score 
title = 'Condel'
scoreCol = 69
tt.condel=table(data[,resNumCol], data[,scoreCol])

## convert table counts of SNVs at each position to dataframe; nonmissing
condel.nonmissData=data[complete.cases(data[,scoreCol]),]

x11()
boxbubbleplot(title,resNumCol,scoreCol,lastCol,condel.nonmissData, cexAxisFontSize=2)

## polyphen ####
## conseq "benign"            NA                  "probably_damaging" "possibly_damaging" "unknown"
## resNum, polyphen score 
title = 'Polyphen'
scoreCol = 72
tt.polyphen=table(data[,resNumCol], data[,scoreCol])

## convert table counts of SNVs at each position to dataframe; nonmissing
polyphen.nonmissData=data[complete.cases(data[,scoreCol]),]

x11()
boxbubbleplot(title,resNumCol,scoreCol,lastCol,polyphen.nonmissData, cexAxisFontSize=2)

## sift ######
## conseq "tolerated"   NA            "deleterious"
## resNum, sift score 
title = 'SIFT'
maxResNum = max(data[,resNumCol])
siftscorecol = 80
tt.sift=table(data[,resNumCol], 1-data[,siftscorecol])

## convert table counts of SNVs at each position to dataframe; nonmissing
sift.nonmissData=data[complete.cases(data[,siftscorecol]),]
sift.nonmissData[,siftscorecol] = 1-sift.nonmissData[,siftscorecol]

## plot ORDERED
## violin plot
sift.ordered = array(list(NULL), c(maxResNum,1))
ctr = 1
for(i in order[,1])
{
  sift.ordered[[ctr]] = sift.nonmissData[,siftscorecol][sift.nonmissData[,lastcol] == i] 
  ctr = ctr + 1
}

x11()
# par(mfrow=c(2,2),tcl=-0.5)
par(cex.axis=1.5,cex.lab=1.5)
do.call(vioplot.las2,c(sift.ordered,col="green",list(names=order[,1])))
mtext("position",side=1, line=2.5, cex=1.5)
mtext("1-SIFT",side=2, line=2.5, cex=1.5)

## count number of SNPs in each pos
counts.pos=as.data.frame(table(sift.nonmissData[,lastcol]))[,2]
m=max(sift.nonmissData[,siftscorecol])
text(1,m,paste(counts.pos,collapse="   "),adj = c(0,0))
                                                         

## boxbubbleplot
# boxbubbleplot(title,resNumCol,siftscorecol,lastCol,sift.nonmissData, cexAxisFontSize=2)

## gerp ####
## resNum, gerp score 
title = 'GERP'
gerpscoreCol = 71
tt.gerp=table(data[,resNumCol], data[,gerpscoreCol])

## in sequence
gerp.nonmissData=data[complete.cases(data[,gerpscoreCol]),]
gerp.NSData=data[complete.cases(data[,siftscorecol]),]

## ordered
gerp.nonmissData.ordered = array(list(NULL), c(maxResNum,1))
gerp.NSData.ordered = array(list(NULL), c(maxResNum,1))
ctr = 1
for(i in order[,1])
{
  gerp.nonmissData.ordered[[ctr]] = data[,gerpscoreCol][data[,lastcol] == i] 
  gerp.NSData.ordered[[ctr]] = data.ns.only[,gerpscoreCol][data.ns.only[,lastcol] == i] 
  ctr = ctr + 1
}


## plot gerp all S + NS ordered
x11()
par(cex.axis=1.5,cex.lab=1.5,mai=c(0.7,0.7,0.05,0.3))
do.call(vioplot.las2,c(gerp.nonmissData.ordered,col="red",list(names=order[,1])))
mtext("position",side=1, line=2.5, cex=1.5)
mtext("GERP all",side=2, line=2, cex=1.5)

## count number of SNPs in each pos
counts.pos=sapply(gerp.nonmissData.ordered,length)
m=max(data[,gerpscoreCol])
text(0.7,m,paste(counts.pos,collapse=" "),adj = c(0,0))

## bubbleboxplot
# boxbubbleplot(title,resNumCol,gerpscoreCol,lastCol,gerp.nonmissData, cexAxisFontSize=2)

## gerp NS only ordered
x11()
par(cex.axis=1.5,cex.lab=1.5,mai=c(0.7,0.7,0.05,0.3))
do.call(vioplot.las2,c(gerp.NSData.ordered,col="red",list(names=order[,1])))
mtext("position",side=1, line=2.5, cex=1.5)
mtext("GERP NS-only",side=2, line=2, cex=1.5)

## count number of SNPs in each pos
counts.pos=sapply(gerp.NSData.ordered,length)
m=max(data[,gerpscoreCol])
text(0.7,m,paste(counts.pos,collapse=" "),adj = c(0,0))

## bubbleboxplot
# boxbubbleplot(title,resNumCol,gerpscoreCol,lastCol,gerp.NSData, cexAxisFontSize=2)


## correlation scatter plots ##
### gerp vs blosum30 ####
x11()
par(xpd=TRUE,cex.axis=1, cex.lab=1.5, cex.main=1.5)
range.x = range(gerp.NSData$V76)
range.y = range(gerp.NSData$V77)
plot(gerp.NSData$V76,gerp.NSData$V77,main=paste("NSonly"),xlab='GERP-NSonly',ylab='BLOSUM30-NSonly', xlim=range.x, ylim=range.y)

rareVar=gerp.NSData[(gerp.NSData$V20 <= 0.005),]
commVar=gerp.NSData[(gerp.NSData$V20 > 0.005),]
# par(new=T)
# plot(rareVar$V76,rareVar$V77,col='red',xaxt='n',yaxt='n',xlab='',ylab='',pch=3, xlim=range.x, ylim=range.y)
par(new=T)
plot(commVar$V76,commVar$V77,col='red',xaxt='n',yaxt='n',xlab='',ylab='',pch=3, xlim=range.x, ylim=range.y, lwd=2)
legend(-12,-4,"comm",col="red",pch=3, cex=1.5, lwd=2)

## by position
resNum = 20
NSData = gerp.NSData[gerp.NSData[,resNumCol]==resNum,]
x11()
range.x = range(NSData$V76)
range.y = range(NSData$V77)
plot(NSData$V76,NSData$V77,main=paste("NSonly at pos",resNum),xlab='GERP',ylab='BLOSUM30',xlim=range.x, ylim=range.y,)
rareVar.pos=NSData[(NSData$V20 <= 0.005),]
commVar.pos=NSData[(NSData$V20 > 0.005),]
# par(new=T)
# plot(rareVar.pos$V76,rareVar.pos$V77,col='red',xaxt='n',yaxt='n',xlab='',ylab='',pch=3, xlim=range.x, ylim=range.y)
par(new=T)
plot(commVar.pos$V76,commVar.pos$V77,col='red',xaxt='n',yaxt='n',xlab='',ylab='',pch=3, xlim=range.x, ylim=range.y, lwd=2)
legend(1,-4,"comm",col="red",pch=3, cex=1.5, lwd=2)

# gerp vs condel
x11()
range.x = range(gerp.NSData$V76)
range.y = range(!is.na(gerp.NSData$V69))
plot(gerp.NSData$V76,gerp.NSData$V69,main=paste("NSonly all pos"),xlab='GERP',ylab='Condel', xlim=range.x, ylim=range.y)
par(new=T)
plot(rareVar$V76,rareVar$V69,col='red',xaxt='n',yaxt='n',xlab='',ylab='',pch=3, xlim=range.x, ylim=range.y)
par(new=T)
plot(commVar$V76,commVar$V69,col='blue',xaxt='n',yaxt='n',xlab='',ylab='',pch=3, xlim=range.x, ylim=range.y)

## by position
resNum = 11
NSData = gerp.NSData[gerp.NSData[,resNumCol]==resNum,]
x11()
plot(NSData$V76,NSData$V69,main=paste("NSonly at pos",resNum),xlab='GERP',ylab='Condel', xlim=range.x, ylim=range.y)


### condel VS blosum30 ####
x11()
plot(gerp.NSData$V69,gerp.NSData$V77,main=paste("NSonly all pos"),xlab='Condel',ylab='BLOSUM30')

par(new=T)
plot(rareVar$V69,rareVar$V77,col='red',xaxt='n',yaxt='n',xlab='',ylab='')
par(new=T)
plot(commVar$V69,commVar$V77,col='blue',xaxt='n',yaxt='n',xlab='',ylab='')


## by position
resNum = 20
NSData = gerp.NSData[gerp.NSData[,resNumCol]==resNum,]
x11()
plot(NSData$V76,NSData$V69,main=paste("NSonly at pos",resNum),xlab='Condel',ylab='BLOSUM30')

### 1-sift VS blosum30 ####
blosum30col = 65
mafcol = 20
x11()
range.x = range(sift.nonmissData[,siftscorecol])
range.y = range(sift.nonmissData[,blosum30col])
plot(sift.nonmissData$V75,sift.nonmissData[,siftscorecol],main=paste("NSonly"),xlab='1-SIFT-NSonly',ylab='BLOSUM30-NSonly', xlim=range.x, ylim=range.y)
rareVar=sift.nonmissData[(gerp.NSData[,mafcol] <= 0.005),]
commVar=sift.nonmissData[(gerp.NSData[,mafcol] > 0.005),]
# par(new=T)
# plot(rareVar$V75,rareVar$V77,col='red',xaxt='n',yaxt='n',xlab='',ylab='')
par(new=T)
plot(commVar[,siftscorecol],commVar[,blosum30col],col='blue',xaxt='n',yaxt='n',xlab='',ylab='',pch=3, xlim=range.x, ylim=range.y)
legend(0.11,-4,"comm",col="blue",pch=3, cex=1.5, lwd=2)

### GERP NS vs 1-sift ####
x11()
range.x = range(sift.nonmissData[,gerpscorecol])
range.y = range(sift.nonmissData[,siftscorecol])
plot(sift.nonmissData[,gerpscorecol],sift.nonmissData[,siftscorecol],main=paste("NSonly"),xlab='GERP-NSonly',ylab='1-SIFT-NSonly', xlim=range.x, ylim=range.y)
rareVar=sift.nonmissData[(gerp.NSData[,mafcol] <= 0.005),]
commVar=sift.nonmissData[(gerp.NSData[,mafcol] > 0.005),]
par(new=T)
plot(commVar[,gerpscorecol],commVar[,siftscorecol],col='blue',xaxt='n',yaxt='n',xlab='',ylab='',pch=3, xlim=range.x, ylim=range.y)
legend(-12.5,0.2,"comm",col="blue",pch=3, cex=1.5, lwd=2)

