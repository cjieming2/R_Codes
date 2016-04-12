#-------------------------
# Script to generate manhattan plot of correlation plot
#-------------------------

#setwd("D:/DataAnalysisWorkingPath/")
#setwd("D:/DataAnalysisWorkingPath/Sonia - Meningococcal Studies/58c-nbs-1.2M-ukmgc1-combined-trend/")
setwd("C:/Documents and Settings/chenjm/Desktop/work_documentation/Mark/amerindians inmegen/4 merge maya tepe")

corr_col = 2
currentPC = "PC1"

#minuslogPlimit = 15
#suggestiveLogPValue = 5
#significantLogPValue = 7

plottitle=paste(resultfile," Manhattan Plot")
color=c("red2","green","orange1","royalblue","yellow3","darkslategrey","purple3","turquoise","hotpink","lightgreen","salmon","skyblue1","goldenrod","slategrey","purple1","maroon","darkgreen","orange4","darkblue","brown","gray","mediumpurple")

#-------------------------
# Read marker file
#-------------------------
mkfile = "cut-22806snps-mayas-tepehuans-affy50kxbaI-500k-illumina550k-inmegen-fwd.mk" 

## change this to get the relevant header
mk_col_snpid = 1
mk_col_chromosome = 2
mk_col_position = 3

marker_temp=read.table(mkfile,header=T)
dim(marker_temp)
head(marker_temp)

marker=marker_temp[,c(mk_col_snpid,mk_col_chromosome,mk_col_position)]		#extract column for snp, chromosome, position

dim(marker)
head(marker)
rm(marker_temp)
gc()

#-------------------------
# Read correlation file
#-------------------------
corrfile = "62samples-22806snps-mayas-tepehuans-affy50kxbaI-500k-illumina550k-inmegen-fwd.cor"
#result_col_snpid = 2
#result_col_pval = 10

corr_col_snpid = 1

corr_temp=read.table(corrfile,header=T)
dim(corr_temp)
head(corr_temp)

corr=corr_temp[,c(corr_col_snpid,corr_col)]

dim(corr)
head(corr)
rm(corr_temp)
gc()

#-------------------------
# Combine two files
#-------------------------

data_temp = merge(marker, corr, by.x = "snp.id", sort=F)
#data_temp$logP = -1 * log10(data_temp[,4])

dim(data_temp)
head(data_temp)

rm(marker, corr)
gc()

data_temp[,"chromosome"]=factor(data_temp$chromosome,levels=c(1:22,"M","X","XY","Y"),labels=1:26)
summary(as.factor(data_temp$chromosome))

data=data_temp[order(data_temp$chromosome,data_temp$position),]
rownames(data)=NULL
dim(data)
head(data)

rm(data_temp)
gc()

#-------------------------
# Calculate relative chromosome position and offset
#-------------------------
chrPosition <- data.frame(chr=numeric(0), min=numeric(0), max=numeric(0), length=numeric(0))

for(i in 1:length(levels(data$chromosome))){
  if(length(which(data$chromosome == i)) > 0) {
    tmp.min =  min(data$position[which(data$chromosome == i)])
    tmp.max =  max(data$position[which(data$chromosome == i)])
    chrPosition <- rbind(chrPosition, data.frame(chr=i, min=tmp.min, max=tmp.max, length=tmp.max-tmp.min))
  }
}

tmp.sum = sum(as.numeric(chrPosition$length))
chrPosition$proportion = (chrPosition$length / tmp.sum) * 100

chrPosition$offset = 0
if (length(chrPosition$chr) > 1) {
  chrPosition$offset[2] =  chrPosition$proportion[1]
}

if (length(chrPosition$chr) > 2) {
  for(i in 3:length(chrPosition$chr)) {
    chrPosition$offset[i] = chrPosition$offset[i-1] + chrPosition$proportion[i-1]
  }
}

#Calculate the relative X axis coordinate across all chromosome
data$plotX = (((data$position - chrPosition$min[data$chromosome]) / chrPosition$length[data$chromosome]) * chrPosition$proportion[data$chromosome]) + chrPosition$offset[data$chromosome]


#-------------------------
# Extract subset of data under and above limit for separate plotting
#-------------------------
#overLimit = which(data$logP > minuslogPlimit)
#belowLimit = which(data$logP <= minuslogPlimit)
#numOverLimit = length(overLimit)


#-------------------------
# Generate plot
#-------------------------
windows()

#windows(width=20,height=10)
#GDD(file=paste(header,"-cor-",colnames(data)[j],".png",sep=""),height=1000,width=1000,ps=18)
#plot(1:nrow(data3),data3[,5],xlim=c(0,nrow(data3)),ylim=c(0,minuslogPlimit),pch=19,col=color[data3[,"chromosome"]],xaxt="n",main="Title",xlab="Chromosome",ylab="-log10 P Value",ann="y")
#plot(belowLimit,data3[belowLimit,5],xlim=c(0,nrow(data3)),ylim=c(0,minuslogPlimit),pch=19,col=color[data3[belowLimit,"chromosome"]],xaxt="n",main="Title",xlab="Chromosome",ylab="-log10 P Value",ann="y")
#plot(belowLimit,data[belowLimit,5],xlim=c(0,nrow(data)),ylim=c(0,minuslogPlimit),pch=19,col=color[data[belowLimit,"chromosome"]],xaxt="n",ann="n")
#title(main=plottitle, xlab="Chromosome", ylab="-log10(P)")

#plot(data[belowLimit,"plotX"],data[belowLimit,"logP"],xlim=c(0,100),ylim=c(0,minuslogPlimit),pch=19,col=color[data[belowLimit,"chromosome"]],xaxt="n",ann="n")

column_desired = corr_col + 2
plot(1:nrow(data),abs(data[,column_desired]),xlim=c(0,nrow(data)),ylim=c(0,max(abs(data[,column_desired]))),pch=19,col=color[data[,"chromosome"]],xaxt="n",xlab="Chromosome",ylab=paste("Correlation with ",currentPC,sep=""))
title(main=plottitle)

#plot(1:nrow(overLimit),minuslogPlimit,xlim=c(0,nrow(data3)),ylim=c(0,minuslogPlimit),pch=24,col=color[data3[,"chromosome"]],xaxt="n",main=plottitle,xlab="Chromosome",ylab="-log10 P Value")
#plot(1:nrow(overLimit),minuslogPlimit,xlim=c(0,nrow(data3)),ylim=c(0,minuslogPlimit),pch=24,col=color[data3[,"chromosome"]],xaxt="n", ann="n")
#plot(1:nrow(overLimit),minuslogPlimit,pch=24,col=color[data3[,"chromosome"]],xaxt="n", ann="n")
#points(overLimit, rep(minuslogPlimit, time=numOverLimit), pch=24, col=color[data3[overLimit,"chromosome"]], bg=color[data3[overLimit,"chromosome"]])

#points(overLimit, rep(minuslogPlimit, time=numOverLimit), pch=24, col=color[data[overLimit,"chromosome"]], bg=color[data[overLimit,"chromosome"]])

points(data[overLimit,"plotX"], rep(minuslogPlimit, time=numOverLimit), pch=19, col=color[data[overLimit,"chromosome"]], bg=color[data[overLimit,"chromosome"]])


#abline(h=significantLogPValue , lwd=2,col="red")
#abline(h=suggestiveLogPValue , lwd=2,col="dodgerblue2")


for(i in 1:length(chrPosition$chr)) {
	axis(1,at=chrPosition$offset[i] + (chrPosition$proportion[i] / 2), chrPosition$chr[i])
}

#axis(1,at=median(as.numeric(rownames(data[data$chromosome==1,]))),1)
#axis(1,at=median(as.numeric(rownames(data[data$chromosome==2,]))),2)
#axis(1,at=median(as.numeric(rownames(data[data$chromosome==3,]))),3)
#axis(1,at=median(as.numeric(rownames(data[data$chromosome==4,]))),4)
#axis(1,at=median(as.numeric(rownames(data[data$chromosome==5,]))),5)
#axis(1,at=median(as.numeric(rownames(data[data$chromosome==6,]))),6)
#axis(1,at=median(as.numeric(rownames(data[data$chromosome==7,]))),7)
#axis(1,at=median(as.numeric(rownames(data[data$chromosome==8,]))),8)
#axis(1,at=median(as.numeric(rownames(data[data$chromosome==9,]))),9)
#axis(1,at=median(as.numeric(rownames(data[data$chromosome==10,]))),10)
#axis(1,at=median(as.numeric(rownames(data[data$chromosome==11,]))),11)
#axis(1,at=median(as.numeric(rownames(data[data$chromosome==12,]))),12)
#axis(1,at=median(as.numeric(rownames(data[data$chromosome==13,]))),13)
#axis(1,at=median(as.numeric(rownames(data[data$chromosome==14,]))),14)
#axis(1,at=median(as.numeric(rownames(data[data$chromosome==15,]))),15)
#axis(1,at=median(as.numeric(rownames(data[data$chr==16,]))),16)
#axis(1,at=median(as.numeric(rownames(data[data$chromosome==17,]))),17)
#axis(1,at=median(as.numeric(rownames(data[data$chromosome==18,]))),18)
#axis(1,at=median(as.numeric(rownames(data[data$chromosome==19,]))),19)
#axis(1,at=median(as.numeric(rownames(data[data$chromosome==20,]))),20)
#axis(1,at=median(as.numeric(rownames(data[data$chromosome==21,]))),21)
#axis(1,at=median(as.numeric(rownames(data[data$chromosome==22,]))),22)


#if (minuslogPlimit < 20) {
#	for(i in 1:minuslogPlimit) {
#		axis(2,at=i,i)
#	}
#}
