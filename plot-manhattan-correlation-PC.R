
#------------------------------------------------------#
#                                                      #
#                                                      #                                    
#  Rscript to plot correlation of PC with SNPs         #
# Date: July 15, 2008                                  #
# Version 1                                            #   
#                                                      #   
#------------------------------------------------------#

## The correlations are absolute valued and are computed between
## the normalized genotype data (SNPs) and the principal components

# library(GDD)

#-------------------------
# Read marker file
#-------------------------
#setwd("C:/Documents and Settings/chenjm/Desktop/work_documentation/JJ/surnames/2-pca/5-pca on 3971sa-373278snps no PC1-cor 0.17 snps")
#setwd("C:/Documents and Settings/chenjm/Desktop/work_documentation/Mark/indian-reich-hapmap-pasnp-ihp/pca")
#setwd("C:/Documents and Settings/chenjm/Desktop/work_documentation/Mark/IBD/colitis-devindri/pca")
#setwd("C:/Documents and Settings/chenjm/Desktop/work_documentation/Affymetrix/axiom-asi/1000samples-evaluation/qc/pca+manhattan")
# setwd("C:/Documents and Settings/chenjm/Desktop/work_documentation/Mark/amerindians inmegen/6 merge mtz")
setwd("C:/Users/JM/Desktop")

marker_temp=read.table("1437snps.mk",header=T) ## includes snp-id, chromosome, position
dim(marker_temp)
head(marker_temp)
marker=marker_temp[,c(1,2,3)] ## change this to get the relevant header
dim(marker)
head(marker)
rm(marker_temp)
gc()


#-------------------------
# Read correlation file
#-------------------------
filecor = "1437snps-posterior.cor"
corr=read.table(filecor,header=T)
dim(corr)
head(corr)

#-------------------------
# Input header for plots
#-------------------------
header=filecor


#---------------------------------------------------------
# Combine two files and order by chromosome and position
#---------------------------------------------------------
#data_temp=merge(marker,corr,by.x="snp.id",by.y="SNP",sort=F)
data_temp=merge(marker,corr,by.x="snp.id",sort=F)
data_temp[,"chromosome"]=factor(data_temp$chromosome,levels=c(1:22,"M","X","XY","Y"),labels=1:26)
dim(data_temp)
head(data_temp)
summary(as.factor(data_temp$chromosome))
#rm(marker,corr)
gc()

data=data_temp[order(data_temp$chromosome,data_temp$position),]
rownames(data)=NULL
dim(data)
head(data)
rm(data_temp)
gc()


#---------------------------------------------------------
# Enter known highly-conserved regions
#---------------------------------------------------------

# conserved=matrix(NA,nrow=5,ncol=3)
# colnames(conserved)=c("chromosome","start","end")
# conserved[1,]=c(8,0,12700000)
# conserved[2,]=c(6,25e+06,34e+06)
# conserved[3,]=c(11,45e+06,57e+06)
# conserved[4,]=c(5,44e+06,51.5e+06)
# conserved[5,]=c(11,84e+06,86e+06)
# conserved


#------------------------------------------------------------------------------
# Map out the row no. of closest SNP for each of the highly conserved region
#------------------------------------------------------------------------------

# For each highly conserved region,
# border_start gives the row number in "data" of the SNP closest to the start of the highly conserved region
# border_end gives the row number in "data" of the SNP closest to the end of the highly conserved region

# border=matrix(NA,ncol=2,nrow=nrow(conserved))
# for(i in 1:nrow(conserved)){
# 
# 	chr=conserved[i,"chromosome"]
# 	start=conserved[i,"start"]
# 	end=conserved[i,"end"]
# 	data2=data[data$chromosome==chr,]
# 
# 	diff_start=abs(data2$position-start)
# 	diff_end=abs(data2$position-end)
# 	data3=cbind(data2,diff_start,diff_end)
# 
# 	border[i,1]=as.numeric(rownames(data3[data3$diff_start==min(abs(data3$position-start)),]))
# 	border[i,2]=as.numeric(rownames(data3[data3$diff_end==min(abs(data3$position-end)),]))
# 	rm(chr,start,end,data2,data3,diff_start,diff_end)
# 
# }
# 
# conserved_regions=cbind(conserved,border)
# colnames(conserved_regions)=c("chromosome","start","end","border_start","border_end")
# conserved_regions
# #rm(conserved,border)
# gc()


#-----------------------------------
# Plot entire genome for each PC
#-----------------------------------

color=c("red2","green","orange1","royalblue","yellow3","darkslategrey","purple3","turquoise","hotpink","lightgreen","salmon","skyblue1","goldenrod","slategrey","purple1","maroon","darkgreen","orange4","darkblue","brown","gray","mediumpurple")

## note that PC1 starts from col4, snpid col1 chr col2 pos col3
for(j in 4:7){

	windows()

	GDD(file=paste(header,"-cor-",colnames(data)[j],".png",sep=""),height=1000,width=1000,ps=18) 

	plot(1:nrow(data),abs(data[,j]),xlim=c(0,nrow(data)),ylim=c(0,max(abs(data[,j]))),pch=19,col=color[data[,"chromosome"]],xaxt="n",main=header,xlab="Chromosome",ylab=paste("Correlation with ",colnames(data)[j],sep=""))

	for(k in 1:nrow(conserved_regions)){

		x1=c(conserved_regions[k,"border_start"]-1,conserved_regions[k,"border_start"]-1)
		x2=c(conserved_regions[k,"border_end"]+1,conserved_regions[k,"border_end"]+1)
		y=c(0,max(abs(data[,j])))
		lines(x1,y,col="black")
		lines(x2,y,col="black")
		abline(h=0.88, col = "red")
		#text(x=-1, y=0.19, labels="0.17", col="black")
		rm(x1,x2,y)

	}

	axis(1,at=median(as.numeric(rownames(data[data$chromosome==1,]))),1)
	axis(1,at=median(as.numeric(rownames(data[data$chromosome==2,]))),2)
	axis(1,at=median(as.numeric(rownames(data[data$chromosome==3,]))),3)
	axis(1,at=median(as.numeric(rownames(data[data$chromosome==4,]))),4)
	axis(1,at=median(as.numeric(rownames(data[data$chromosome==5,]))),5)
	axis(1,at=median(as.numeric(rownames(data[data$chromosome==6,]))),6)
	axis(1,at=median(as.numeric(rownames(data[data$chromosome==7,]))),7)
	axis(1,at=median(as.numeric(rownames(data[data$chromosome==8,]))),8)
	axis(1,at=median(as.numeric(rownames(data[data$chromosome==9,]))),9)
	axis(1,at=median(as.numeric(rownames(data[data$chromosome==10,]))),10)
	axis(1,at=median(as.numeric(rownames(data[data$chromosome==11,]))),11)
	axis(1,at=median(as.numeric(rownames(data[data$chromosome==12,]))),12)
	axis(1,at=median(as.numeric(rownames(data[data$chromosome==13,]))),13)
	axis(1,at=median(as.numeric(rownames(data[data$chromosome==14,]))),14)
	axis(1,at=median(as.numeric(rownames(data[data$chromosome==15,]))),15)
	axis(1,at=median(as.numeric(rownames(data[data$chromosome==16,]))),16)
	axis(1,at=median(as.numeric(rownames(data[data$chromosome==17,]))),17)
	axis(1,at=median(as.numeric(rownames(data[data$chromosome==18,]))),18)
	axis(1,at=median(as.numeric(rownames(data[data$chromosome==19,]))),19)
	axis(1,at=median(as.numeric(rownames(data[data$chromosome==20,]))),20)
	axis(1,at=median(as.numeric(rownames(data[data$chromosome==21,]))),21)
	axis(1,at=median(as.numeric(rownames(data[data$chromosome==22,]))),22)
	
	dev.off()
}
