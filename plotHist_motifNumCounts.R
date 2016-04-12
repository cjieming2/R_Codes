setwd("C:/Users/JM/thesis/lynne_work/motifVar/motifNumCounts")


## read files
ifile       = "ALL_ANK.prot.new"
data        = read.delim(ifile,header=T,sep="\t")

## histograms
numofbins = max(data$motifNum) + 4
oriMotifNum = hist(data$motifNum,breaks=numofbins)
newMotifNum = hist(data$new_motifNum,breaks=numofbins)

## plots bars
# x11()
plot(newMotifNum,col=rgb(0,0,1,1/4),main=ifile,xlab="motifNum") #blue
plot(oriMotifNum,col=rgb(1,0,0,1/4),add=T) #red
# axis(side=1, at=seq(1,100, 1), labels=seq(1,100,1))

# ## plot line plots
# plot(newMotifNum$density,col="red",main=ifile,type="l",xlim=c(0,20),ylim=c(0,0.3)) #blue
# par(new=TRUE)
# plot(oriMotifNum$density,col="blue",type="l",xaxt='n',yaxt='n',ann=FALSE,xlim=c(0,20),ylim=c(0,0.3)) #red

# ## we have same number so we can use probability kernel fxns to compare
# oriMotifNum_d = density(data$motifNum)
# newMotifNum_d = density(data$new_motifNum)
# 
# ## plots kernel densities
# library(sm)
# 
# # combine the 2 vectors
# motifNum_c = c(oriMotifNum$density,newMotifNum$density)
# labels_c   = c(rep(1,times=length(oriMotifNum$density)),rep(2,times=length(newMotifNum$density)))
# # create value labels
# lc.f = factor(labels_c, levels= c(1,2),
#                 labels = c("oriMotifNum", "newMotifNum")) 
# 
# # plot densities 
# sm.density.compare(motifNum_c, labels_c, xlab="numMotifs")
# title(main=ifile)
# 
# # add legend via mouse click
# colfill<-c(2:(2+length(levels(lc.f)))) 
# legend(locator(1), levels(lc.f), fill=colfill)