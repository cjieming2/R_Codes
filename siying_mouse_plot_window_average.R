setwd('C:/Users/Jieming/Documents/thesis/siying/')


#### ko_433.txt ######
file = 'ko_433.txt'
data = read.table(file,header=F,stringsAsFactors=F)
data = as.data.frame(data[4:(nrow(data)-10),])

## plot
xlt=c(0,nrow(data))
ylt=c(0,0.23)
x11()
par(mar=c(5,5,1,1),cex.axis=1, cex.lab=2, cex.main=2)
# plot(seq(1,nrow(data)),data[,1], xlab="time",ylab="amplitude",pch=19, type="b", xlim=xlt, ylim=ylt, col="black")


## define number of points in one window size
wsize=35
## define the points in each window to subset
wstart=seq(1,nrow(data),wsize)

## calc the average in each window
# aves = as.data.frame(sapply(wstart, function(x) mean(data[x:(x+wsize-1),])))
aves = as.data.frame(sapply(wstart, function(x) max(data[x:(x+wsize-1),]) - min(data[x:(x+wsize-1),])))

## normalize averages
aves_norm = aves / max(aves)

## plot averages
# x11()
# par(mar=c(5,5,1,1),cex.axis=1, cex.lab=2, cex.main=2)
# par(new=T)
# plot(seq(1,nrow(aves)),aves[,1], xaxt="n",yaxt="n", xlab="", ylab="",pch=19, type="b", col="red", ylim=ylt) 
# plot(seq(1,nrow(aves)),aves[,1], xlab="time",ylab="amplitude",pch=19, type="b", col="red", ylim=ylt) 
plot(seq(1,nrow(aves_norm)),aves_norm[,1], xlab="time",ylab="amplitude",pch=19, type="b", col="red", ylim=c(0,1)) 



## calc amplitude using the average
# datanew = as.data.frame(unlist(lapply(seq(1,length(wstart)),
#                                function(x) abs((data[wstart[x]:(wstart[x]+wsize-1),] - aves[x,])))))
# plot(seq(1,nrow(datanew)),datanew[,1], xlab="time",ylab="amplitude",pch=19, type="b",xlim=xlt)

write.table(as.matrix(cbind(aves,aves_norm)),file=paste(file,"2xampnorm.txt",sep=""), 
            row.names=FALSE, col.names=FALSE, append=FALSE)


#### fBLACK WT ######
file_b = 'fBLACK.txt'
data_b = read.table(file_b,header=F,stringsAsFactors=F)
data_b = as.data.frame(data_b[(55021):(nrow(data_b)),])

## plot
# xlt_b=c(0,nrow(data))
# ylt_b=c(0,0.23)
x11()
par(mar=c(5,5,1,1),cex.axis=1, cex.lab=2, cex.main=2)
plot(seq(1,nrow(data_b)),data_b[,1], xlab="time",ylab="amplitude",pch=19, type="b", col="black")

## define number of points in one window size
wsize_b=3003
## define the points in each window to subset
wstart_b=seq(1,nrow(data_b),wsize_b)

## calc the 2xamplitude in each window
aves_b = as.data.frame(sapply(wstart_b, function(x) max(data_b[x:(x+wsize_b-1),]) - min(data_b[x:(x+wsize_b-1),])))

## normalize 2xamplitude by the max value
aves_norm_b = aves_b / max(aves_b)

## plot 2xamplitude
x11()
plot(seq(1,nrow(aves_norm_b)),aves_norm_b[,1], xlab="time",ylab="amplitude",pch=19, 
     type="b", col="red", ylim=c(0,1)) 

write.table(as.matrix(cbind(aves_b,aves_norm_b)),file=paste(file_b,"2xampnorm.txt",sep=""), 
            row.names=FALSE, col.names=FALSE, append=FALSE)