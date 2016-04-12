setwd('C:/Users/Jieming/Documents/thesis/siying/JimHokanson-adinstruments_sdk_matlab-17a5368')


# file = 'siying433.txt'
file = 'fBLACK.txt'
# file = 'fWHITE.txt'
data = read.table(file,header=T,stringsAsFactors=F)
data_new = as.data.frame(data[40000:nrow(data),])
newseq100 = seq(1,nrow(data_new),100)
newseq1000 = seq(1,nrow(data_new),1000)

data_new100 = as.data.frame(data_new[newseq100,])
data_new1000 = as.data.frame(data_new[newseq1000,])

## plot
x11()
par(mar=c(5,5,1,1),cex.axis=1, cex.lab=2, cex.main=2)
plot(seq(1,nrow(data)),data[,1], xlab="time",ylab="amplitude",pch=15)
x11()
par(mar=c(5,5,1,1),cex.axis=1, cex.lab=2, cex.main=2)
plot(seq(1,nrow(data_new)),data_new[,1], xlab="time",ylab="amplitude",pch=19)
x11()
par(mar=c(5,5,1,1),cex.axis=1, cex.lab=2, cex.main=2)
plot(seq(1,nrow(data_new100)),data_new100[,1], xlab="time",ylab="amplitude",pch=19, type="b")
x11()
par(mar=c(5,5,1,1),cex.axis=1, cex.lab=2, cex.main=2)
plot(seq(1,nrow(data_new1000)),data_new1000[,1], xlab="time",ylab="amplitude",pch=19, type="b")


write(as.matrix(data_new100),file=paste(file,"_new100.txt",sep=""), ncolumns=1, append=FALSE)
write(as.matrix(data_new1000),file=paste(file,"_new1000.txt",sep=""), ncolumns=1, append=FALSE)