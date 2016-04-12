#####
#THIS PROGRAM CALCULATES VALUES FOR DETERMINING NUMBER OF PRINCIPLE
#COMPONENTS TO RETAIN.  IT COMPUTES THE FIRST 10 CUT OFF VALUES FOR
#TWO PARALLAL METHODS (LONGMAN AND ALLEN) AND THE BROKEN STICK METHOD
#(FRONTIER).  ;
######

#n is the number of observations in the dataset
#p is the number of variables in the dataset

parallel<-function(n,p){

if (n>1000 || p>100) {
    print ("Sorry, this only works for n<1000 and p<100")
    stop()
}

coefs<-matrix(
c(0.0316, 0.7611, -0.0979, -0.3138, 0.9794, -.2059, .1226, 0, 0.1162, 
0.8613, -0.1122, -0.9281, -0.3781, 0.0461, 0.0040, 1.0578, 0.1835, 
0.9436, -0.1237, -1.4173, -0.3306, 0.0424, .0003, 1.0805 , 0.2578, 
1.0636, -0.1388, -1.9976, -0.2795, 0.0364, -.0003, 1.0714, 0.3171, 
1.1370, -0.1494, -2.4200, -0.2670, 0.0360, -.0024, 1.08994, 0.3809, 
1.2213, -0.1619, -2.8644, -0.2632, 0.0368, -.0040, 1.1039, 0.4492, 
1.3111, -0.1751, -3.3392, -0.2580, 0.0360, -.0039, 1.1173, 0.5309, 
1.4265, -0.1925, -3.8950, -0.2544, 0.0373, -.0064, 1.1421, 0.5734, 
1.4818, -0.1986, -4.2420, -0.2111, 0.0329, -.0079, 1.1229, 0.6460, 
1.5802, -0.2134, -4.7384, -0.1964, 0.0310, -.0083, 1.1320),ncol=8, byrow=TRUE)

calclim<-p
if (p>10) calclim=10
coefsred<-coefs[1:calclim,]
temp<-c(p:1)
stick<-sort(cumsum(1/temp),decreasing=TRUE)[1:calclim]
multipliers<-matrix(c(log(n),log(p),log(n)*log(p),1),nrow=1)
longman<-exp(multipliers%*%t(coefs[,1:4]))[1:calclim]

allen<-rep(NA,calclim)
leig0<-0
newlim=calclim
if (calclim+2<p) newlim=newlim+2
for ( i in 1:(newlim-2)){
    leig1<-coefsred[i,5:8]%*%matrix(c(1,log(n-1),log((p-i-1)*(p-i+2)/2),leig0))
    leig0<-leig1
    allen[i]<-exp(leig1)
}
pcompnum<-c(1:calclim)
data.frame(cbind(pcompnum,stick,longman,allen))

}


#########
#this function makes a nice plot if given the input from a PCA analysis
#created by princomp()
##
#arguments are
#    n=number of observations

parallelplot=function(comp){
 if (comp$n.obs>1000 || length(comp$sdev)>100) {
    print ("Sorry, this only works for n<1000 and p<100")
    stop()
 }
#if (round(length(comp$sdev)) < round(sum(comp$sdev^2))) {
#    print ("Sorry, this only works for analyses using the correlation matrix")
#    stop()
# }

parallelanal=parallel(comp$n.obs,length(comp$sdev))
calclim=min(10,length(comp$sdev))
eigenvalues=(comp$sdev^2)[1:calclim]
limits=as.matrix(parallelanal[,2:4])
limits=limits[complete.cases(limits)]
ymax=range(c(eigenvalues),limits)
plot(parallelanal$pcompnum,eigenvalues,xlab="Principal Component Number",ylim=c(ymax),ylab="Eigenvalues and Thresholds",main="Scree plot with Parallel Analysis Limits",type="b",pch=15,lwd=2)
lines(parallelanal$pcompnum,parallelanal[,2], type="b",col="red",pch=16,lwd=2)
lines(parallelanal$pcompnum,parallelanal[,3], type="b",col="green",pch=17,lwd=2)
lines(parallelanal$pcompnum,parallelanal[,4], type="b",col="blue",pch=18,lwd=2)
legend((calclim/2),ymax[2],legend=c("Eigen-values","Stick Method","Longman Method","Allen Method"),  pch=c(15:18), col=c("black","red","green","blue"),lwd=2)

}





