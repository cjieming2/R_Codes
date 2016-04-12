#The Input is the gct file and the cls and cls2 come from the cls file
#The genes, p-values, q-values, and fwer values are output in order to the terminal

Input <- read.delim("all_aml_train.gct", header=TRUE,
                  sep="\t", skip=2,                    	
                  blank.lines.skip=TRUE)

cls<-scan("all_aml_train.cls", nlines = 1)                  
cls2<- scan("all_aml_train.cls", skip = 2, nlines =1)

dims <- dim(Input)
M <- dims[1]
N <- dims[2]
n<-cls[1]
sk <- N - n

##Determine where tumor starts
for (i in 1:n)
{
	if(cls2[i]==1)
	{
		tumor <- as.numeric(i)
		break
	}
}

A <- tumor-1
B <- tumor

#allocate array for T-test p-values
ttest <- array(0,dim=c(M))

#x: no cancer, y: cancer
x<- Input[,(sk+1):(A+sk)]
y<-Input[,(B+sk):(n+sk)]

#compute p-values for t-test
for (i in 1:M)
{
	ttest[i]<- t.test(x[i,],y[i,])$p.value
}

#Find order of t-tests
porder<-order(abs(ttest))

#compute Q value
#citation("qvalue")
source("http://bioconductor.org/biocLite.R")
biocLite("qvalue")
library(qvalue)
qvals<-qvalue(abs(ttest))

#Compute FWER
fwer <- array(0,dim=c(M))
fwer<-p.adjust(ttest, method = "bonferroni")

#print gene names then stats
print("Top 20 Genes (according to p-value): ")
for (i in 1:20)
{
	print(Input$Name[porder[i]], max.levels =0)
}
print("P-values: ")
for (i in 1:20)
{
	print(ttest[porder[i]])
}
print("Q-values: ")
for (i in 1:20)
{
	print(qvals$qvalue[porder[i]])
}
print("FWER values: ")
for (i in 1:20)
{
	print(fwer[porder[i]])
}
