

### fdr from description from https://support.bioconductor.org/p/18967/
### why getting uniform adjusted values despite ranking
# In the BH (fdr) adjustment, you
# basically take the smallest raw p-value and
# multiply it by n, the number of genes, to get a
# temporary adjusted p-value, then you take the
# second smallest raw p and multiply it by n/rank,
# etc., up the line creating temporary adjusted
# p-values.  When a higher ranking temporary
# p-value is lower than ones below it, all the
# p-values below it get replaced with that value as
# the final adjusted p-value. The largest raw
# p-value is only multiplied by 1 and in this data
# set of yours,  this is the lowest temporary
# adjusted p-value for ALL the genes, which is why
# you're seeing only one uniform "adjusted p-value".



# load cls file from working directory, then same with res file
cls = read.delim("all_aml_train.cls", header=F,skip=2,sep=" ")
cls = cls[,-ncol(cls)]
cls = sapply(cls,function(x) { if(x==1) {"AML"} else {"ALL"}})
data = read.delim("all_aml_train.res", header=F,quote="",row.names=2,skip=3)
#
# throw out the description and absent/present flag
# transpose data to make each gene a column and add column for cls
data = data[,-seq(from=1, to=77, by=2)]
tdata = data.frame(t(data))
data = cbind(cls,tdata)
colnames(data)[1] = "cls"
#
# split data into two frames, one for ALL, one for AML
data.all = data[data$cls=="ALL",-1]
data.aml = data[data$cls=="AML",-1]
#
# store the number of genes
ngenes = ncol(data.all)
#
# create results data frame that will store (adjusted) p-values
results = data.frame("p.value"=numeric(ngenes),
"FWER"=numeric(ngenes),"FDR" = numeric(ngenes))
rownames(results) = colnames(data.all)
#
# calculate two-sided t-test p-value for each gene
for(i in 1:ngenes) {
results$p.value[i] = t.test(data.all[,i],data.aml[,i])$p.value
}
#
# sort results data-frame by p-value
results = results[order(results$p.value),]
#
# calculate adjusted p-values to control FWER, then FDR
for(i in 1:ngenes) {
pval.fwer = results$p.value[i]*ngenes
results$FWER[i] = min(pval.fwer,1)
results$FDR[i] = min(pval.fwer/i,1)
}
#
results[1:20,]