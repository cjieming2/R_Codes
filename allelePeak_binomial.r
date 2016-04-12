setwd("C:/Users/JM/thesis/mark_work/allele_specificity/alleleSeqPeak")

ctfile = 'counts.aPeak'
data = read.delim(ctfile,header=T,sep='\t',row.names=1,stringsAsFactors=F)
p = 0.5
results = data.frame(row.names=rownames(data),
                     'cMat'=numeric(nrow(data)),
                     'cPat'=numeric(nrow(data)),
                     'wins'=character(nrow(data)),
                     'snpid'=character(nrow(data)),
                     'snpcounts'=character(nrow(data)),
                     'p.value'=numeric(nrow(data)),
                     'FDR'=numeric(nrow(data)),'BH'=numeric(nrow(data)),
                     stringsAsFactors=F);

for (i in 1:nrow(data))
{
  n = data[i,]$cMat + data[i,]$cPat
  
  if(n == 0)
  {
    results[i,] = data.frame(data[i,],'NA','NA')
    next;
  }
  else
  {
    k = data[i,]$cPat
    x = binom.test(k,n,p,"two.sided");
    results[i,] = data.frame(data[i,],x$p.value,0)
  }
}
results = results[order(results$p.value),]

for (i in 1:nrow(data))
{
  pval.fwer = results$p.value[i]*nrow(data)
  results$FDR[i] = min(pval.fwer/i,1) ## minimum of qvalue and 1 (ensure not >1)
}

results$BH = p.adjust(results$p.value, method = "fdr")

write.table(results, file = "counts.aPeak.fdr",sep="\t",quote=F)