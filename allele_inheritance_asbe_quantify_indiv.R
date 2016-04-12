setwd("C:/Users/JM/thesis/mark_work/allele_specificity/compare trios ratios inheritance/asb")

#####################################################################
## this function parses the data into ratios for plotting

ratioParse <- function(uniq.tf.pairs, data)
{
  mybiglist <- list();
  tflist    <- vector();

  for (i in 1:length(uniq.tf.pairs))
  {
    # split the pair
    tf = unlist(strsplit(uniq.tf.pairs[i], "-"));
    tmp = c();
    
    # find this pair in the data file
    for (j in 1:nrow(data))
    {
      if(tf[1] == data[j,5] & tf[2] == data[j,7])
      {
        tmp = rbind(c(data[j,6],data[j,8]),tmp);
      }
    }
    mybiglist[[uniq.tf.pairs[i]]] = tmp;
    tflist = rbind(tf[1],tf[2],tflist);
  }
  
  tflist  = unique(tflist);
  newList = list("mybiglist"=mybiglist, "tflist"=tflist);
  return(newList);
}


#####################################################################

## parameters
qnumDataPointsCutOff = 10;


## hets 2 and 3 combined
# title = 'combined'
# na12878.91.file <- "ratio.3_2Hets.12878.12891.ase.bed";
# na12878.92.file <- "ratio.3_2Hets.12878.12892.ase.bed";
# na12891.92.file <- "ratio.3_2Hets.12891.12892.ase.bed";

title = 'combined'
na12878.91.file <- "ratio.3_2Hets.12878.12891.asb.bed";
na12878.92.file <- "ratio.3_2Hets.12878.12892.asb.bed";
na12891.92.file <- "ratio.3_2Hets.12891.12892.asb.bed";

# title = 'test'
# na12878.91.file <- "test.78.91.bed";
# na12878.92.file <- "test.78.92.bed";
# na12891.92.file <- "test.91.92.bed";


## hets snp.calls chr-pos chr pos refallele mat pat child phase
na12878.91.data <- read.table(na12878.91.file, header=F, sep = "\t", stringsAsFactors = F)
na12878.92.data <- read.table(na12878.92.file, header=F, sep = "\t", stringsAsFactors = F)
na12891.92.data <- read.table(na12891.92.file, header=F, sep = "\t", stringsAsFactors = F)


## all TF names to upper case and PU.1 to PU1
na12878.91.data <- data.frame(lapply(na12878.91.data, function(v) 
                        { if (is.character(v)) return(toupper(v))
                          else return(v) }),stringsAsFactors=FALSE)
na12878.92.data <- data.frame(lapply(na12878.92.data, function(v) 
                        { if (is.character(v)) return(toupper(v))
                          else return(v) }),stringsAsFactors=FALSE)
na12891.92.data <- data.frame(lapply(na12891.92.data, function(v) 
                        { if (is.character(v)) return(toupper(v))
                          else return(v) }),stringsAsFactors=FALSE)
na12878.91.data <- data.frame(lapply(na12878.91.data, function(v) 
                        { if (is.character(v)|is.factor(v)) gsub("PU.1","PU1",v) else v }), stringsAsFactors=FALSE)
na12878.92.data <- data.frame(lapply(na12878.92.data, function(v) 
                        { if (is.character(v)|is.factor(v)) gsub("PU.1","PU1",v) else v }), stringsAsFactors=FALSE)
na12891.92.data <- data.frame(lapply(na12891.92.data, function(v) 
                        { if (is.character(v)|is.factor(v)) gsub("PU.1","PU1",v) else v }), stringsAsFactors=FALSE)


## get only common TFs
common.78.91 <- subset(na12878.91.data, V5 == V7)
common.78.92 <- subset(na12878.92.data, V5 == V7)
common.91.92 <- subset(na12891.92.data, V5 == V7)


## find unique tf pairs in common TFs
uniq.tf.pairs.78.91 = unique(paste(common.78.91[,5],common.78.91[,7],sep="-"));
uniq.tf.pairs.78.92 = unique(paste(common.78.92[,5],common.78.92[,7],sep="-"));
uniq.tf.pairs.91.92 = unique(paste(common.91.92[,5],common.91.92[,7],sep="-"));



## create matrices for each pair
newlist.78.91   <- ratioParse(uniq.tf.pairs.78.91, common.78.91);
mybiglist.78.91 <- newlist.78.91[['mybiglist']];
tflist.78.91    <- newlist.78.91[['tflist']];

newlist.78.92 <- ratioParse(uniq.tf.pairs.78.92, common.78.92);
mybiglist.78.92 <- newlist.78.92[['mybiglist']];
tflist.78.92    <- newlist.78.92[['tflist']];

newlist.91.92 <- ratioParse(uniq.tf.pairs.91.92, common.91.92);
mybiglist.91.92 <- newlist.91.92[['mybiglist']];
tflist.91.92    <- newlist.91.92[['tflist']];


## split into 4 quadrants
# 1) add up parent-child take away parent-parent
# 2) fisher's test
# 3) ratio of parent-child to parent-parent

## for each TF (CEU) ; do not APPEND for this
write(paste("TF","diff","diff.t","diff.78.91","diff.78.92",
            "ratio","ratio.t","ratio.r","ratio.78.91","ratio.78.92",
            "fisher.pval","fisher.OR",
              sep="\t"),"asb_inheritance_quantify.txt",sep="\t",append = FALSE)

for (i in 1:length(names(mybiglist.78.91)))
{
  tfpair = names(mybiglist.78.91)[i];
  tf = unlist(strsplit(tfpair, "-"));
  name=paste(tfpair,"-","78.91",sep='');
  
  if(tfpair == "POU2F2-POU2F2"){ next;  }
  
  mybiglist.p2c = rbind(mybiglist.78.91[[tfpair]],mybiglist.78.92[[tfpair]])
  mybiglist.p2p = rbind(mybiglist.91.92[[tfpair]])
  
  ## count total for each TF
  total.78.91 = nrow(subset(common.78.91, 
                            common.78.91$V5 == tf[1] & common.78.91$V7 == tf[2]))
  total.78.92 = nrow(subset(common.78.92, 
                            common.78.92$V5 == tf[1] & common.78.92$V7 == tf[2]))
  total.91.92 = nrow(subset(common.91.92,
                            common.91.92$V5 == tf[1] & common.91.92$V7 == tf[2]))
  total       = total.78.91 + total.78.92 + total.91.92
  
  
  ## parent2child p2c quads
  quadA.p2c = subset(mybiglist.p2c,
                   (mybiglist.p2c[,1]<=0.5) & (mybiglist.p2c[,2]>=0.5))
  quadB.p2c = subset(mybiglist.p2c,
                   (mybiglist.p2c[,1]>=0.5) & (mybiglist.p2c[,2]>=0.5))
  quadC.p2c = subset(mybiglist.p2c,
                   (mybiglist.p2c[,1]<=0.5) & (mybiglist.p2c[,2]<=0.5))
  quadD.p2c = subset(mybiglist.p2c,
                   (mybiglist.p2c[,1]>=0.5) & (mybiglist.p2c[,2]<=0.5))
  
  quadA.p2c.78.91 = subset(mybiglist.78.91[[tfpair]],
                     (mybiglist.78.91[[tfpair]][,1]<=0.5) & (mybiglist.78.91[[tfpair]][,2]>=0.5))
  quadB.p2c.78.91 = subset(mybiglist.78.91[[tfpair]],
                     (mybiglist.78.91[[tfpair]][,1]>=0.5) & (mybiglist.78.91[[tfpair]][,2]>=0.5))
  quadC.p2c.78.91 = subset(mybiglist.78.91[[tfpair]],
                     (mybiglist.78.91[[tfpair]][,1]<=0.5) & (mybiglist.78.91[[tfpair]][,2]<=0.5))
  quadD.p2c.78.91 = subset(mybiglist.78.91[[tfpair]],
                     (mybiglist.78.91[[tfpair]][,1]>=0.5) & (mybiglist.78.91[[tfpair]][,2]<=0.5))
  
  quadA.p2c.78.92 = subset(mybiglist.78.92[[tfpair]],
                           (mybiglist.78.92[[tfpair]][,1]<=0.5) & (mybiglist.78.92[[tfpair]][,2]>=0.5))
  quadB.p2c.78.92 = subset(mybiglist.78.92[[tfpair]],
                           (mybiglist.78.92[[tfpair]][,1]>=0.5) & (mybiglist.78.92[[tfpair]][,2]>=0.5))
  quadC.p2c.78.92 = subset(mybiglist.78.92[[tfpair]],
                           (mybiglist.78.92[[tfpair]][,1]<=0.5) & (mybiglist.78.92[[tfpair]][,2]<=0.5))
  quadD.p2c.78.92 = subset(mybiglist.78.92[[tfpair]],
                           (mybiglist.78.92[[tfpair]][,1]>=0.5) & (mybiglist.78.92[[tfpair]][,2]<=0.5))
  
  
  ## parent2parent p2p quads
  quadA.p2p = subset(mybiglist.p2p, 
                     (mybiglist.p2p[,1]<=0.5) & (mybiglist.p2p[,2]>=0.5))
  quadB.p2p = subset(mybiglist.p2p, 
                     (mybiglist.p2p[,1]>=0.5) & (mybiglist.p2p[,2]>=0.5))
  quadC.p2p = subset(mybiglist.p2p, 
                     (mybiglist.p2p[,1]<=0.5) & (mybiglist.p2p[,2]<=0.5))
  quadD.p2p = subset(mybiglist.p2p, 
                     (mybiglist.p2p[,1]>=0.5) & (mybiglist.p2p[,2]<=0.5))
  
  ## quads
  quadA = (nrow(quadA.p2c)-nrow(quadA.p2p))
  quadB = (nrow(quadB.p2c)-nrow(quadB.p2p))
  quadC = (nrow(quadC.p2c)-nrow(quadC.p2p))
  quadD = (nrow(quadD.p2c)-nrow(quadD.p2p))
  
  quadA.78.91 = (nrow(quadA.p2c.78.91)-nrow(quadA.p2p))
  quadB.78.91 = (nrow(quadB.p2c.78.91)-nrow(quadB.p2p))
  quadC.78.91 = (nrow(quadC.p2c.78.91)-nrow(quadC.p2p))
  quadD.78.91 = (nrow(quadD.p2c.78.91)-nrow(quadD.p2p))
  
  quadA.78.92 = (nrow(quadA.p2c.78.92)-nrow(quadA.p2p))
  quadB.78.92 = (nrow(quadB.p2c.78.92)-nrow(quadB.p2p))
  quadC.78.92 = (nrow(quadC.p2c.78.92)-nrow(quadC.p2p))
  quadD.78.92 = (nrow(quadD.p2c.78.92)-nrow(quadD.p2p))
  
  ## quantity 1 : difference of differences
  q1   = (( quadB + quadC ) - ( quadA + quadD ))
  q1.t = (( quadB + quadC ) - ( quadA + quadD )) / total
  q1.78.91 = (( quadB.78.91 + quadC.78.91 ) - ( quadA.78.91 + quadD.78.91 )) / total.78.91
  q1.78.92 = (( quadB.78.92 + quadC.78.92 ) - ( quadA.78.92 + quadD.78.92 )) / total.78.92
  
  
  ## quantity 2 : ratio of differences and ratio of ratios (r)
  q2   = round(( quadB + quadC ) / ( quadA + quadD + 0.0001 ))
  q2.t = 1 - 1/(round(( quadB + quadC ) / ( quadA + quadD + 0.0001 )))
  q2.r = round(( nrow(quadB.p2c) + nrow(quadC.p2c) ) / ( nrow(quadA.p2c) + nrow(quadD.p2c) + 0.0001 )) /
            round(( nrow(quadB.p2p) + nrow(quadC.p2p) ) / ( nrow(quadA.p2p) + nrow(quadD.p2p) + 0.0001 ))
  q2.r.78.91 = round(( nrow(quadB.p2c.78.91) + nrow(quadC.p2c.78.91) ) / ( nrow(quadA.p2c.78.91) + nrow(quadD.p2c.78.91) + 0.0001 )) /
            round(( nrow(quadB.p2p) + nrow(quadC.p2p) ) / ( nrow(quadA.p2p) + nrow(quadD.p2p) + 0.0001 ))
  q2.r.78.92 = round(( nrow(quadB.p2c.78.92) + nrow(quadC.p2c.78.92) ) / ( nrow(quadA.p2c.78.92) + nrow(quadD.p2c.78.92) + 0.0001 )) /
            round(( nrow(quadB.p2p) + nrow(quadC.p2p) ) / ( nrow(quadA.p2p) + nrow(quadD.p2p) + 0.0001 ))
  
  ## quantity 3 : hypergeometric direct difference of p2c and p2p
  q3 = matrix(c( abs( quadB+quadC ) , abs( quadA+quadC ) , 
                 abs( quadB+quadD ) , abs( quadA+quadD )), 2,2,
                 dimnames = list(upper = c("B","A"), lower = c("C","D")))
    
  # default is 2-sided, 95% CI
  q3.x = fisher.test(q3,alternative="two.sided")
  
  write(paste(tfpair,q1,q1.t,q1.78.91,q1.78.92,
                     q2,q2.t,q2.r,q2.r.78.91,q2.r.78.92,
                     q3.x$p.value,q3.x$estimate,sep="\t"),"asb_inheritance_quantify.txt",sep="\t",append = TRUE)
}
