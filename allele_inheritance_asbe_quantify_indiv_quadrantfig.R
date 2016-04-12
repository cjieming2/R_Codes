setwd("C:/Users/Jieming/Documents/thesis/mark_work/allele_specificity/compare trios ratios inheritance")

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
na12878.91.file <- "allelicRatio.CTCF.trio.3_2intHets.inherit.78_91.asb.bed";
na12878.92.file <- "allelicRatio.CTCF.trio.3_2intHets.inherit.78_92.asb.bed";
na12891.92.file <- "allelicRatio.CTCF.trio.3_2intHets.inherit.91_92.asb.bed";

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

## manual; PU.1=>PU1; Pol2=>POL2
TF    = 'SA1' ; comp2 = common.91.92
# length(which(comp2$V5 == TF & comp2$V7 == TF))
qa=length(which(comp2$V5 == TF & comp2$V7 == TF & comp2$V6 <= 0.5 & comp2$V8 >= 0.5)) # quadA
qb=length(which(comp2$V5 == TF & comp2$V7 == TF & comp2$V6 <= 0.5 & comp2$V8 <= 0.5)) # quadB
qc=length(which(comp2$V5 == TF & comp2$V7 == TF & comp2$V6 >= 0.5 & comp2$V8 >= 0.5)) # quadC
qd=length(which(comp2$V5 == TF & comp2$V7 == TF & comp2$V6 >= 0.5 & comp2$V8 <= 0.5)) # quadD
matrix(c(qa,qb,qc,qd),ncol=2)

# x=binom.test(192,199,0.5,"two.sided");  x$p.value
# x=binom.test(237,250,0.5,"two.sided");	x$p.value
# x=binom.test(110,113,0.5,"two.sided");	x$p.value
# x=binom.test(56,77,0.5,"two.sided");	x$p.value
# x=binom.test(84,112,0.5,"two.sided");	x$p.value
# x=binom.test(26,39,0.5,"two.sided");	x$p.value
# x=binom.test(251,258,0.5,"two.sided");	x$p.value
# x=binom.test(146,153,0.5,"two.sided");	x$p.value
# x=binom.test(93,96,0.5,"two.sided");	x$p.value
# x=binom.test(28,31,0.5,"two.sided");	x$p.value
# x=binom.test(68,75,0.5,"two.sided");	x$p.value
# x=binom.test(9,10,0.5,"two.sided");	x$p.value
# x=binom.test(53,60,0.5,"two.sided");	x$p.value
# x=binom.test(54,63,0.5,"two.sided");	x$p.value
# x=binom.test(38,52,0.5,"two.sided");	x$p.value
# x=binom.test(441,800,0.5,"two.sided");	x$p.value
# x=binom.test(450,767,0.5,"two.sided");	x$p.value
# x=binom.test(315,577,0.5,"two.sided");	x$p.value
# x=binom.test(1051,1315,0.5,"two.sided");	x$p.value
# x=binom.test(1028,1267,0.5,"two.sided");	x$p.value
# x=binom.test(1068,1479,0.5,"two.sided");	x$p.value
