setwd("C:/Users/JM/thesis/mark_work/allele_specificity/tf-coassoc")
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
## this function plots the ratios and prints to pdfs
## prints Spearman correlation

ratioPlot <- function(tflist, mybiglist, filename, numDataPointsCutOff,
                      indiv1, indiv2)
{
  # initialize count matrix
  mymatrix = matrix(0,length(tflist),length(tflist));
  dimnames(mymatrix) = list(tflist,tflist);
  
  # initialize spearman correlation matrix
  #mycorrel = matrix(0,length(tflist),length(tflist));
  #dimnames(mycorrel) = list(tflist,tflist);
  
  for (i in 1:length(names(mybiglist)))
  {
    tfpair = names(mybiglist)[i];
    tf = unlist(strsplit(tfpair, "-"));
    name=paste(tfpair,"-",filename,sep='');
    
    
    ## these are not symmetric matrices since x and y are also composites indiv+TFs
    ## combine counts and correlation matrices
    mycounts = nrow(mybiglist[[tfpair]]);
    mycorrel = cor(mybiglist[[tfpair]][,1],mybiglist[[tfpair]][,2],
                   method="spearman");
    #     mycortest= cor.test(mybiglist[[tfpair]][,1],mybiglist[[tfpair]][,2],
    #                         alternative="two.sided",
    #                         method="spearman")
    
    #     myentry  = paste(mycounts,";",mycorrel,";",mycortest$p.value);
    myentry  = paste(mycounts,";",mycorrel);
    
    mymatrix[tf[1],tf[2]] = myentry;
    
    #     mymatrix[tf[1],tf[2]] = nrow(mybiglist[[tfpair]]);
    #     mycorrel[tf[1],tf[2]] = cor(mybiglist[[tfpair]][,1],mybiglist[[tfpair]][,2],
    #                                 method="spearman");
    #     cortest               = cor.test(mybiglist[[tfpair]][,1],mybiglist[[tfpair]][,2],
    #                                      alternative="two.sided",
    #                                      method="spearman")
    
    
    if(mycounts >= numDataPointsCutOff)
    {
      pdf(paste(name,".pdf",sep=''));
      plot(mybiglist[[tfpair]], main=name, 
           xlab=paste('ref/total for',indiv1,'for',tf[1],sep=' '),
           ylab=paste('ref/total for',indiv2,'for',tf[2],sep=' '),
           ylim=c(0,1),xlim=c(0,1));
      #       text(0.3,0.75,paste("p val=",signif(cortest$p.value,3)))
      text(0.3,0.8,paste("Spearman's cor=",signif(mycorrel,3)))
      text(0.3,0.85,paste("numPoints=",mycounts))
      dev.off();
    }
    
  }
  
  ## print count matrix
  write.table(mymatrix, file=paste("matrix-",filename,".txt",sep=''), 
              col.names=TRUE, row.names=TRUE,
              sep="\t")
  
  #   ## print correlation matrix
  #   write.table(mycorrel, file=paste("correl-",filename,".txt",sep=''), 
  #               col.names=TRUE, row.names=TRUE,
  #               sep="\t")
}

#####################################################################
## this gives you just the correlation matrix
## prints Spearman correlation

corMat <- function(tflist, mybiglist, filename, numDataPointsCutOff,
                   indiv1, indiv2)
{
  # initialize count matrix
  mymatrix = matrix(0,length(tflist),length(tflist));
  dimnames(mymatrix) = list(tflist,tflist);
  
  # initialize spearman correlation matrix
  mycorrel = matrix(0,length(tflist),length(tflist));
  dimnames(mycorrel) = list(tflist,tflist);
  
  for (i in 1:length(names(mybiglist)))
  {
    tfpair = names(mybiglist)[i];
    tf = unlist(strsplit(tfpair, "-"));
    name=paste(tfpair,"-",filename,sep='');
    
    
    ## these are not symmetric matrices since x and y are also composites indiv+TFs
    mycounts = nrow(mybiglist[[tfpair]]);
    
    if(mycounts > numDataPointsCutOff)
    {
      
      mycorrel[tf[1],tf[2]] = cor(mybiglist[[tfpair]][,1],mybiglist[[tfpair]][,2],
                                  method="spearman");
      cortest               = cor.test(mybiglist[[tfpair]][,1],mybiglist[[tfpair]][,2],
                                       alternative="two.sided",
                                       method="spearman")
      mymatrix[tf[1],tf[2]] = paste(mycounts,";",cortest$p.value)
      #       mymatrix[tf[1],tf[2]] = paste(mycounts)
    }
    
  }
  
  
  ## print count matrix
  write.table(mymatrix, file=paste("zmatrix-",filename,".txt",sep=''), 
              col.names=TRUE, row.names=TRUE,
              sep="\t")
  
  ## print correlation matrix
  write.table(mycorrel, file=paste("zcorrel-",filename,".txt",sep=''), 
              col.names=TRUE, row.names=TRUE,
              sep="\t")
}

#####################################################################
#####################################################################
## this function plots the ratios and prints to pdfs
## prints Spearman correlation FOR one plot

ratioPlotOne <- function(mybiglist, filename, indiv1, indiv2, tfpair)
{  
  tf = unlist(strsplit(tfpair, "-"));
  name=paste(tfpair,"-",filename,sep='');
  
  
  counts = nrow(mybiglist[[tfpair]]);
  
  correl = cor(mybiglist[[tfpair]][,1],mybiglist[[tfpair]][,2],
               method="spearman");
  cortest = cor.test(mybiglist[[tfpair]][,1],mybiglist[[tfpair]][,2],
                     alternative="two.sided",
                     method="spearman")  
  
  
  
  
  pdf(paste(name,".pdf",sep=''));
  plot(mybiglist[[tfpair]], main=name, 
       xlab=paste('ref/total for',indiv1,'for',tf[1],sep=' '),
       ylab=paste('ref/total for',indiv2,'for',tf[2],sep=' '),
       ylim=c(0,1),xlim=c(0,1));
  #     text(0.5,0.45,paste("p val=",signif(cortest$p.value,3)))
  #     text(0.5,0.5,paste("Spearman's cor=",signif(correl,3)))
  #     text(0.5,0.55,paste("numPoints=",counts))
  
  text(0.3,0.75,paste("p val=",signif(cortest$p.value,3)))
  text(0.3,0.8,paste("Spearman's cor=",signif(correl,3)))
  text(0.3,0.85,paste("numPoints=",counts))
  dev.off();
  
  
}

#####################################################################
## this function plots the ratios and prints to pdfs
## prints Spearman correlation and homo 00_11 colored

ratioColoredPlot <- function(tflist, mybiglist, filename, numDataPointsCutOff,
                             indiv1, indiv2, color00, color11)
{
  
  for (i in 1:length(names(mybiglist)))
  {
    tfpair = names(mybiglist)[i];
    tf = unlist(strsplit(tfpair, "-"));
    name=paste(tfpair,"-",filename,sep='');
    
    
    ## these are not symmetric matrices since x and y are also composites indiv+TFs
    ## combine counts and correlation matrices
    mycounts = nrow(mybiglist[[tfpair]]);
    mycorrel = cor(mybiglist[[tfpair]][,1],mybiglist[[tfpair]][,2],
                   method="spearman");
    
    if(mycounts >= numDataPointsCutOff)
    {
      pdf(paste(name,"-colored.pdf",sep=''));
      plot(mybiglist[[tfpair]], main=name, 
           xlab=paste('ref/total for',indiv1,'for',tf[1],sep=' '),
           ylab=paste('ref/total for',indiv2,'for',tf[2],sep=' '),
           ylim=c(0,1),xlim=c(0,1));
      points(color00[[tfpair]],col='red');
      points(color11[[tfpair]],col='blue');
      
      text(0.3,0.8,paste("Spearman's cor=",signif(mycorrel,3)))
      text(0.3,0.85,paste("numPoints=",mycounts))
      dev.off();
    }
    
  }
  
}

#####################################################################

# within individual 78
# na12878.78.file <- "ratio.allHets.12878.12878.asb.bed";
na12878.78.file <- "ratio.2Hets.12878.12878.asb.bed";

na12878.78.data <- read.table(na12878.78.file, header=F, sep = "\t", stringsAsFactors = F)
uniq.tf.pairs.78.78 = unique(paste(na12878.78.data[,5],na12878.78.data[,7],sep="-"));

newlist.78.78   <- ratioParse(uniq.tf.pairs.78.78, na12878.78.data);
mybiglist.78.78 <- newlist.78.78[['mybiglist']];
tflist.78.78    <- newlist.78.78[['tflist']];

numDataPointsCutOffcore = 10
corMat(tflist.78.78, mybiglist.78.78, na12878.78.file, numDataPointsCutOffcore,
       'na12878', 'na12891')

numDataPointsCutOff = 5
ratioPlot(tflist.78.78, mybiglist.78.78, na12878.78.file, numDataPointsCutOff,
          'na12878', 'na12878')

## plot everything in one plot
pdf(paste('ALL-',na12878.78.file,'.pdf',sep=''))
plot(na12878.78.data[,6],na12878.78.data[,8],
     main=paste('ALL ',na12878.78.file),
     xlab='ref/total for na12878',
     ylab='ref/total for na12878')
dev.off();