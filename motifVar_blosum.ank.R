setwd("C:/Users/Jieming/Documents/thesis/lynne_work/motifVar/ank/blosum_sift_polyphen2_gerp")

library(vioplot)
library(gplots)

### modified vioplot las=2 function ####
### this a modified function of vioplot from the package vioplot
### the only difference is las=2, so the x labels are vertical
vioplot.las2 <- function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, 
          horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
          lwd = 1, rectCol = "black", colMed = "white", pchMed = 19, 
          at, add = FALSE, wex = 1, drawRect = TRUE) 
{
  datas <- list(x, ...)
  n <- length(datas)
  if (missing(at)) 
    at <- 1:n
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)
  args <- list(display = "none")
  if (!(is.null(h))) 
    args <- c(args, h = h)
  for (i in 1:n) {
    data <- datas[[i]]
    data.min <- min(data)
    data.max <- max(data)
    q1[i] <- quantile(data, 0.25)
    q3[i] <- quantile(data, 0.75)
    med[i] <- median(data)
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, data.max)
    lower[i] <- max(q1[i] - range * iqd, data.min)
    est.xlim <- c(min(lower[i], data.min), max(upper[i], 
                                               data.max))
    smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
                                     args))
    hscale <- 0.4/max(smout$estimate) * wex
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
  }
  if (!add) {
    xlim <- if (n == 1) 
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  }
  else {
    label <- names
  }
  boxwidth <- 0.05 * wex
  if (!add) 
    plot.new()
  if (!horizontal) {
    if (!add) {
      plot.window(xlim = xlim, ylim = ylim)
      axis(2)
      axis(1, at = at, label = label, las=2)
    }
    box()
    for (i in 1:n) {
      polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
              c(base[[i]], rev(base[[i]])), col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
              lty = lty)
        rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, 
             q3[i], col = rectCol)
        points(at[i], med[i], pch = pchMed, col = colMed)
      }
    }
  }
  else {
    if (!add) {
      plot.window(xlim = ylim, ylim = xlim)
      axis(1)
      axis(2, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], 
                                              rev(at[i] + height[[i]])), col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
              lty = lty)
        rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + 
               boxwidth/2, col = rectCol)
        points(med[i], at[i], pch = pchMed, col = colMed)
      }
    }
  }
  invisible(list(upper = upper, lower = lower, median = med, 
                 q1 = q1, q3 = q3))
}

### boxbubbleplot function ####
### this function takes in the title, residue number col, blosum score col
### and the last col
### make 2 combined plots (box plus bubble) of sequential and entropy-ordered positions

boxbubbleplot <- function(title,resNumCol,scoreCol,lastCol,data, cexAxisFontSize=2)
{
  title = title
  resNum = data[,resNumCol]
  blosum = data[,scoreCol]
  m=max(blosum)
  tt.blosum=table(resNum, blosum)
  range01 = function(x){(x-min(x))/(max(x)-min(x))}
# #   x11()
#   par(mar=c(5,6,4,4),xpd=TRUE,mar=c(1, 4, 1, 1) + 1)
#   boxplot(blosum~resNum,xlab="position",ylab=paste(title,sep=" "),cex.lab=2, cex.axis=cexAxisFontSize, las=2)
#   
#   ## convert table counts of SNVs at each position to dataframe
#   a.blosum=as.data.frame(table(resNum))[,2]
#   text(0,m+0.5,paste(a.blosum,collapse=" "),adj = c(0,0))
#   
#   ###### bubble plot
#   df = expand.grid(sort(unique(resNum)),sort(unique(blosum)))
#   aa = as.vector(table(resNum, blosum))
#   df$value = aa
#   
#   # x11()
#   area = pi * (df$value/150)^2
#   par(mar=c(5,6,4,4),xpd=TRUE)
#   area = (range01(df$value)/2)^2 * pi
#   symbols(df$Var1,df$Var2,circles=area, inches=F, bg=rgb(1,0,0,0.5), fg="black", yaxt="n", 
#           xaxt="n",add=T)
#   
#   # axis(1,at=seq(1,34), las=2)
#   # axis(2,at=seq(-4,9))
#   # text(0,3,paste(a.blosum,collapse=" "),adj = c(0,0))
  
  ## ordered
#   x11()
  par(mar=c(5,6,4,4),xpd=TRUE,mar=c(1, 4, 1, 1) + 1)
  boxplot(blosum~data[,lastcol],xlab="position",ylab=paste(title,sep=" "),cex.lab=2, cex.axis=cexAxisFontSize, las=2)
  a.blosumordered=as.data.frame(table(data[,lastcol]))[,2]
  text(0,m+0.5,paste(a.blosumordered,collapse=" "),adj = c(0,0))
  
  #### bubble
  # x11()
  par(mar=c(5,6,4,4),xpd=TRUE)
  df.ordered = expand.grid(sort(unique(data[,lastcol])),sort(unique(blosum)))
  aa.ordered = as.vector(table(data[,lastcol], blosum))
  df.ordered$value = aa.ordered
  area.ordered = (range01(df.ordered$value)/2)^2 * pi
  symbols(df.ordered$Var1,df.ordered$Var2,circles=area.ordered, inches=F, bg=rgb(1,0,0,0.5), fg="black", yaxt="n", 
          xaxt="n",xlab='position',ylab=paste(title, "scores",sep=" "),add=T)
  # axis(1,at=seq(1,34),labels=as.vector(sort(unique(data[,lastcol]))), las=2)
  # axis(2,at=seq(-4,9))
  
#   return(blosum);
}


#### main setup #####
filename = 'final.merged.ExAC.ens73.ANK.30aa.vepoutput.coding.canonical.noORF.sorted.auto.blosum-30-35-40-45-50-62.gerp.sift.bed'
data1 = read.table(filename, header=F, stringsAsFactors = F, comment.char="", sep="\t")

## remove AC=1 V8
data = data1[data1$V8 > 1,]

## create order
order = read.table('order_relentropy-ank.txt', header=T, stringsAsFactors = F, comment.char="", sep="\t")
lastcol=ncol(data)+1
resNumCol=45
## we insert ordered info into lastCol so order info is lastCol
data[,lastcol] = factor(data[,resNumCol], levels=levels(data[,resNumCol])<-order$pos)

## remove synonymous
data.ns.only <- data[(data$V56 == "NON_SYNONYMOUS_CODING" | data$V56 == "NON_SYNONYMOUS_CODING,SPLICE_SITE"),]


#### BLOSUMs ####
maxResNum = max(data[,resNumCol])
blos30 = matrix(list(NULL), c(maxResNum,1))
blos30.ns = matrix(list(NULL), c(maxResNum,1))
blos62.ns = matrix(list(NULL), c(maxResNum,1))

blos30col = 65
blos62col = 70

## create lookup
seqrange30 = as.data.frame(seq(min(data[,blos30col]),max(data[,blos30col])),row.names=seq(min(data[,blos30col]),max(data[,blos30col]))); names(seqrange30) = "blosum30"
seqrange30.ns = as.data.frame(seq(min(data.ns.only[,blos30col]),max(data.ns.only[,blos30col])),row.names=seq(min(data.ns.only[,blos30col]),max(data.ns.only[,blos30col]))); names(seqrange30.ns) = "blosum30"
# note that this is initializing according to blosum30
seqrange62.ns = as.data.frame(seq(min(data.ns.only[,blos30col]),max(data.ns.only[,blos30col])),row.names=seq(min(data.ns.only[,blos30col]),max(data.ns.only[,blos30col]))); names(seqrange62.ns) = "blosum62"
blos30.freq = matrix(0,ncol=maxResNum,nrow=nrow(seqrange30))
blos30.ns.freq = matrix(0,ncol=maxResNum,nrow=nrow(seqrange30.ns))
blos62.ns.freq = matrix(0,ncol=maxResNum,nrow=nrow(seqrange62.ns))

for(i in 1:maxResNum)
{  
  ## create freq table from each set of scores for each TPR pos i
  ## blosum30
  blos30[[i]] = data[,blos30col][data[,resNumCol] == i]
  temp30 = as.data.frame(prop.table(table(blos30[[i]])))
  names(temp30) = c("blosum30","freq") ## set the name to the same as seqrange30 for lookup
  temp30 = merge(seqrange30,temp30,by="blosum30",all.x=TRUE)
  temp30[is.na(temp30)] = 0
  temp30 = as.data.frame(temp30)
  blos30.freq[,i] = temp30$freq
  
  ## blosum30-NS
  blos30.ns[[i]] = data.ns.only[,blos30col][data.ns.only[,resNumCol] == i]
  temp30.ns = as.data.frame(prop.table(table(blos30.ns[[i]])))
  names(temp30.ns) = c("blosum30","freq") ## set the name to the same as seqrange30 for lookup
  temp30.ns = merge(seqrange30.ns,temp30.ns,by="blosum30",all.x=TRUE)
  temp30.ns[is.na(temp30.ns)] = 0
  temp30.ns = as.data.frame(temp30.ns)
  blos30.ns.freq[,i] = temp30.ns$freq
  
  ## blosum62-NS
  blos62.ns[[i]] = data.ns.only[,blos62col][data.ns.only[,resNumCol] == i]
  temp62.ns = as.data.frame(prop.table(table(blos62.ns[[i]])))
  names(temp62.ns) = c("blosum62","freq") ## set the name to the same as seqrange62 for lookup
  temp62.ns = merge(seqrange62.ns,temp62.ns,by="blosum62",all.x=TRUE)
  temp62.ns[is.na(temp62.ns)] = 0
  temp62.ns = as.data.frame(temp62.ns)
  blos62.ns.freq[,i] = temp62.ns$freq
}


## plot BLOSUM30

# heatmap (ordered)
# color map ; check min and max
# min(blos30.freq) ; max(blos30.freq)
cols = c(colorRampPalette(c("white","cornflowerblue"))(10), colorRampPalette(c("yellow","red"))(32))
bks = c(0,0,seq(0,0.5,by=0.0125)) ## only 0 will be white, the rest are split between 0.5/39 (since we have 39 colors - need to have 1 more break than color)
x11()
par(cex.axis=1,cex.lab=1.5)
heatmap.2(blos30.freq[c(nrow(blos30.freq):1),c(order[,1])],trace="none",key=TRUE,Rowv=NULL,Colv=NULL,dendrogram="none",
          labCol=order[,1],labRow=rev(seqrange30$blosum30), col=cols, breaks=bks,
          colsep=1:34,rowsep=1:25,sepcolor="white",cexRow=1.5,cexCol=1.5,
          lhei=c(1,5),lwid=c(1,5),
          key.title=NA,key.xlab=NA,key.ylab=NA,
          xlab="position", ylab="BLOSUM30 all", cex=1.5)

## plot BLOSUM30-NS ordered
x11()
heatmap.2(blos30.ns.freq[c(nrow(blos30.ns.freq):1),c(order[,1])],trace="none",key=TRUE,Rowv=NULL,Colv=NULL,dendrogram="none",
          labCol=order[,1],labRow=rev(seqrange30.ns$blosum30), col=cols, breaks=bks,
          colsep=1:34,rowsep=1:25,sepcolor="white",cexRow=1.5,cexCol=1.5,
          lhei=c(1,5),lwid=c(1,5),
          key.title=NA,key.xlab=NA,key.ylab=NA,
          xlab="position", ylab="BLOSUM30 NS", cex=1.5)

## plot BLOSUM62-NS ordered
x11()
heatmap.2(blos62.ns.freq[c(nrow(blos62.ns.freq):1),c(order[,1])],trace="none",key=TRUE,Rowv=NULL,Colv=NULL,dendrogram="none",
          labCol=order[,1],labRow=rev(seqrange62.ns$blosum62), col=cols, breaks=bks,
          colsep=1:34,rowsep=1:25,sepcolor="white",cexRow=1.5,cexCol=1.5,
          lhei=c(1,5),lwid=c(1,5),
          key.title=NA,key.xlab=NA,key.ylab=NA,
          xlab="position", ylab="BLOSUM62 NS", cex=1.5)

#vioplot
# blos30 = array(list(NULL), c(maxResNum,1))
# blos62 = array(list(NULL), c(maxResNum,1))
# for(i in 1:maxResNum)
# {  
#   blos30[[i]] = data[,blos30col][data[,resNumCol] == i]
#   blos62[[i]] = data[,blos62col][data[,resNumCol] == i] 
# }
# x11()
# do.call(vioplot,c(blos30,col="yellow"))
# mtext("position",side=1, line=2, cex=1.5)
# mtext("BLOSUM30",side=2, line=2, cex=1.5)
# 
# ## plot BLOSUM62
# x11()
# par(cex.axis=1,cex.lab=1.5)
# do.call(vioplot,c(blos62,col="yellow"))
# mtext("position",side=1, line=2, cex=1.5)
# mtext("BLOSUM62",side=2, line=2, cex=1.5)

## bubbleboxplot
# par(mfrow=c(4,1),tcl=-0.5)
# boxbubbleplot('BLOSUM30',resNumCol,77,lastCol,data,1.5)
# # boxbubbleplot('BLOSUM35',resNumCol,78,lastCol,data,1.5)
# boxbubbleplot('BLOSUM40',resNumCol,79,lastCol,data,1.5)
# # boxbubbleplot('BLOSUM45',resNumCol,80,lastCol,data,1.5)
# boxbubbleplot('BLOSUM50',resNumCol,81,lastCol,data,1.5)
# boxbubbleplot('BLOSUM62',resNumCol,66,lastCol,data,1.5)


## plot vioplots ####
blos30.ns = array(list(NULL), c(maxResNum,1))
blos62.ns = array(list(NULL), c(maxResNum,1))

for(i in 1:maxResNum)
{  
  blos30.ns[[i]] = data.ns.only[,blos30col][data.ns.only[,resNumCol] == i] 
  blos62.ns[[i]] = data.ns.only[,blos62col][data.ns.only[,resNumCol] == i] 
}



x11()
par(mfrow=c(2,2),tcl=-0.5)
par(cex.axis=1.5,cex.lab=1.5,mai=c(0.5,0.6,0.05,0.3))

## plot BLOSUM30
do.call(vioplot,c(blos30,col="yellow"))
mtext("position",side=1, line=2, cex=1.5)
mtext("BLOSUM30",side=2, line=2, cex=1.5)

## plot BLOSUM30-NS
do.call(vioplot,c(blos30.ns,col="orange"))
mtext("position",side=1, line=2, cex=1.5)
mtext("BLOSUM30-NS",side=2, line=2, cex=1.5)

## plot BLOSUM62
do.call(vioplot,c(blos62,col="yellow"))
mtext("position",side=1, line=2, cex=1.5)
mtext("BLOSUM62",side=2, line=2, cex=1.5)

## plot BLOSUM62-NS
do.call(vioplot,c(blos62.ns,col="orange"))
mtext("position",side=1, line=2, cex=1.5)
mtext("BLOSUM62-NS",side=2, line=2, cex=1.5)

## vioplots plot ORDERED #######
# blos30.ordered = array(list(NULL), c(maxResNum,1))
# blos62.ordered = array(list(NULL), c(maxResNum,1))
# blos30.ns.ordered = array(list(NULL), c(maxResNum,1))
# blos62.ns.ordered = array(list(NULL), c(maxResNum,1))
# ctr = 1
# for(i in order[,1])
# {
#   blos30.ordered[[ctr]] = data[,blos30col][data[,lastcol] == i] 
#   blos62.ordered[[ctr]] = data[,blos62col][data[,lastcol] == i]
#   blos30.ns.ordered[[ctr]] = data.ns.only[,blos30col][data.ns.only[,lastcol] == i] 
#   blos62.ns.ordered[[ctr]] = data.ns.only[,blos62col][data.ns.only[,lastcol] == i] 
#   ctr = ctr + 1
# }
# 
# x11()
# par(mfrow=c(2,2),tcl=-0.5)
# par(cex.axis=1.5,cex.lab=1.5,mai=c(0.7,0.6,0.05,0.3))
# 
# ## plot BLOSUM30
# do.call(vioplot.las2,c(blos30.ordered,col="yellow",list(names=order[,1])))
# mtext("position",side=1, line=2.5, cex=1.5)
# mtext("BLOSUM30",side=2, line=2, cex=1.5)
# 
# ## plot BLOSUM30-NS
# do.call(vioplot.las2,c(blos30.ns.ordered,col="orange",list(names=order[,1])))
# mtext("position",side=1, line=2.5, cex=1.5)
# mtext("BLOSUM30-NS",side=2, line=2, cex=1.5)
# 
# ## plot BLOSUM62
# do.call(vioplot.las2,c(blos62.ordered,col="yellow",list(names=order[,1])))
# mtext("position",side=1, line=2.5, cex=1.5)
# mtext("BLOSUM62",side=2, line=2, cex=1.5)
# 
# ## plot BLOSUM62-NS
# do.call(vioplot.las2,c(blos62.ns.ordered,col="orange",list(names=order[,1])))
# mtext("position",side=1, line=2.5, cex=1.5)
# mtext("BLOSUM62-NS",side=2, line=2, cex=1.5)

## boxbubbleplot
# par(mfrow=c(2,2),tcl=-0.5)
# boxbubbleplot('BLOSUM30',resNumCol,77,lastCol,data,1.5)
# boxbubbleplot('BLOSUM30-NS',resNumCol,77,lastCol,data.ns.only,1.5)
# boxbubbleplot('BLOSUM62',resNumCol,66,lastCol,data,1.5)
# boxbubbleplot('BLOSUM62-NS',resNumCol,66,lastCol,data.ns.only,1.5)



## vioplots #########
#### plot ordered blosum30, 40, 50, 62 all and NS-only
# blos40.ordered = array(list(NULL), c(maxResNum,1))
# blos50.ordered = array(list(NULL), c(maxResNum,1))
# blos40.ns.ordered = array(list(NULL), c(maxResNum,1))
# blos50.ns.ordered = array(list(NULL), c(maxResNum,1))
# ctr = 1
# for(i in order[,1])
# {
#   blos40.ordered[[ctr]] = data$V79[data[,lastcol] == i] 
#   blos50.ordered[[ctr]] = data$V81[data[,lastcol] == i]
#   blos40.ns.ordered[[ctr]] = data.ns.only$V79[data.ns.only[,lastcol] == i] 
#   blos50.ns.ordered[[ctr]] = data.ns.only$V81[data.ns.only[,lastcol] == i] 
#   ctr = ctr + 1
# }
# 
# x11()
# par(mfrow=c(4,2),tcl=-0.5)
# par(cex.axis=1.5,cex.lab=1.5,mai=c(0.7,0.5,0,0.3))

# ## plot BLOSUM30
# do.call(vioplot.las2,c(blos30.ordered,col="yellow",list(names=order[,1])))
# mtext("position",side=1, line=2.5, cex=1.2)
# mtext("BLOSUM30",side=2, line=2, cex=1.2)
# 
# ## plot BLOSUM30-NS
# do.call(vioplot.las2,c(blos30.ns.ordered,col="orange",list(names=order[,1])))
# mtext("position",side=1, line=2.5, cex=1.2)
# mtext("BLOSUM30-NS",side=2, line=2, cex=1.2)
# 
# ## plot BLOSUM40
# do.call(vioplot.las2,c(blos40.ordered,col="yellow",list(names=order[,1])))
# mtext("position",side=1, line=2.5, cex=1.2)
# mtext("BLOSUM40",side=2, line=2, cex=1.2)
# 
# ## plot BLOSUM40-NS
# do.call(vioplot.las2,c(blos40.ns.ordered,col="orange",list(names=order[,1])))
# mtext("position",side=1, line=2.5, cex=1.2)
# mtext("BLOSUM40-NS",side=2, line=2, cex=1.2)
# 
# ## plot BLOSUM50
# do.call(vioplot.las2,c(blos50.ordered,col="yellow",list(names=order[,1])))
# mtext("position",side=1, line=2.5, cex=1.2)
# mtext("BLOSUM50",side=2, line=2, cex=1.2)
# 
# ## plot BLOSUM50-NS
# do.call(vioplot.las2,c(blos50.ns.ordered,col="orange",list(names=order[,1])))
# mtext("position",side=1, line=2.5, cex=1.2)
# mtext("BLOSUM50-NS",side=2, line=2, cex=1.2)
# 
# ## plot BLOSUM62
# do.call(vioplot.las2,c(blos62.ordered,col="yellow",list(names=order[,1])))
# mtext("position",side=1, line=2.5, cex=1.2)
# mtext("BLOSUM62",side=2, line=2, cex=1.2)
# 
# ## plot BLOSUM62-NS
# do.call(vioplot.las2,c(blos62.ns.ordered,col="orange",list(names=order[,1])))
# mtext("position",side=1, line=2.5, cex=1.2)
# mtext("BLOSUM62-NS",side=2, line=2, cex=1.2)

### boxbubbleplot
# boxbubbleplot('BLOSUM30',resNumCol,77,lastCol,data,1.5)
# boxbubbleplot('BLOSUM30-NS',resNumCol,77,lastCol,data.ns.only,1.5)
# boxbubbleplot('BLOSUM40',resNumCol,79,lastCol,data,1.5)
# boxbubbleplot('BLOSUM40-NS',resNumCol,79,lastCol,data.ns.only,1.5)
# boxbubbleplot('BLOSUM50',resNumCol,81,lastCol,data,1.5)
# boxbubbleplot('BLOSUM50-NS',resNumCol,81,lastCol,data.ns.only,1.5)
# boxbubbleplot('BLOSUM62',resNumCol,66,lastCol,data,1.5)
# boxbubbleplot('BLOSUM62-NS',resNumCol,66,lastCol,data.ns.only,1.5)