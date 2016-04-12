library(reshape)
library(ggplot2)

setwd('C:/Users/JM/thesis/lynne_work/motifVar_pNets/tpr/')

## create lines connecting boxes
lines_df <- structure(
  list(x =    c(1, 1, 2, 2, 2, 3), 
       y =    c(0.815, 0.83, 0.815, 0.92, 0.93, 0.92), 
       xend = c(1, 2, 2, 2, 3, 3), 
       yend = c(0.83, 0.83, 0.83, 0.93, 0.93, 0.93)), 
  .Names = c("x", "y", "xend", "yend"),
  row.names = c(NA, -6L), class = "data.frame"
)

## information about asterisk position

astpos_df1 <- structure(
  list(x = 1.5, y = 0.85), 
  .Names = c("x", "y"), 
  row.names = c(NA, -1L),
  class = "data.frame"
)

astpos_df2 <- structure(
  list(x = 2.5, y = 0.95), 
  .Names = c("x", "y"), 
  row.names = c(NA, -1L),
  class = "data.frame"
)

## relative entropy data
sorted.rel.ent = read.table('1KG.snps.nonmono.smartDomain2gPos.TPR.34aa.sorted.sorting.oldconserved', 
                            header=T, stringsAsFactors = F)

## conserved=1, none=2, hyperv=3
catnum = 0
sorted.rel.ent = cbind(sorted.rel.ent,catnum)
sorted.rel.ent$catnum[sorted.rel.ent$cat == "conserved"] = 1
sorted.rel.ent$catnum[sorted.rel.ent$cat == "none"] = 3
sorted.rel.ent$catnum[sorted.rel.ent$cat == "hypervar"] = 2

cons = subset(sorted.rel.ent, sorted.rel.ent$cat == "conserved")
none = subset(sorted.rel.ent, sorted.rel.ent$cat == "none")
hypv = subset(sorted.rel.ent, sorted.rel.ent$cat == "hypervar")

## boxplot P.NS.noS
# column = "P.NS.noS"
# x = wilcox.test(cons$P.NS.noS, hypv$P.NS.noS)
# y = wilcox.test(cons$P.NS.noS, none$P.NS.noS)
# z = wilcox.test(none$P.NS.noS, hypv$P.NS.noS)
# x11()
# par(mar=c(5,8,4,4),xpd=TRUE)
# par(cex.axis=1.5, cex.lab=1.5)
# boxplot(P.NS.noS~catnum,sorted.rel.ent, 
#         names=c('conser','hyper','none'),
#         col  =c('red','blue','grey'),
#         xlab="category", ylab=column, main=column)
# text(1.5,0.7,round(x$p.value,2), cex=2)
column = "P.NS.noS"
x11()
p <- ggplot(sorted.rel.ent, aes(x=cat,y=P.NS.noS)) + 
          geom_boxplot(outlier.shape=NA, aes(fill=factor(cat))) 

a = wilcox.test(cons$P.NS.noS, hypv$P.NS.noS)
c = wilcox.test(cons$P.NS.noS, none$P.NS.noS)
b = wilcox.test(none$P.NS.noS, hypv$P.NS.noS)

p2 <- p + geom_segment(data = lines_df, size = 1, aes(x=x, y=y, xend=xend, yend=yend)) 
p3 <- p2 + geom_text(data = astpos_df1, aes(x=x, y=y), label=paste("p=",round(a$p.value,2)), size = 8)
p4 <- p3 + geom_text(data = astpos_df2, aes(x=x, y=y), label=paste("p=",round(b$p.value,2)), size = 8)
p4 + geom_point() + guides(fill = FALSE) + 
                    labs(title=column, x="category",y=column) + 
                    theme(text=element_text(size=14), axis.text=element_text(size=14)) + 
                    scale_fill_manual(values=c("red","blue","grey"))

## boxplot P.rare.noS
column = "P.rare.noS"
x11()
p <- ggplot(sorted.rel.ent, aes(x=cat,y=P.rare.noS)) + 
  geom_boxplot(outlier.shape=NA, aes(fill=factor(cat))) 

a = wilcox.test(cons$P.rare.noS, hypv$P.rare.noS)
c = wilcox.test(cons$P.rare.noS, none$P.rare.noS)
b = wilcox.test(none$P.rare.noS, hypv$P.rare.noS)

p2 <- p + geom_segment(data = lines_df, size = 1, aes(x=x, y=y, xend=xend, yend=yend)) 
p3 <- p2 + geom_text(data = astpos_df1, aes(x=x, y=y), label=paste("p=",round(a$p.value,2)), size = 8)
p4 <- p3 + geom_text(data = astpos_df2, aes(x=x, y=y), label=paste("p=",round(b$p.value,2)), size = 8)
p4 + geom_point() + guides(fill = FALSE) + 
  labs(title=column, x="category",y=column) + 
  theme(text=element_text(size=14), axis.text=element_text(size=14)) + 
  scale_fill_manual(values=c("red","blue","grey"))

## boxplot P.rare
column = "P.rare"
x11()
p <- ggplot(sorted.rel.ent, aes(x=cat,y=P.rare)) + 
  geom_boxplot(outlier.shape=NA, aes(fill=factor(cat))) 

a = wilcox.test(cons$P.rare, hypv$P.rare)
c = wilcox.test(cons$P.rare, none$P.rare)
b = wilcox.test(none$P.rare, hypv$P.rare)

p2 <- p + geom_segment(data = lines_df, size = 1, aes(x=x, y=y, xend=xend, yend=yend)) 
p3 <- p2 + geom_text(data = astpos_df1, aes(x=x, y=y), label=paste("p=",round(a$p.value,2)), size = 8)
p4 <- p3 + geom_text(data = astpos_df2, aes(x=x, y=y), label=paste("p=",round(b$p.value,2)), size = 8)
p4 + geom_point() + guides(fill = FALSE) + 
  labs(title=column, x="category",y=column) + 
  theme(text=element_text(size=14), axis.text=element_text(size=14)) + 
  scale_fill_manual(values=c("red","blue","grey"))

## boxplot ratio.comm2rare.noS
column = "ratio.comm2rare.noS"
x11()
p <- ggplot(sorted.rel.ent, aes(x=cat,y=ratio.comm2rare.noS)) + 
  geom_boxplot(outlier.shape=NA, aes(fill=factor(cat))) 

a = wilcox.test(cons$ratio.comm2rare.noS, hypv$ratio.comm2rare.noS)
c = wilcox.test(cons$ratio.comm2rare.noS, none$ratio.comm2rare.noS)
b = wilcox.test(none$ratio.comm2rare.noS, hypv$ratio.comm2rare.noS)

p2 <- p + geom_segment(data = lines_df, size = 1, aes(x=x, y=y, xend=xend, yend=yend)) 
p3 <- p2 + geom_text(data = astpos_df1, aes(x=x, y=y), label=paste("p=",round(a$p.value,2)), size = 8)
p4 <- p3 + geom_text(data = astpos_df2, aes(x=x, y=y), label=paste("p=",round(b$p.value,2)), size = 8)
p4 + geom_point() + guides(fill = FALSE) + 
  labs(title=column, x="category",y=column) + 
  theme(text=element_text(size=14), axis.text=element_text(size=14)) + 
  scale_fill_manual(values=c("red","blue","grey"))

## boxplot ratio.comm2rare
column = "ratio.comm2rare"
x11()
p <- ggplot(sorted.rel.ent, aes(x=cat,y=ratio.comm2rare)) + 
  geom_boxplot(outlier.shape=NA, aes(fill=factor(cat))) 

a = wilcox.test(cons$ratio.comm2rare, hypv$ratio.comm2rare)
c = wilcox.test(cons$ratio.comm2rare, none$ratio.comm2rare)
b = wilcox.test(none$ratio.comm2rare, hypv$ratio.comm2rare)

p2 <- p + geom_segment(data = lines_df, size = 1, aes(x=x, y=y, xend=xend, yend=yend)) 
p3 <- p2 + geom_text(data = astpos_df1, aes(x=x, y=y), label=paste("p=",round(a$p.value,2)), size = 8)
p4 <- p3 + geom_text(data = astpos_df2, aes(x=x, y=y), label=paste("p=",round(b$p.value,2)), size = 8)
p4 + geom_point() + guides(fill = FALSE) + 
  labs(title=column, x="category",y=column) + 
  theme(text=element_text(size=14), axis.text=element_text(size=14)) + 
  scale_fill_manual(values=c("red","blue","grey"))

## boxplot P.RareNS.noS
column = "P.RareNS.noS"
x11()
p <- ggplot(sorted.rel.ent, aes(x=cat,y=P.RareNS.noS)) + 
  geom_boxplot(outlier.shape=NA, aes(fill=factor(cat))) 

a = wilcox.test(cons$P.RareNS.noS, hypv$P.RareNS.noS)
c = wilcox.test(cons$P.RareNS.noS, none$P.RareNS.noS)
b = wilcox.test(none$P.RareNS.noS, hypv$P.RareNS.noS)

p2 <- p + geom_segment(data = lines_df, size = 1, aes(x=x, y=y, xend=xend, yend=yend)) 
p3 <- p2 + geom_text(data = astpos_df1, aes(x=x, y=y), label=paste("p=",round(a$p.value,2)), size = 8)
p4 <- p3 + geom_text(data = astpos_df2, aes(x=x, y=y), label=paste("p=",round(b$p.value,2)), size = 8)
p4 + geom_point() + guides(fill = FALSE) + 
  labs(title=column, x="category",y=column) + 
  theme(text=element_text(size=14), axis.text=element_text(size=14)) + 
  scale_fill_manual(values=c("red","blue","grey"))

## boxplot P.RareNS
column = "P.RareNS"
x11()
p <- ggplot(sorted.rel.ent, aes(x=cat,y=P.RareNS)) + 
  geom_boxplot(outlier.shape=NA, aes(fill=factor(cat))) 

a = wilcox.test(cons$P.RareNS, hypv$P.RareNS)
c = wilcox.test(cons$P.RareNS, none$P.RareNS)
b = wilcox.test(none$P.RareNS, hypv$P.RareNS)

p2 <- p + geom_segment(data = lines_df, size = 1, aes(x=x, y=y, xend=xend, yend=yend)) 
p3 <- p2 + geom_text(data = astpos_df1, aes(x=x, y=y), label=paste("p=",round(a$p.value,2)), size = 8)
p4 <- p3 + geom_text(data = astpos_df2, aes(x=x, y=y), label=paste("p=",round(b$p.value,2)), size = 8)
p4 + geom_point() + guides(fill = FALSE) + 
  labs(title=column, x="category",y=column) + 
  theme(text=element_text(size=14), axis.text=element_text(size=14)) + 
  scale_fill_manual(values=c("red","blue","grey"))