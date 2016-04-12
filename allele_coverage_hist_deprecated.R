library(graphics)

setwd('C:/Users/JM/thesis/mark_work/allele_specificity/datasets')
data = read.table('380_intHets_counts_gCov.txt',header=T,stringsAsFactors = F)
ordered = order(data$rnaseq.intHets)

## reorder sample as factors/levels (for plotting)
x=factor(data$sample)
x=factor(x,levels(x)[ordered])

## gCoverage histogram plot1
x11()
hist(data$gCov, 380, main=paste("Coverage for ",nrow(data),"individuals"),
              xlab='average coverage', col="black")
# axis(side=1, at=seq(0,25,1))

# boxplot(data$gCov)

##rnaseq intHets histogram plot2
x11()
hist(data$rnaseq.intHets, 380, main=paste("ACounts for ",nrow(data),"individuals"),
     xlab='counts', col="black")
# axis(side=1, at=seq(0,25,1))

## barplot for intHets and dotplot for gCov plot3
## right axis faulty
x11()
ordered.labels.1 = x
ordered.labels.2 = round(seq(0,max(data$gCov)))
ordered.labels.3 = seq(0,max(data$rnaseq.intHets),by=1000)

# barplot
par(mar=c(5,4,4,4),xpd=TRUE)
m = numeric(nrow(data))
# color by pop
for(i in 1:nrow(data))
{
  if(data$pop[i] == "GBR"){m[i] = 1} else
    if(data$pop[i] == "FIN"){m[i] = 2} else
      if(data$pop[i] == "CEU"){m[i] = 3} else
        if(data$pop[i] == "YRI"){m[i] = 4} else
          if(data$pop[i] == "CHB"){m[i] = 5} else
            if(data$pop[i] == "JPT"){m[i] = 6} else
              if(data$pop[i] == "TSI"){m[i] = 7}
}
data = cbind(data,m)
# GBR FIN CEU YRI CHB JPT TSI
cols = c('red','purple','blue','green',
           'orange','yellow','grey')[data$m[ordered]]
ctplot=barplot(data$rnaseq.intHets[ordered], col=cols,border=cols, yaxt="n")
axis(4, at=seq(0,max(data$rnaseq.intHets),by=1000),
     labels=ordered.labels.3, pos=458) ## put the axis right after the plot

# dot plot
par(new=TRUE) ## new plot
plot(x[ordered],data$gCov[ordered],axes=FALSE,bty = "n")
axis(2, at=seq(0,max(data$gCov)),
     labels=ordered.labels.2, pos=c(-2,-5)) 
axis(1, at=x, labels=ordered.labels.1)
mtext("coverage",side=2,line=2)
mtext("intHets.rnaseq",side=4,line=1.5)
mtext("sample",side=1,line=2)

# CEU CHB FIN GBR JPT TSI YRI run this only after the plotting
legend(200, 25, sort(unique(factor(data$pop))), 
       col=c('blue','orange','purple','red','yellow','grey','green'), 
       text.col = "black", pch = 20, bg = 'white', horiz = 1)