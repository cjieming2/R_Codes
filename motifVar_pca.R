setwd('C:/Users/JM/thesis/lynne_work/motifVar_pNets/tpr/')

## relative entropy data
sorted.rel.ent = read.table('1KG.snps.nonmono.smartDomain2gPos.TPR.34aa.sorted.enrich2', row.names=1, header=T, stringsAsFactors = F)

## pca
pca.data = sorted.rel.ent[,-1:-2]
pca.fit = princomp(pca.data, scale=TRUE, center=TRUE)
pca.sum = summary(pca.fit) ## variance accounted for
pca.load= loadings(pca.fit) ## pc loadings
x11(); plot(pca.fit) ## scree plot
#x11(); biplot(pca.fit, xlabs=row.names(pca.data))
x11(); biplot(pca.fit, xlabs=row.names(pca.data), choices=c(1,2))
pca.fit$scores # pcs

## plot pc 1 and 2 for samples PC scores
colors <- c("green", "red", "blue","purple","magenta","pink","cyan") ##Set of colours to use for all plots
categories = sorted.rel.ent[,1]
x11()
xaxis <- pca.fit$scores[,1]
yaxis <- pca.fit$scores[,2]
plot(xaxis,yaxis,main="Plot of PCA scores", bg = colors[unclass(factor(categories))],
     pch=21, xlab="PC1",ylab="PC2", cex = 3)
text(xaxis-0.48,yaxis+0.48,labels=row.names(sorted.rel.ent))
legend(locator(1), sort(matrix((unique(categories)))), 
       pt.bg = colors,text.col = "black", pch = 21, bg = 'white')
##identify(xaxis,yaxis,labels=as.character(row.names(sorted.rel.ent)))

# find rship between variables loadings
x11()
xaxis <- pca.fit$loadings[,1]
yaxis <- pca.fit$loadings[,2]
plot(xaxis,yaxis,main="Plot of PCA scores", bg = colors[unclass(factor(categories))],
     pch=21, xlab="PC1",ylab="PC2")
text(xaxis,yaxis,labels=row.names(t(pca.data)))
legend(locator(1), sort(matrix((unique(categories)))), 
       pt.bg = colors,text.col = "black", pch = 21, bg = 'white')
##identify(xaxis,yaxis,labels=as.character(row.names(sorted.rel.ent)))

## pairwise
pairs(pca.fit$scores, bg = colors[unclass(factor(categories))],pch=21,)

## to check the correlatiom matrix
cov.wt(pca.data)

##################prcomp transpose -- correlate the residues -- plot loadings
## pca
pca.data = t(sorted.rel.ent[,-1:-2])
pca.fit = prcomp(pca.data, scale=TRUE, center=TRUE)
pca.sum = summary(pca.fit) ## variance accounted for
pca.load= loadings(pca.fit) ## pc loadings
x11(); plot(pca.fit) ## scree plot
#x11(); biplot(pca.fit, xlabs=row.names(pca.data))
x11(); biplot(pca.fit, xlabs=row.names(pca.data), choices=c(1,2))
pca.fit$scores # pcs

## plot pc 1 and 2 loadings rship between samples
colors <- c("green", "red", "blue","purple","magenta","pink","cyan") ##Set of colours to use for all plots
categories = sorted.rel.ent[,1]
x11()
xaxis <- pca.fit$rotation[,1] ## loadings
yaxis <- pca.fit$rotation[,2]
plot(xaxis,yaxis,main="Plot of PCA scores", bg = colors[unclass(factor(categories))],
     pch=21, xlab="PC1",ylab="PC2")
text(xaxis,yaxis,labels=row.names(sorted.rel.ent))
legend(locator(1), sort(matrix((unique(categories)))), 
       pt.bg = colors,text.col = "black", pch = 21, bg = 'white')
##identify(xaxis,yaxis,labels=as.character(row.names(sorted.rel.ent)))

# plot scores rship between variables
x11()
xaxis <- pca.fit$x[,1] ## scores
yaxis <- pca.fit$x[,2]
plot(xaxis,yaxis,main="Plot of PCA scores", bg = colors[unclass(factor(categories))],
     pch=21, xlab="PC1",ylab="PC2")
text(xaxis,yaxis,labels=row.names(sorted.rel.ent))
legend(locator(1), sort(matrix((unique(categories)))), 
       pt.bg = colors,text.col = "black", pch = 21, bg = 'white')


## pairwise
pairs(pca.fit$rotation, bg = colors[unclass(factor(categories))],pch=21,)

## to check the correlatiom matrix
cov.wt(pca.data)