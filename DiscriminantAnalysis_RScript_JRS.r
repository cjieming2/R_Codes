###########################################################
###  Discriminant Analysis in R 
###  Multivariate Statistics
###  J. Reuning-Scherer
###########################################################


###########################################################
##  Land Cover Data - Reich.csv
##  N - concentration of nitrogen
##  AMASS - - mass-based net photosynthetic capacity
##  AAREA- area-based net photosynthetic capacity
##  GS - leaf diffuse conductance at photosynthetic capacity
##  LSLA - log10 transformation of specific leaf area
##  FUNCTION - 1=Forbes, 2=Shrub
##  LOCATION - 1=Colorado, 2=Wisconsin
##  GROUP - a combination of Function and Locatio
##    
###########################################################


#load MASS package
library(MASS)
#read data
forbs=read.csv("http://reuningscherer.net/stat660/data/reich.csv",header=T)


#see if data is multivariate normal in EACH GROUP
#get online function
source("http://www.reuningscherer.net/STAT660/R/CSQPlot.r.txt")

#examine multivariate normality within each group
reschisqplot(forbs[forbs[,2]==1,4:8],label="Forbs")
reschisqplot(forbs[forbs[,2]==2,4:8],label="Shrubs")

#make matrix plot to look at differences between groups
plot(forbs[,4:8],col=forbs[,2]+2,pch=forbs[,2]+15,cex=1.2)

#compare standard deviations in each group
sumstats=round(sqrt(aggregate(forbs[,4:8],by=list(forbs[,2]),FUN=var)),2)[,-1]
rownames(sumstats)=c("Forbes","Shrubs")
sumstats

#it appears that the covariances matrices are NOT the same between groups - however,
#we ignore this go being with and fix later (i.e. use quadratic discrim analysis)

#run linear discriminant analysis
forb.disc=lda(forbs[,4:8],grouping=forbs[,2],prior=c(.5,.5))
summary(forb.disc)
#get univarite and multivariate comparisons
forb.manova=manova(as.matrix(forbs[,4:8])~forbs[,2])
summary.manova(forb.manova,test="Wilks")
summary.aov(forb.manova)




###########################################################
##  Depression Data - Depression.csv.  297 individuals in UCLA study
##  CASES - 1=Clinical Depression, 0=No Clinical Depression
##  Education (1-7 scale with 7 being the most)
##  Income (thousands dollars per year)
##  Health (1-4 scale from Excellent=1 to Poor=4)
##  Age (in years)
##    
###########################################################


#load MASS package
library(MASS)
#read data
depress=read.csv("http://reuningscherer.net/stat660/data/depression.csv",header=T)
#keep subset of columns
depress=depress[,c(3,5,7,32,30)]

#see if data is multivariate normal in EACH GROUP
#get online function
source("http://www.reuningscherer.net/STAT660/R/CSQPlot.r.txt")

#examine multivariate normality within each group
reschisqplot(depress[depress$CASES==1,1:4],label="Depression")
reschisqplot(depress[depress$CASES==0,1:4],label="No Depression")

#make matrix plot to look at differences between groups
plot(depress[,1:4],col=depress[,5]+3,pch=depress[,5]+15,cex=1.2)

#compare standard deviations in each group
sumstats=round(sqrt(aggregate(depress[,1:4],by=list(depress$CASES),FUN=var)),2)[,-1]
rownames(sumstats)=c("No Depression","Depression")
sumstats

#run linear discriminant analysis
depress.disc=lda(depress[,1:4],
   grouping=depress$CASES,prior=c(.5,.5))
summary(depress.disc)
depress.manova=manova(as.matrix(depress[,1:4])~depress$CASES)
summary.manova(depress.manova,test="Wilks")
summary.aov(depress.manova)

