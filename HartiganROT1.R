

#-----------------------------------------------------
# Function: "harti.rule"
# This function compute the Hartigan's rule of thumb 
# in K-means analysis (clustering analysis).
# Also makes scree-plot like plot of SSquares Within for first 20 clusters.
# The syntax is:>  harti.rule(dataframe)
# It is assumed the first column of the dataframe contains row labels (such as country names)
# Creates Screeplot of SSWithin
# Returns database with cluster membership, SSwithin, Hartigan ROT for various values
# of k in K-means
#
# Date: Revised 2.15.12
#-----------------------------------------------------

harti.rule <-function(complete.dbase){
     dbase=complete.dbase[,-1]
     labs=as.character(complete.dbase[,1])
     num.groups.table=nrow(dbase);kmax=nrow(dbase)-1;
     SSwithin=c();HartiganROT=c();
     tab.matrix <- matrix(nrow=num.groups.table,ncol=kmax);SS.prev=0
     tab.matrix[,1]<- kmeans(dbase,1)$cluster  #clusters by objects
         SSwithin[1]<-sum(kmeans(dbase,1)$withinss)#ths SSwithin clusters
         HartiganROT[1]<-c(0)
         SS.prev <-SSwithin[1]
       for (i in 2:kmax){
         tab.matrix[,i]<- kmeans(dbase,i)$cluster  #clusters by objects
         SSwithin[i]<-sum(kmeans(dbase,i)$withinss)#ths SSwithin clusters
         HartiganROT[i]<-((SS.prev/SSwithin[i])-1)*(kmax-i)     #this is Hartigan rule of thumb
         SS.prev <-SSwithin[i]}
        rownames(tab.matrix)<-rownames(dbase)
        colnames(tab.matrix)<-paste("k=", 1:kmax, sep = "")
        #make screeplot for internal sum of squares
        plot(c(1:20),SSwithin[1:20],type="b",col="red",xlab="K",ylab="Total Within Cluster SS",lwd=2,pch=16,cex=1.3)
        outdata=round(rbind(tab.matrix,SSwithin,HartiganROT),1) #out.table
        rownames(outdata)=c(labs,"SSwithin","HartiganROT")
        outdata
}


###Example
#Poverty Dataset : 80 most populous countries and transformed variables.
#Remove # from next three lines to run example

povclust=read.csv("http://reuningscherer.net/stat660/data/pov2000transredclust.csv",header=T)
povnorm<-data.frame(povclust[,1],scale(povclust[,2:14]))
harti.rule(povnorm)





