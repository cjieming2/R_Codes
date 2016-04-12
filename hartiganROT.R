# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
harti.rule <-function(dbase,have.rownames){
#-----------------------------------------------------
# Function: "harti.rule"
# This function compute the Hartigan's rule of thumb # in K-means analysis (clustering analysis).
# How this function works?:
# The syntax is:>  harti.rule(data,have.rownames=1) # Where:
# - Data: is the dataframe, which is assumed to have in the first colum # the labels for the object (species, plots, etc) to cluster. Also # the dataframe can already have the labels of the objects to cluster # as rownames # - have.rownames: is an indicator variable: 1 indicates that # the dataframe already have the rownames with the objects labels (i.e # the first colum of the data frame have the first variable to be analyzed), # and a value of 0 indicates that you have the objects' labels in # the first column (if so you also need to load the function "ext.obj.names") # # Note that the Hartigan ratio for k=2 is just set to cero, but it should # not be considered in the analysis.
#
# Author: Christian Salas
# Date: February 20, 2007
# New Haven, CT, USA
#-----------------------------------------------------
 if (have.rownames==0) {dbase=ext.obj.names(dbase)}
     num.groups.table=nrow(dbase);kmax=nrow(dbase)-1;
     SSwithin=c();HartiganROT=c();
     tab.matrix <- matrix(nrow=num.groups.table,ncol=kmax-1);SS.prev=0
       for (i in 2:kmax){
         tab.matrix[,i-1]<- kmeans(dbase,i)$cluster  #clusters by objects
         SSwithin[i-1]<-sum(kmeans(dbase,i)$withinss)#ths SSwithin clusters
         HartiganROT[i-1]<-(SS.prev/SSwithin[i-1]-1)*(kmax-i) #ths Hartigan rule of thumb
         SS.prev <-SSwithin[i-1]}
        rownames(tab.matrix)<-rownames(dbase)
        colnames(tab.matrix)<-paste("k=", 2:kmax, sep = "")
        round(rbind(tab.matrix,SSwithin,HartiganROT),1) #out.table
               }
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
ext.obj.names <-function(dbase1){
#-----------------------------------------------------
# Function: "ext.obj.names"
# This function extract the first column (which has label for # the rows) and assign those labels as rownames for the data frame.
# This function is used for the function "harti.rule" (for conducting # k-means analysis) # # Author: Christian Salas # Date: February 20, 2007 # New Haven, CT, USA
#-----------------------------------------------------
          dbase<-dbase1[,2:ncol(dbase1)]
          rownames(dbase)<-(dbase1[,1])
          dbase}
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


