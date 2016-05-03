## my own library
source("/Users/jiemingchen/R_codes/jmRlib.R")

library(RImmPort)
library(DBI)
library(sqldf)
library(plyr)
library(RMySQL)

## studies i want
study_ids=c('SDY131','SDY132','SDY133','SDY134','SDY289','SDY290','SDY291',
            'SDY292','SDY294','SDY352','SDY354','SDY355','SDY356','SDY357',
            'SDY358','SDY546','SDY567','SDY662','SDY668','SDY671','SDY674')

recip_ids=c('SDY131','SDY132','SDY133','SDY134','SDY352','SDY354','SDY355',
            'SDY356','SDY357','SDY358','SDY546','SDY567','SDY662','SDY668','SDY671','SDY674')

donor_ids=c('SDY289','SDY290','SDY291','SDY292','SDY294','SDY356')

## where the serialized data is locally
studies_dir = "/Users/jiemingchen/Documents/immport/data"

####--------------------------------------------####
#### starting up, download and serialize specific studies - do once ####
# ## set up tunnel
# system("/Users/jiemingchen/.start_tunnel.sh")
# 
# ## a script with db connect 
# source("/Users/jiemingchen/R_codes/immport_connect.R")
# 
# ## serialize specific studies on local directory
# serialzeStudyData(study_ids, data_dir=studies_dir)
# 
# ## close tunnel
# system("/Users/jiemingchen/.stop_tunnel.sh")
####--------------------------------------------####

####---------------####
#### studies stats ####

# list.files(studies_dir)

## refer to Word doc for the following: ####
## recipients v donors
## organs
## Interventional/observational


####--------------####
#### subject info ####

## loop through the studies and compile
## checked that the total is just 1000 - small - safe to compile
#study_id = "SDY356"

## recipients
ids = recip_ids
for (i in 1:length(ids))
{
  ## load serialized study data into R 
  data.study = as.data.frame(loadSerializedStudyData(studies_dir, ids[i], "Demographics")[1])
  
  ## gender and samplesize stats
  if(!exists("data.recip"))
  {
    data.recip = data.study
  } else
  {
    data.recip = rbind(data.recip,data.study)
  }

}

## donors
ids = donor_ids
for (i in 1:length(ids))
{
  ## load serialized study data into R 
  data.study = as.data.frame(loadSerializedStudyData(studies_dir, ids[i], "Demographics")[1])
  
  ## gender and samplesize stats
  if(!exists("data.donor"))
  {
    data.donor = data.study
  } else
  {
    data.donor = rbind(data.donor,data.study)
  }
  
}





